"""主模拟引擎: 13步 Euler-forward 时间步循环

V4 核心修正:
- 氧账本闭合: 冲击区 O2 全额分配 (C直接/Si/Mn/逸散份额, 其余全部 Fe→FeO),
  每步记账可审计, 不再有 39% O2 凭空消失
- Fe 作溶剂 (a_Fe=1): FeO 界面通量由渣侧传质控制 J_FeO = F_FeO·(w_FeO^b − w_FeO_eq(a_O*)),
  w_FeO_eq = 100·K_Fe·a_O*·M_FeO/(γ_FeO·M̄), 取代把 Fe 塞进式(7)的错误约定
- 删除 J_O 项 (原实现用渣 FeO 平衡氧势冒充钢液本体氧, 与 J_Fe 重复计氧且不守恒)
- 热平衡符号化: 氧化 +ΔH / 还原 −ΔH, FeO 分解吸热自动抵扣,
  如 (FeO)+[C]→[Fe]+CO 净热 = 140−240 = −100 kJ/mol (吸热, 符合实际)
- 废钢不再双计: 初始熔池 = 铁水, 熔化后才进入 (w_C(0)=铁水碳含量)
- q_Ar 单位修正 (m³/min → Nm³/s 除以 60, 原除以 60000 错 1000 倍)
- 石灰/白云石按实际分批 (0/3/6 min) 加入, 烧结矿保留 O2 进度窗口 (冷却剂操作)
"""
import numpy as np
from config import (RHO_STEEL, RHO_SLAG, CP_STEEL, CP_SLAG, KITAMURA_BASE,
                    TIME_STEP_INTENSE, TIME_STEP_NORMAL, TIME_STEP_END, SOLVER_TOL,
                    SLAG_REF_MASS, O2_SHARE, SCRAP_CP, SCRAP_LATENT,
                    CALIBRATION_DEFAULTS, O2_RAMP_TIME)
from thermo.constants import (OXIDATION_THERMO, MOLAR_MASS, OXIDE_MOLAR_MASS,
                              OXIDE_UNIT_MOLAR_MASS, HEAT_C3P_FORMATION,
                              HEAT_LIME_DECOMPOSITION, HEAT_DOLOMITE_DECOMPOSITION,
                              CP_LIME_SOLID, LIME_LOI_FRACTION)
from thermo.equilibrium import calc_K_M
from thermo.wagner import calc_all_f_i
from thermo.ion_theory import calc_slag_activity_coeff
from thermo.slag_phase import (calc_CaO_saturation, calc_solid_fraction,
                               calc_gamma_FeO, calc_gamma_MnO, calc_L_P,
                               calc_slag_viscosity)
from model.static_model import StaticModel
from model.mixing_energy import MixingEnergy
from model.modified_mass_transfer import ModifiedMassTransfer
from model.interface import InterfaceConcentration
from model.oxygen_balance import OxygenBalance
from model.root_solver import RootSolver
from model.lime_melting import LimeMelting
from model.scrap_melting import ScrapMelting
from model.mass_balance import MassBalance
from model.heat_balance import HeatBalance

# 石灰/白云石分批加入时刻 s (张强/设计文档 §4.7: 0/3/6 min 各 1/3)
FLUX_BATCH_TIMES = [0.0, 180.0, 360.0]

# 渣金反应元素及其氧化学计量数 n ([M]+n[O]=MOn)
# P 不在耦合残差中: 主吹期界面 a_O* 被碳钉在 ~1e-4, 纯质量作用式下脱磷必然反向,
# 而真实脱磷发生在氧势接近渣侧的乳化液滴上 —— 故 P 用 Kitamura(2009)式(1)速率方程
# + Suito-Inoue 分配比 L_P 单独处理 (其耗氧仅占总量 ~2%, 氧由渣中 FeO 显式扣减)
COUPLED_ELEMENTS = [('C', 1.0), ('Si', 2.0), ('Mn', 1.0),
                    ('V', 1.5), ('Ti', 2.0)]

# 元素→渣中氧化物质量转换 (kg氧化物/kg元素)
ELEM_TO_OXIDE = {
    'Si': ('SiO2', 0.06009 / 0.02809),
    'Mn': ('MnO', 0.07094 / 0.05494),
    'P':  ('P2O5', 0.14194 / 0.03097 / 2.0),
}


def _actual_or(actual, fallback):
    """实测值优先, 缺失(None)时回退到静态模型估算值 (0 是合法实测值, 不回退)"""
    return fallback if actual is None else actual


class Engine:
    def __init__(self, params: dict = None):
        """params: 标定参数覆盖 (键见 config.CALIBRATION_DEFAULTS)。
        一组参数共享于全部炉次 (设备/工艺级常数), 炉次差异由 run(inputs) 的输入承载。"""
        self.p = dict(CALIBRATION_DEFAULTS)
        if params:
            unknown = set(params) - set(CALIBRATION_DEFAULTS)
            if unknown:
                raise ValueError(f"未知标定参数: {unknown}")
            self.p.update(params)
        self.static = StaticModel()
        self.mixing = MixingEnergy()
        self.transfer = ModifiedMassTransfer()
        self.interface = InterfaceConcentration()
        self.oxy = OxygenBalance()
        self.solver = RootSolver(tol=SOLVER_TOL)
        self.lime = LimeMelting()
        self.scrap = ScrapMelting()
        self.mass_bal = MassBalance()
        self.heat_bal = HeatBalance()

    def select_time_step(self, w_C):
        if w_C > 0.8: return TIME_STEP_INTENSE
        elif w_C > 0.3: return TIME_STEP_NORMAL
        else: return TIME_STEP_END

    def run(self, inputs, max_time=1500.0):
        static_result = self.static.calculate(inputs)

        # --- State initialization: 初始熔池 = 铁水 (废钢熔化后逐步进入) ---
        m_hm = inputs['铁水质量'] * 1000.0
        m_scrap_init = inputs['废钢加入量'] * 1000.0
        T = inputs['铁水温度'] + 273.15

        steel_kg = {}
        for elem in ['C', 'Si', 'Mn', 'P', 'S', 'V', 'Ti']:
            steel_kg[elem] = inputs[f'w([{elem}])'] / 100.0 * m_hm
        steel_kg['Fe'] = m_hm - sum(steel_kg.values())
        m_steel = m_hm
        m_scrap_remaining = m_scrap_init
        T_scrap = 298.0   # 未熔废钢温度 (常温入炉, 一阶弛豫升温)

        # 辅料量: 有实测值用实测 (此时碱度/MgO目标列不再被消费), 缺失回退静态模型估算
        lime_total = _actual_or(inputs.get('实际石灰kg'), static_result['石灰加入量'])
        dolomite_total = _actual_or(inputs.get('实际白云石kg'), static_result['白云石加入量'])
        # Initial slag: hot metal carry-over (~0.8% of HM mass)
        # P2O5 取 0.2%: 高炉渣磷几乎全部进入铁水, 渣中 P2O5 实际 <0.5%
        # (旧值 7% 无依据且直接给回磷供料)
        m_slag_initial = m_hm * 0.008
        slag_kg = {
            'SiO2': m_slag_initial * 0.218, 'CaO': m_slag_initial * 0.40,
            'FeO': m_slag_initial * 0.20, 'MnO': m_slag_initial * 0.08,
            'MgO': m_slag_initial * 0.10, 'P2O5': m_slag_initial * 0.002,
        }
        lime_pool = 0.0      # 已入炉未熔解石灰
        dolo_pool = 0.0      # 已入炉未熔解白云石
        batches_added = 0
        sinter_total = _actual_or(inputs.get('实际烧结矿kg'), inputs['铁水质量'] * 18.0)
        sinter_remaining = sinter_total

        ts = {'time_s': [], 'T_C': [], 'T_K': [],
              'w_C': [], 'w_Si': [], 'w_Mn': [], 'w_P': [], 'w_S': [], 'w_Fe': [],
              'w_FeO': [], 'w_SiO2': [], 'w_CaO': [], 'w_MgO': [], 'w_MnO': [], 'w_P2O5': [],
              'm_steel_kg': [], 'm_slag_kg': [], 'a_O_star': [],
              # FeO 分通道 (kg FeO/步, 佐证前中末期 生成 vs 消耗 谁占上风):
              'feo_gen_impact': [],   # 冲击区 Fe+O2→FeO 生成 (放热, 前/末期主导)
              'feo_net_iface': [],    # 界面耦合净 (+生成/−C还原消耗, 中期消耗为主)
              'feo_cons_deP': []}     # 脱磷 2[P]+5(FeO) 消耗

        q_CO_prev = 0.0; O2_consumed = 0.0; t = 0.0
        self._sinter_deferred = 0.0
        self._a_O_bulk_prev = 1e-4   # 上一步碳氧缓冲本体氧 (冲击区 Mn 门控用, 前向欧拉滞后)
        # 氧量目标: 有实测总耗氧量用实测 (终止条件/加料窗口/废钢周期都以真实氧账为锚)
        O2_total_target = _actual_or(inputs.get('实际耗氧量Nm3'), static_result['总耗氧量'])
        # 终止模式:
        #   氧锚定 (默认, 提供实际耗氧量时): 吹完真实氧量读终点 —— 最严苛的检验框架
        #   目标碳 (inputs['终止模式']='目标碳', 或未提供实际耗氧量): 吹到 C≤目标停,
        #     与张强同口径 —— 模型的使用框架 (过程展示+机理终点); 此时实际耗氧量
        #     (若有) 仍用于加料窗口/废钢周期的氧账锚定, 只是不作为终止条件
        oxygen_anchored = (inputs.get('实际耗氧量Nm3') is not None
                           and inputs.get('终止模式') != '目标碳')
        # 氧账本 (mol O 原子): 每步审计 O_direct + O_inj <= O_blown·(1-escape)+微量
        ledger = {'O_blown': 0.0, 'O_direct': 0.0, 'O_interface_inj': 0.0, 'O_escape': 0.0}
        # 热账本 (J, 全炉累计): 逐项审计热收支, 定位随氧缩放的敏感项
        hledger = {'Q_oxid': 0.0, 'Q_FeO_iface': 0.0, 'Q_deP': 0.0, 'Q_postcomb': 0.0,
                   'Q_C2S': 0.0, 'Q_scrap': 0.0, 'Q_flux_decomp': 0.0, 'Q_slag_heat': 0.0,
                   'Q_offgas': 0.0, 'Q_rad': 0.0, 'Q_sinter': 0.0}

        while t < max_time:
            w_C_pct = steel_kg.get('C', 0) / m_steel * 100 if m_steel > 0 else 0
            dt = self.select_time_step(w_C_pct)
            # 点火爬坡: 开吹 O2 流量从 0 线性升至额定 (真实点火期 ~60s)
            flow_factor = min(1.0, (t + 0.5 * dt) / O2_RAMP_TIME) if O2_RAMP_TIME > 0 else 1.0
            flow_eff = inputs['供氧流量'] * flow_factor
            if oxygen_anchored:
                # 末步裁剪: 精确吹到实际耗氧量为止
                remaining_O2 = O2_total_target - O2_consumed
                dt = min(dt, remaining_O2 / (max(flow_eff, 1e-6) / 60.0))
                if dt <= 1e-9:
                    break

            m_slag = max(sum(slag_kg.values()), 1.0)
            slag_wt = {k: v / m_slag * 100 for k, v in slag_kg.items()}
            f_solid = calc_solid_fraction(slag_wt.get('CaO', 45), slag_wt.get('SiO2', 15), slag_wt.get('FeO', 20))
            steel_wt_pct = {e: steel_kg.get(e, 0) / m_steel * 100 if m_steel > 0 else 0
                            for e in ['C', 'Si', 'Mn', 'P', 'S', 'V', 'Ti', 'Fe']}

            # === Steps 1-2: Mixing energy + dynamic mass transfer ===
            eps_T = self.mixing.calc_epsilon_T(flow_eff, m_steel, inputs['枪位高度'],
                                               inputs['喷嘴直径'], 300.0, inputs['喷嘴夹角'])
            eps_B = self.mixing.calc_epsilon_B(inputs['底吹氩气流量'], T, m_steel, inputs['熔池深度'])
            eps_total = eps_T + eps_B
            if t == 0.0:
                # 归一基准用额定流量 (爬坡期流量近零, 若以其为基准则后续比值爆炸)
                eps_T_full = self.mixing.calc_epsilon_T(inputs['供氧流量'], m_steel, inputs['枪位高度'],
                                                        inputs['喷嘴直径'], 300.0, inputs['喷嘴夹角'])
                self._eps_total0 = max(eps_T_full + eps_B, 0.01)

            km = self.transfer.calc_k_m(KITAMURA_BASE['k_m0'] * self.p['alpha_steel'],
                                        eps_total, self._eps_total0, gamma=self.p['gamma'])
            q_Ar = inputs['底吹氩气流量'] / 60.0   # m³/min → Nm³/s
            stir = min(1.0 + 3.0 * q_CO_prev / max(q_Ar, 1e-10), 10.0)
            ks_prime = KITAMURA_BASE['k_s0'] * self.p['alpha_slag'] * stir * (SLAG_REF_MASS / m_slag)
            ks = self.transfer.calc_k_s(ks_prime, f_solid, KITAMURA_BASE['n_solid'])
            F_dict = self.transfer.calc_all_F(km, ks, RHO_STEEL, RHO_SLAG, MOLAR_MASS, OXIDE_MOLAR_MASS)

            # === Step 3: Thermodynamics (K_M, 活度, E_M) ===
            K_M = {elem: calc_K_M(elem, T) for elem in OXIDATION_THERMO}
            f_steel = calc_all_f_i(steel_wt_pct, T)
            basicity = slag_wt.get('CaO', 45) / max(slag_wt.get('SiO2', 15), 0.1)
            gamma_FeO = calc_gamma_FeO(basicity)
            gamma_MnO = calc_gamma_MnO(basicity)
            slag_pct = {k: slag_wt[k] for k in ['FeO', 'SiO2', 'CaO', 'MnO', 'MgO', 'P2O5'] if k in slag_wt}
            f_slag = {'FeO': gamma_FeO, 'MnO': gamma_MnO}
            for ox in ['SiO2', 'P2O5']:
                try: f_slag[ox] = calc_slag_activity_coeff(ox, slag_pct, T)
                except Exception: f_slag[ox] = 1.0

            # 炉渣实时平均摩尔质量 M̄ 与体积摩尔数 C (取代固定 60 g/mol)
            inv_M = sum((slag_wt.get(ox, 0.0) / 100.0) / OXIDE_MOLAR_MASS[ox]
                        for ox in slag_kg if ox in OXIDE_MOLAR_MASS)
            M_bar = 1.0 / max(inv_M, 1e-6)
            C_slag = RHO_SLAG / M_bar

            E_M = {}
            ox_map = {'Si': 'SiO2', 'Mn': 'MnO', 'P': 'P2O5', 'V': 'V2O3', 'Ti': 'TiO2'}
            for elem, ox in ox_map.items():
                fM = f_steel.get(elem, 1.0)
                fMOn = f_slag.get(ox, 1.0)
                N_unit = OXIDE_UNIT_MOLAR_MASS[ox]   # 每 mol 元素的氧化物单元 (PO2.5, VO1.5)
                E_M[elem] = 100.0 * C_slag * N_unit * fM * K_M[elem] / (RHO_SLAG * max(fMOn, 1e-6))
            K_C = K_M['C']
            f_C = f_steel.get('C', 1.0)
            K_Fe = K_M['Fe']

            O2_mol_per_s = flow_eff / (60.0 * 0.0224)   # mol O2/s (含点火爬坡)
            O_total_mol = O2_mol_per_s * 2.0 * dt                  # mol O 原子/步
            ledger['O_blown'] += O_total_mol
            ledger['O_escape'] += O_total_mol * O2_SHARE['escape']
            # 末期"入铁氧保留率" retain (A方向): 高碳时=1 (剩余氧全部生成FeO),
            # 碳降到 c_ref=1% 以下线性降到 eta_O2_late。物理: 碳降低→冲击区碳饥饿→
            # 部分过剩氧来不及反应, 以未反应O2/炉膛后燃逸出而非全部生成FeO。
            # 只作用于"入铁"通道 (下方), 碳/Si/Mn 供氧不动 → 既压末期FeO又不伤脱碳。
            retain = self.p['eta_O2_late'] + (1.0 - self.p['eta_O2_late']) * min(w_C_pct / 1.0, 1.0)

            # === Step 4a: 冲击区直接氧化 (氧账本分配, 剩余全部 Fe→FeO) ===
            # 冲击区暴露的主要是 Fe(~95%), 大部分 O2 先生成 FeO 入渣, 再由渣-金界面
            # 耦合反应还原 (间接氧化/传氧机制); 元素选择性由耦合求解器的热力学决定
            direct_oxidized = {}
            si_frac = O2_SHARE['Si'] if steel_kg.get('Si', 0) > 0.01 else 0.0
            # 冲击区 Mn 门控 (13.22 修复之二): 冲击区直接氧化 Mn 必须尊重渣金平衡——
            # 熔池 Mn 已低于与当前渣 (MnO)/氧势 的平衡值时, 烧掉的 Mn 会被渣立刻还回,
            # 净通道关闭。旧实现无热力学背压, 末期把 Mn 生吃到 0 (0.078 实测 vs 0.002)。
            w_Mn_eq_gate = slag_wt.get('MnO', 0.0) / max(E_M['Mn'] * self._a_O_bulk_prev, 1e-9)
            mn_frac = O2_SHARE['Mn'] if (steel_kg.get('Mn', 0) > 0.01
                                         and steel_wt_pct.get('Mn', 0) > w_Mn_eq_gate) else 0.0
            c_frac = self.p['C_direct']
            unused_O = 0.0

            o2_Si = O_total_mol * si_frac
            si_mol = 0.0
            if steel_kg.get('Si', 0) > 0.001:
                si_mol = min(steel_kg['Si'] / MOLAR_MASS['Si'], o2_Si / 2.0)
                dm = si_mol * MOLAR_MASS['Si']
                steel_kg['Si'] -= dm
                slag_kg['SiO2'] += si_mol * OXIDE_MOLAR_MASS['SiO2']
                direct_oxidized['Si'] = dm
            unused_O += o2_Si - 2.0 * si_mol

            o2_Mn = O_total_mol * mn_frac
            mn_mol = 0.0
            if steel_kg.get('Mn', 0) > 0.001 and mn_frac > 0:
                # 单步也不烧穿平衡地板
                mn_room = max(0.0, (steel_wt_pct.get('Mn', 0) - w_Mn_eq_gate)) / 100.0 * m_steel / MOLAR_MASS['Mn']
                mn_mol = min(steel_kg['Mn'] / MOLAR_MASS['Mn'], o2_Mn, mn_room)
                dm = mn_mol * MOLAR_MASS['Mn']
                steel_kg['Mn'] -= dm
                slag_kg['MnO'] += mn_mol * OXIDE_MOLAR_MASS['MnO']
                direct_oxidized['Mn'] = dm
            unused_O += o2_Mn - mn_mol

            o2_C = O_total_mol * c_frac
            c_mol = 0.0
            if steel_kg.get('C', 0) > 0.001:
                c_mol = min(steel_kg['C'] / MOLAR_MASS['C'], o2_C)
                dm = c_mol * MOLAR_MASS['C']
                steel_kg['C'] -= dm
                direct_oxidized['C'] = dm
            unused_O += o2_C - c_mol

            # Fe: 除逸散外的全部剩余 O2 (含 Si/Mn/C 未用完的份额)。
            # retain 决定这部分"入铁氧"中真正生成 FeO 的比例, 其余 (1-retain) 逸出:
            # 高碳时 retain=1 全部成 FeO; 低碳时 retain→eta_O2_late, 过剩氧逸出封顶末期 FeO。
            o2_Fe_avail = O_total_mol * (1.0 - O2_SHARE['escape'] - c_frac - si_frac - mn_frac) + unused_O
            o2_Fe = o2_Fe_avail * retain
            o2_Fe_escaped = o2_Fe_avail * (1.0 - retain)
            fe_mol = 0.0
            if steel_kg.get('Fe', 0) > 0.001:
                fe_mol = min(steel_kg['Fe'] / MOLAR_MASS['Fe'], o2_Fe)
                dm = fe_mol * MOLAR_MASS['Fe']
                steel_kg['Fe'] -= dm
                slag_kg['FeO'] += fe_mol * OXIDE_MOLAR_MASS['FeO']
                direct_oxidized['Fe'] = dm
            feo_gen_impact_step = fe_mol * OXIDE_MOLAR_MASS['FeO']   # 冲击区FeO生成(埋点)
            o2_Fe_escaped += (o2_Fe - fe_mol)   # 钢中Fe不足时的余量也计逸出 (极少发生)
            ledger['O_direct'] += 2.0 * si_mol + mn_mol + c_mol + fe_mol
            ledger['O_escape'] += o2_Fe_escaped

            m_steel = sum(steel_kg.values())
            steel_wt_pct = {e: steel_kg.get(e, 0) / m_steel * 100 if m_steel > 0 else 0
                            for e in ['C', 'Si', 'Mn', 'P', 'S', 'V', 'Ti', 'Fe']}
            m_slag = max(sum(slag_kg.values()), 1.0)
            slag_wt = {k: v / m_slag * 100 for k, v in slag_kg.items()}

            # === Step 4b: 渣-金界面耦合反应 (C,Si,Mn,P,V,Ti 竞争氧; FeO 为供氧方) ===
            steel_b = {f'w([{e}])b': steel_wt_pct.get(e, 0) for e in ['C', 'Si', 'Mn', 'P', 'V', 'Ti']}
            slag_b = {f'w(({o}))b': slag_wt.get(o, 0) for o in ['SiO2', 'MnO', 'P2O5', 'V2O3', 'TiO2']}

            G_CO = max(self.p['G_CO0'] * stir, 1e-6)
            # 反应面积 = 平静渣金界面 + 乳化液滴面积 (脱碳/FeO还原/脱磷主要在液滴上)
            A = (KITAMURA_BASE['A0'] * 75.0 + self.p['A_emulsion']) * (1.0 - f_solid) ** (2.0 / 3.0)
            F_FeO = F_dict.get('F_FeO', 1e-6)
            w_FeO_b = slag_wt.get('FeO', 0.0)

            def w_FeO_eq(a_O_star):
                """与界面氧活度平衡的 FeO 质量分数: a_FeO = K_Fe·a_O* (a_Fe=1)"""
                return 100.0 * K_Fe * a_O_star * OXIDE_MOLAR_MASS['FeO'] / (gamma_FeO * M_bar)

            # 冲击区溶解 [O] 迁移至渣-金界面的直接供氧 (≤3%, 计在 escape 余量内, 标定参数)
            o2_interface_frac = 0.03 * (1.0 - min(w_C_pct / 5.0, 0.9))
            J_O2_inj = O2_mol_per_s * 2.0 / max(A, 1e-6) * o2_interface_frac

            # FeO 库存约束进入残差 (V4.2): 供氧通量在求解时就不得超过渣中 FeO 存量,
            # 否则元素氧化量按未钳位的解施加会凭空产氧 (破坏守恒), 且把 FeO 硬打到 0
            J_FeO_cap = slag_kg.get('FeO', 0.0) / OXIDE_MOLAR_MASS['FeO'] / max(A * dt, 1e-9)

            def residual_func(a_O_star):
                """耗氧需求(a_O*) − 供氧能力(a_O*): 需求单调增/供给单调减 → 唯一根"""
                w_star = self.interface.calc_all(F_dict, steel_b, slag_b, E_M,
                                                 a_O_star, G_CO, K_C, f_C)
                demand = 0.0
                for e, n_e in COUPLED_ELEMENTS:
                    F_M = F_dict.get(f'F_{e}', 0.0)
                    J = F_M * (steel_b.get(f'w([{e}])b', 0.0) - w_star.get(f'w([{e}])*', 0.0))
                    demand += n_e * J
                supply = min(F_FeO * (w_FeO_b - w_FeO_eq(a_O_star)), J_FeO_cap) + J_O2_inj
                return demand - supply

            try:
                a_O_star = self.solver.solve(residual_func, bracket=(1e-10, 0.5))
            except Exception:
                a_O_star = 1e-8
            # 碳氧缓冲本体氧 (13.23 修复③-refined): Mn/P 应与熔池"本体"氧平衡, 而本体氧
            # 被 C-O 反应缓冲 (熔池有碳时氧被碳消耗、不能过度偏离碳氧平衡):
            # a_O_bulk = mult/(K_C·f_C·w_C), 约为界面 a_O* 的 1/3。用它算 Mn 平衡地板,
            # 而 **不动全局 a_O***(全局压 a_O* 会掐断 FeO 间接氧化→FeO/温度爆炸, 已验证)。
            # w_C→0 时 a_O_bulk 自然放开。
            a_O_bulk = min(self.p['aO_cap_CO_mult'] / max(K_C * f_C * max(w_C_pct, 1e-4), 1e-9),
                           a_O_star)   # 本体氧不高于界面氧
            self._a_O_bulk_prev = a_O_bulk   # 供下一步冲击区 Mn 门控 (前向欧拉滞后)

            # === Steps 5-8: 反应通量施加 + 物质平衡 ===
            w_star_all = self.interface.calc_all(F_dict, steel_b, slag_b, E_M,
                                                 a_O_star, G_CO, K_C, f_C)
            # Mn 耦合地板 (13.23-refined): Mn 不被耦合反应氧化到低于与渣 (MnO)/本体氧 的
            # 平衡值 —— 用 a_O_bulk (碳氧缓冲) 而非界面 a_O*, 修末期 Mn 过氧化, 不动全局氧平衡。
            # w([Mn])* 抬高 → 氧化驱动力 (w_b−w*) 减小 → Mn 保留。
            w_Mn_floor = slag_wt.get('MnO', 0.0) / max(E_M['Mn'] * a_O_bulk, 1e-9)
            w_star_all['w([Mn])*'] = max(w_star_all.get('w([Mn])*', 0.0), w_Mn_floor)

            delta_kg = {}
            for elem in ['C', 'Si', 'Mn']:
                w_b = steel_wt_pct.get(elem, 0)
                w_s = w_star_all.get(f'w([{elem}])*', 0)
                F_M = F_dict.get(f'F_{elem}', 0)
                d = A * F_M * MOLAR_MASS[elem] * (w_b - w_s) * dt
                d = min(d, steel_kg.get(elem, 0.0))               # 氧化不超过钢中存量
                if elem in ELEM_TO_OXIDE and d < 0:               # 还原不超过渣中存量
                    ox, factor = ELEM_TO_OXIDE[elem]
                    d = max(d, -slag_kg.get(ox, 0.0) / factor)
                delta_kg[elem] = d

            for elem in ['C', 'Si', 'Mn']:
                steel_kg[elem] = max(0.0, steel_kg.get(elem, 0) - delta_kg.get(elem, 0))
            for elem, (ox, factor) in ELEM_TO_OXIDE.items():
                slag_kg[ox] = max(0.0, slag_kg.get(ox, 0) + delta_kg.get(elem, 0) * factor)

            # FeO 界面通量: J_FeO>0 → FeO 分解供氧, Fe 入钢液; J_FeO<0 → 界面生成 FeO
            # 与残差使用完全相同的库存钳位表达式, 保证求解与施加一致 (氧严格守恒)
            J_FeO_val = min(F_FeO * (w_FeO_b - w_FeO_eq(a_O_star)), J_FeO_cap)
            n_FeO = A * J_FeO_val * dt   # mol (signed)
            n_FeO = min(n_FeO, slag_kg.get('FeO', 0.0) / OXIDE_MOLAR_MASS['FeO'])
            n_FeO = max(n_FeO, -steel_kg.get('Fe', 0.0) / MOLAR_MASS['Fe'])
            slag_kg['FeO'] = max(0.0, slag_kg.get('FeO', 0.0) - n_FeO * OXIDE_MOLAR_MASS['FeO'])
            steel_kg['Fe'] = max(0.0, steel_kg.get('Fe', 0.0) + n_FeO * MOLAR_MASS['Fe'])
            m_steel = sum(steel_kg.values())
            feo_net_iface_step = -n_FeO * OXIDE_MOLAR_MASS['FeO']   # 界面净(+生成/−消耗, 埋点)

            # === 脱磷: Kitamura(2009)式(1)速率方程 + Suito-Inoue 分配比 L_P ===
            # -d[%P]/dt = (A·k_eff)·{[%P] − (%P)/L_P};  氧由渣中 FeO 供给:
            # 2[P] + 5(FeO) + 3(CaO) → (3CaO·P2O5) + 5[Fe]
            # 氧势由 L_P 中的 2.5·log(T.Fe) 项承载 (乳化液滴局部氧势接近渣侧的图像),
            # 不受碳钉死的界面 a_O* 限制 —— 这是脱磷能在主吹期发生的机理所在
            m_slag = max(sum(slag_kg.values()), 1.0)
            slag_wt = {k: v / m_slag * 100 for k, v in slag_kg.items()}
            L_P = calc_L_P(T, slag_wt)
            w_P_b = steel_kg.get('P', 0.0) / m_steel * 100 if m_steel > 0 else 0.0
            w_P_slag = slag_wt.get('P2O5', 0.0) * 0.4364        # %P2O5 → 渣中%P
            w_P_eq = w_P_slag / max(L_P, 1e-6)                  # 与当前渣平衡的钢中P
            F_P = F_dict.get('F_P', 1e-9)
            F_P_slag = ks * RHO_SLAG / (100.0 * MOLAR_MASS['P'])
            R_P = 1.0 / max(F_P, 1e-12) + 1.0 / max(L_P * F_P_slag, 1e-12)  # 钢侧+渣侧串联阻力
            J_P = (w_P_b - w_P_eq) / R_P                        # mol P/(m²·s), >0 为脱磷
            dm_P = A * J_P * MOLAR_MASS['P'] * dt               # kg P (signed)
            cap = 0.5 * abs(w_P_b - w_P_eq) / 100.0 * m_steel   # 单步≤差额一半, 防震荡
            dm_P = max(-cap, min(dm_P, cap))
            dm_P = min(dm_P, steel_kg.get('P', 0.0))
            dm_P = max(dm_P, -slag_kg.get('P2O5', 0.0) * 0.4364)
            n_P = dm_P / MOLAR_MASS['P']                        # mol P (signed)
            n_FeO_forP = 2.5 * n_P
            if n_FeO_forP > 0:                                  # FeO 不足则按比例缩减脱磷
                n_FeO_avail = slag_kg.get('FeO', 0.0) / OXIDE_MOLAR_MASS['FeO']
                if n_FeO_forP > n_FeO_avail:
                    scale = n_FeO_avail / max(n_FeO_forP, 1e-12)
                    dm_P *= scale; n_P *= scale; n_FeO_forP = 2.5 * n_P
            steel_kg['P'] = max(0.0, steel_kg.get('P', 0.0) - dm_P)
            slag_kg['P2O5'] = max(0.0, slag_kg.get('P2O5', 0.0) + dm_P * (0.14194 / 0.03097 / 2.0))
            slag_kg['FeO'] = max(0.0, slag_kg.get('FeO', 0.0) - n_FeO_forP * OXIDE_MOLAR_MASS['FeO'])
            steel_kg['Fe'] = max(0.0, steel_kg.get('Fe', 0.0) + n_FeO_forP * MOLAR_MASS['Fe'])
            m_steel = sum(steel_kg.values())
            feo_cons_deP_step = -n_FeO_forP * OXIDE_MOLAR_MASS['FeO']   # 脱磷消耗(埋点)

            ledger['O_interface_inj'] += A * J_O2_inj * dt

            C_mol_coupled = max(delta_kg.get('C', 0), 0.0) / MOLAR_MASS['C']
            C_mol_direct = direct_oxidized.get('C', 0) / MOLAR_MASS['C']
            q_CO_prev = (C_mol_coupled + C_mol_direct) / max(dt, 1e-6) * 0.0224 if dt > 0 else 0

            # === Step 9: 辅料 — 石灰/白云石分批入炉, Kitamura 溶解; 烧结矿按 O2 进度 (冷却剂) ===
            while batches_added < len(FLUX_BATCH_TIMES) and t >= FLUX_BATCH_TIMES[batches_added]:
                lime_pool += lime_total / len(FLUX_BATCH_TIMES)
                dolo_pool += dolomite_total / len(FLUX_BATCH_TIMES)
                batches_added += 1

            # 渣黏度: Riboud(1981)+Roscoe 固相修正 (D修复, 取代自造指数式)
            eta_slag = calc_slag_viscosity(slag_wt, T, f_solid)
            D_CaO = 1e-9 * np.exp(-5000.0 / T)   # CaO扩散系数, 量级与文献(1e-10~1e-9)相符, 待换源
            w_CaO_sat_lim = calc_CaO_saturation(slag_wt.get('FeO', 20), slag_wt.get('SiO2', 15))
            if lime_pool > 0.01 and slag_wt.get('CaO', 45) < w_CaO_sat_lim:
                kp = self.lime.calc_k_prime(RHO_SLAG, eta_slag, D_CaO, d=0.01, U=0.5)
                dc = w_CaO_sat_lim - slag_wt.get('CaO', 45)
                nl = lime_pool / (3300.0 * 4 / 3 * np.pi * 0.005 ** 3)
                Al = nl * 4 * np.pi * 0.005 ** 2
                dl = self.lime.calc_mass_change_rate(nl, Al, kp, RHO_SLAG, 3300.0, dc, x=1.0) * dt
                dl = min(dl, lime_pool)
            else:
                dl = 0.0
            lime_pool -= dl
            slag_kg['CaO'] += dl * 0.90; slag_kg['MgO'] += dl * 0.03; slag_kg['SiO2'] += dl * 0.02

            dd = min(dolo_pool * dt / 120.0, dolo_pool) if dolo_pool > 0.01 else 0.0
            dolo_pool -= dd
            slag_kg['CaO'] += dd * 0.30; slag_kg['MgO'] += dd * 0.21; slag_kg['SiO2'] += dd * 0.02

            O2_progress = O2_consumed / max(O2_total_target, 1.0)
            flux_frac = max(0.0, min(1.0, (O2_progress - 0.10) / 0.80))
            sinter_target = flux_frac * sinter_total
            ds_target = max(0.0, sinter_target - (sinter_total - sinter_remaining))
            dsinter = min(ds_target * dt / 60.0, sinter_remaining) if sinter_remaining > 0.01 and ds_target > 0 else 0.0
            sinter_remaining -= dsinter
            slag_kg['FeO'] += dsinter * (0.09 + 0.64 * 0.9); slag_kg['CaO'] += dsinter * 0.13
            slag_kg['SiO2'] += dsinter * 0.07; slag_kg['MgO'] += dsinter * 0.02; slag_kg['MnO'] += dsinter * 0.015

            # === Step 10: 废钢熔化 (质量: 唯象抛物线; 热: 升温显热弛豫 + 熔化潜热) ===
            # 前期冷却主要来自固态废钢升温显热 (cp·ΔT≈0.9 MJ/kg, 一阶弛豫前置),
            # 熔化潜热 (272 kJ/kg) 跟随熔化质量 —— 取代原 1.3 MJ/kg 集总值,
            # 使 0-1 min 的废钢吸热效应在时间上出现在正确的位置
            t_end_scrap = O2_total_target / max(inputs['供氧流量'], 1.0) * 60.0 * 0.80
            if t < t_end_scrap and m_scrap_remaining > 0.01:
                frac = t / t_end_scrap
                ds = m_scrap_init * 6.0 * frac * (1.0 - frac) / t_end_scrap * dt
                ds = min(ds, m_scrap_remaining)
            else:
                ds = min(m_scrap_remaining, m_scrap_init / t_end_scrap * dt)
            Q_scrap_heating = 0.0
            if m_scrap_remaining > 0.01 and T > T_scrap:
                dT_sc = (T - T_scrap) * min(dt / self.p['scrap_tau_heat'], 1.0)
                Q_scrap_heating = m_scrap_remaining * SCRAP_CP * dT_sc
                T_scrap += dT_sc
            Q_scrap_melt = ds * (SCRAP_CP * max(T - T_scrap, 0.0) + SCRAP_LATENT)
            m_scrap_remaining -= ds
            steel_kg['Fe'] += ds
            m_steel = sum(steel_kg.values())

            # === Steps 11-12: 热量平衡 (符号化: 氧化 +ΔH, 还原 −ΔH) ===
            Q_oxid = 0.0
            for elem, dm in direct_oxidized.items():
                Q_oxid += dm / MOLAR_MASS[elem] * OXIDATION_THERMO[elem]['delta_H']
            for elem in ['C', 'Si', 'Mn']:
                Q_oxid += delta_kg.get(elem, 0.0) / MOLAR_MASS[elem] * OXIDATION_THERMO[elem]['delta_H']
            # FeO 界面分解吸热 (n_FeO>0) / 界面生成放热 (n_FeO<0)
            Q_FeO_iface = -n_FeO * OXIDATION_THERMO['Fe']['delta_H']
            # 脱磷净热: ΔH_P + C3P形成热/2 − 2.5·ΔH_Fe (FeO供氧), 符号随 n_P (回磷时反号)
            Q_deP = n_P * (OXIDATION_THERMO['P']['delta_H'] + 0.5 * HEAT_C3P_FORMATION
                           - 2.5 * OXIDATION_THERMO['Fe']['delta_H'])
            # 二次燃烧 (eta_post·CO→CO2, 其 O2 计入 escape 份额; 标定参数)
            Q_C_ox = (max(delta_kg.get('C', 0), 0.0) + direct_oxidized.get('C', 0.0)) \
                     / MOLAR_MASS['C'] * OXIDATION_THERMO['C']['delta_H']
            Q_postcomb = self.p['eta_post'] * Q_C_ox
            # C2S 形成热 (净氧化的 Si)
            dm_Si_net = max(delta_kg.get('Si', 0), 0.0) + direct_oxidized.get('Si', 0.0)
            Q_C2S = dm_Si_net / MOLAR_MASS['Si'] * 120000.0
            Q_HG = Q_oxid + Q_FeO_iface + Q_deP + Q_postcomb + Q_C2S

            Q_scrap = (Q_scrap_heating + Q_scrap_melt) * self.p['scrap_heat_scale']
            # 石灰为煅烧CaO: 付"298K→熔池温度"显热 + 烧损部分(3.95%)分解热,
            # 不再错收全额分解热 (D修复, 修正前多收约 1.2 MJ/kg ≈ 每炉 5~6 GJ)
            q_lime_per_kg = CP_LIME_SOLID * (T - 298.0) + LIME_LOI_FRACTION * HEAT_LIME_DECOMPOSITION
            Q_CM = dl * q_lime_per_kg + dd * HEAT_DOLOMITE_DECOMPOSITION
            Q_slag_heat = 0.0
            # 白云石分解产物(0.53氧化物)从分解温1100K升到熔池温度; 石灰升温已含在 q_lime_per_kg (298K起)
            if T > 1100.0: Q_slag_heat += dd * 0.53 * CP_SLAG * (T - 1100.0)
            if T > 298.0: Q_slag_heat += dsinter * 1000.0 * (T - 298.0)
            Q_offgas = q_CO_prev / 0.0224 * 33.0 * (T - 298.0) * dt if q_CO_prev > 0 and T > 298 else 0.0
            Q_rad_lining = self.p['Q_rad_W'] * dt
            Q_sinter_now = dsinter * 1800000.0 if w_C_pct < 2.5 else 0.0
            self._sinter_deferred += dsinter * 1800000.0 - Q_sinter_now
            if self._sinter_deferred > 0 and w_C_pct < 1.0:
                release = min(self._sinter_deferred, self._sinter_deferred * dt / 90.0)
                Q_sinter_now += release; self._sinter_deferred -= release

            Q_consumed = Q_scrap + Q_CM + Q_slag_heat + Q_offgas + Q_rad_lining + Q_sinter_now
            # 热账本累计
            hledger['Q_oxid'] += Q_oxid; hledger['Q_FeO_iface'] += Q_FeO_iface
            hledger['Q_deP'] += Q_deP; hledger['Q_postcomb'] += Q_postcomb; hledger['Q_C2S'] += Q_C2S
            hledger['Q_scrap'] += Q_scrap; hledger['Q_flux_decomp'] += Q_CM
            hledger['Q_slag_heat'] += Q_slag_heat; hledger['Q_offgas'] += Q_offgas
            hledger['Q_rad'] += Q_rad_lining; hledger['Q_sinter'] += Q_sinter_now
            delta_T = self.heat_bal.calc_delta_T(Q_HG, Q_consumed, 0, CP_STEEL, m_steel, CP_SLAG, m_slag)
            T += max(-50.0, min(delta_T, 100.0))

            # === Record ===
            m_slag = sum(slag_kg.values())
            ts['time_s'].append(t); ts['T_C'].append(T - 273.15); ts['T_K'].append(T)
            ts['w_C'].append(steel_kg.get('C', 0) / m_steel * 100 if m_steel > 0 else 0)
            ts['w_Si'].append(steel_kg.get('Si', 0) / m_steel * 100 if m_steel > 0 else 0)
            ts['w_Mn'].append(steel_kg.get('Mn', 0) / m_steel * 100 if m_steel > 0 else 0)
            ts['w_P'].append(steel_kg.get('P', 0) / m_steel * 100 if m_steel > 0 else 0)
            ts['w_S'].append(steel_kg.get('S', 0) / m_steel * 100 if m_steel > 0 else 0)
            ts['w_Fe'].append(steel_kg.get('Fe', 0) / m_steel * 100 if m_steel > 0 else 0)
            ts['w_FeO'].append(slag_kg.get('FeO', 0) / m_slag * 100 if m_slag > 0 else 0)
            ts['w_SiO2'].append(slag_kg.get('SiO2', 0) / m_slag * 100 if m_slag > 0 else 0)
            ts['w_CaO'].append(slag_kg.get('CaO', 0) / m_slag * 100 if m_slag > 0 else 0)
            ts['w_MgO'].append(slag_kg.get('MgO', 0) / m_slag * 100 if m_slag > 0 else 0)
            ts['w_MnO'].append(slag_kg.get('MnO', 0) / m_slag * 100 if m_slag > 0 else 0)
            ts['w_P2O5'].append(slag_kg.get('P2O5', 0) / m_slag * 100 if m_slag > 0 else 0)
            ts['m_steel_kg'].append(m_steel); ts['m_slag_kg'].append(m_slag)
            ts['a_O_star'].append(a_O_star)
            ts['feo_gen_impact'].append(feo_gen_impact_step)
            ts['feo_net_iface'].append(feo_net_iface_step)
            ts['feo_cons_deP'].append(feo_cons_deP_step)

            t += dt
            O2_consumed += flow_eff * dt / 60.0
            w_C_pct = steel_kg.get('C', 0) / m_steel * 100 if m_steel > 0 else 0

            if oxygen_anchored:
                if O2_consumed >= O2_total_target - 1e-6: break
            else:
                if O2_consumed >= O2_total_target * 0.5 and w_C_pct <= inputs['终点碳目标']: break

        O_used = ledger['O_direct'] + ledger['O_interface_inj']
        return {
            '终点C': round(ts['w_C'][-1], 4) if ts['w_C'] else 0.0,
            '终点T': round(ts['T_C'][-1], 1) if ts['T_C'] else 0.0,
            '时序数据': ts,
            '步数': len(ts['time_s']),
            '氧账本': {
                'O_blown_mol': round(ledger['O_blown'], 0),
                'O_direct_mol': round(ledger['O_direct'], 0),
                'O_interface_inj_mol': round(ledger['O_interface_inj'], 0),
                'O_escape_mol': round(ledger['O_escape'], 0),
                '化学利用率': round(O_used / max(ledger['O_blown'], 1e-9), 4),
            },
            '热账本_MJ': {k: round(v / 1e6, 1) for k, v in hledger.items()},
        }
