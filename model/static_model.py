"""转炉静态模型：物料平衡 + 热量平衡反算合理废钢量

复现张强(2024)两段式架构的静态段：给定始态(铁水条件)与终点目标(终点碳、终点温度)，
用物料平衡算配料(石灰/白云石)与耗氧量，再由热平衡差额反算合理废钢量——废钢即冷却剂，
热量富余多则能化更多废钢。与动态预测器分工：静态段吹炼前配料，动态段吹炼中预测终点。

热平衡常数复用 config / thermo.constants，与动态模型同源，保证首尾衔接一致
(石灰分解热等待 D 项统一修正后两段自动同步)。
"""
import numpy as np
from config import CP_STEEL, CP_SLAG, SCRAP_CP, SCRAP_LATENT, CALIBRATION_DEFAULTS
from thermo.constants import (MOLAR_MASS, OXIDATION_THERMO, HEAT_C2S_FORMATION,
                              HEAT_C3P_FORMATION, HEAT_LIME_DECOMPOSITION,
                              HEAT_DOLOMITE_DECOMPOSITION, CP_LIME_SOLID, LIME_LOI_FRACTION)

LIME_COMP = {'CaO': 90.0, 'MgO': 3.0, 'SiO2': 2.0, 'Al2O3': 1.0, 'S': 0.05, '烧损': 3.95}
DOLOMITE_COMP = {'CaO': 30.0, 'MgO': 21.0, 'SiO2': 2.0, 'Al2O3': 1.0, '烧损': 46.0}
SINTER_COMP = {'FeO': 9.0, 'Fe2O3': 64.0, 'CaO': 13.0, 'SiO2': 7.0, 'MgO': 2.0, 'Al2O3': 2.0, 'MnO': 1.5, '其他': 1.5}

T_AMB = 298.0     # 环境基准温度 K (废钢/辅料入炉温度)
CP_CO = 33.0      # CO 摩尔比热 J/(mol·K)，与动态模型 Q_offgas 一致
Q_RAD_W = 2.0e6   # 辐射+炉衬散热 W (标称值，与动态默认一致)
ETA_POST = CALIBRATION_DEFAULTS['eta_post']   # 二次燃烧比例 (标称 0.10)


class StaticModel:
    def __init__(self):
        self.cp_steel = CP_STEEL
        self.cp_slag = CP_SLAG

    def calculate(self, inputs: dict) -> dict:
        m_hm = inputs['铁水质量'] * 1000.0
        T_hm = inputs['铁水温度'] + 273.15
        T_end = inputs['终点温度目标'] + 273.15
        R = inputs['炉渣碱度']
        hm = {k: inputs[k] for k in ['w([C])', 'w([Si])', 'w([Mn])', 'w([P])',
                                     'w([S])', 'w([V])', 'w([Ti])']}

        # ---------- 物料平衡：元素氧化量 → 耗氧 → 渣量 → 辅料 ----------
        oxided = {}
        for elem in ['Si', 'Mn', 'P', 'V', 'Ti']:
            oxided[elem] = hm[f'w([{elem}])'] / 100.0 * m_hm
        dC = hm['w([C])'] - inputs['终点碳目标']
        oxided['C'] = max(0.0, dC / 100.0 * m_hm)
        # 铁氧化率(占铁水%): 决定烧铁产热/铁损, 是废钢容量的关键杠杆。
        # 默认 1%(张强低废钢比); 大废钢比高FeO操作靠多烧铁化废钢, 约 4~6%
        # (本数据集实测出钢量反推约 4.8%)。有实测出钢量时优先用它反推真实铁氧化。
        fe_ox_pct = inputs.get('预期铁氧化率', 1.0)
        oxided['Fe'] = fe_ox_pct / 100.0 * m_hm

        O2_kg = 0.0
        for elem, mass_kg in oxided.items():
            O2_kg += mass_kg / MOLAR_MASS[elem] * OXIDATION_THERMO[elem]['n'] * 0.016
        O2_actual = O2_kg / 0.95 / 1.429   # Nm³ (95%利用率, O₂密度1.429 kg/Nm³)

        slag_oxides = {'SiO2': 0.0, 'MnO': 0.0, 'P2O5': 0.0, 'FeO': 0.0}
        conv = {'Si': ('SiO2', 0.06009), 'Mn': ('MnO', 0.07094),
                'P': ('P2O5', 0.14194 / 2.0), 'Fe': ('FeO', 0.07185)}
        for elem, (ox, M_ox) in conv.items():
            if oxided.get(elem, 0) > 0:
                slag_oxides[ox] += oxided[elem] / MOLAR_MASS[elem] * M_ox

        CaO_need = R * slag_oxides['SiO2']
        lime_kg = max(0.0, CaO_need / (LIME_COMP['CaO'] / 100.0))
        slag_mass_est = sum(slag_oxides.values()) + lime_kg * (LIME_COMP['CaO'] + LIME_COMP['MgO']) / 100.0
        MgO_need = inputs['MgO目标含量'] / 100.0 * slag_mass_est
        dolomite_kg = max(0.0, (MgO_need - lime_kg * LIME_COMP['MgO'] / 100.0) / (DOLOMITE_COMP['MgO'] / 100.0))
        m_slag = sum(slag_oxides.values()) + lime_kg * 0.96 + dolomite_kg * 0.54   # 渣总量(扣烧损)

        # ---------- 热量平衡：反算合理废钢量 ----------
        # 热收入(化学热): 元素氧化放热 + C₂S/C₃P 形成热 + 二次燃烧
        Q_ox = sum(oxided[e] / MOLAR_MASS[e] * OXIDATION_THERMO[e]['delta_H'] for e in oxided)
        Q_C2S = oxided['Si'] / MOLAR_MASS['Si'] * HEAT_C2S_FORMATION
        Q_C3P = oxided['P'] / MOLAR_MASS['P'] / 2.0 * HEAT_C3P_FORMATION
        Q_postcomb = ETA_POST * oxided['C'] / MOLAR_MASS['C'] * OXIDATION_THERMO['C']['delta_H']
        Q_chem = Q_ox + Q_C2S + Q_C3P + Q_postcomb

        # 热支出(与废钢量无关的部分):
        Q_hm_heating = m_hm * self.cp_steel * (T_end - T_hm)          # 铁水升温到终点
        # 炉渣升温: 仅白云石分解产物(0.53)从1100K升温; 石灰升温含在 q_lime_per_kg (298K起);
        # 氧化产物(SiO2/FeO等)由熔池内高温元素原位生成、无需从冷加热 → 与动态模型一致
        Q_slag = 0.0
        if T_end > 1100.0:
            Q_slag = dolomite_kg * 0.53 * self.cp_slag * (T_end - 1100.0)
        # 石灰为煅烧CaO: 显热(298K→终点)+烧损部分分解热; 白云石为生料付全额分解热 (D修复)
        q_lime_per_kg = CP_LIME_SOLID * (T_end - T_AMB) + LIME_LOI_FRACTION * HEAT_LIME_DECOMPOSITION
        Q_flux = lime_kg * q_lime_per_kg + dolomite_kg * HEAT_DOLOMITE_DECOMPOSITION
        n_CO = oxided['C'] / MOLAR_MASS['C']
        # CO 全程陆续离炉, 显热按释放时刻熔池温度计 —— 取全程平均 (T_hm+T_end)/2 近似
        # (动态模型逐步按当时温度计, 此为静态端的一致化; 修正前按终点温度计全部 CO,
        # 高估 ~4-5 GJ, 是 13.13/13.18 两次诊断都指向的"烟气偏高"尾巴)
        T_offgas_avg = 0.5 * (T_hm + T_end)
        Q_offgas = n_CO * CP_CO * (T_offgas_avg - T_AMB)             # 烟气带走显热
        blow_time_s = O2_actual / max(inputs['供氧流量'], 1.0) * 60.0
        Q_loss = Q_RAD_W * blow_time_s + 0.01 * Q_chem               # 辐射炉衬 + 粉尘(~1%化学热)

        # 单位废钢吸热 (从环境升温到终点 + 熔化潜热), 与动态模型 SCRAP_CP/LATENT 一致
        q_scrap_per_kg = SCRAP_CP * (T_end - T_AMB) + SCRAP_LATENT
        # 平衡: Q_chem = Q_hm_heating + Q_slag + Q_flux + Q_offgas + Q_loss + m_scrap·q_scrap_per_kg
        surplus = Q_chem - Q_hm_heating - Q_slag - Q_flux - Q_offgas - Q_loss
        m_scrap_reasonable = max(0.0, surplus / q_scrap_per_kg)      # kg

        # ---------- 出钢量与铁耗 (用反算废钢) ----------
        m_ox_total = sum(oxided.values())
        steel_output_t = (m_hm + m_scrap_reasonable - m_ox_total) / 1000.0
        scrap_reasonable_t = m_scrap_reasonable / 1000.0
        iron_loss = (inputs['铁水质量'] + scrap_reasonable_t) / max(steel_output_t, 1e-6) * 1000.0

        result = {
            '总耗氧量': round(O2_actual, 0),
            '石灰加入量': round(lime_kg, 0),
            '白云石加入量': round(dolomite_kg, 0),
            '合理废钢量': round(scrap_reasonable_t, 2),
            '出钢量': round(steel_output_t, 2),
            '铁耗': round(iron_loss, 2),
            '热平衡明细_MJ': {
                '化学热收入': round(Q_chem / 1e6, 1),
                '铁水升温': round(Q_hm_heating / 1e6, 1),
                '炉渣升温': round(Q_slag / 1e6, 1),
                '辅料分解': round(Q_flux / 1e6, 1),
                '烟气显热': round(Q_offgas / 1e6, 1),
                '辐射炉衬粉尘': round(Q_loss / 1e6, 1),
                '可化废钢的富余热': round(surplus / 1e6, 1),
            },
        }

        # ---------- 校验：计算值 vs 真实值 (若数据提供) ----------
        chk = {}
        if inputs.get('废钢加入量') is not None:
            chk['废钢_实际t'] = round(inputs['废钢加入量'], 2)
            chk['废钢_反算t'] = round(scrap_reasonable_t, 2)
            chk['废钢_偏差t'] = round(scrap_reasonable_t - inputs['废钢加入量'], 2)
        if inputs.get('实际耗氧量Nm3') is not None:
            chk['耗氧_实际Nm3'] = round(inputs['实际耗氧量Nm3'], 0)
            chk['耗氧_计算Nm3'] = round(O2_actual, 0)
            chk['耗氧_偏差Nm3'] = round(O2_actual - inputs['实际耗氧量Nm3'], 0)
        if inputs.get('实际石灰kg') is not None:
            chk['石灰_实际kg'] = round(inputs['实际石灰kg'], 0)
            chk['石灰_计算kg'] = round(lime_kg, 0)
        if chk:
            result['校验'] = chk
        return result
