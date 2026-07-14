"""全局配置与物理常数"""
import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT_DIR, 'data')

# 钢液/炉渣物理常数
RHO_STEEL = 7000.0       # 钢液密度 kg/m³
RHO_SLAG = 3500.0        # 炉渣密度 kg/m³
CP_STEEL = 820.0         # 钢液比热容 J/(kg·K)
CP_SLAG = 1200.0         # 炉渣比热容 J/(kg·K)
R_GAS = 8.314            # 气体常数 J/(mol·K)
T_REF = 1873.0           # 参考温度 K

# 炉渣工艺参数
SLAG_BASICITY_TARGET = 3.0   # 目标碱度 CaO/SiO₂
MGO_TARGET = 8.0             # 目标 MgO %

# Kitamura基准参数 (Kitamura 2009, Table 2; k_s0/G_CO0 为工业尺度标定参数)
KITAMURA_BASE = {
    'k_m0': 2.0e-4,       # 基准钢液传质系数 m/s
    'k_s0': 1.0e-5,       # 基准炉渣传质系数 m/s (配合 SLAG_REF_MASS 归一, 标定参数)
    'q_Ar0': 0.00033,     # Kitamura实验基准氩气流量 Nm³/s (70kg炉, 仅存档)
    'W_s0': 2.0,          # Kitamura实验基准炉渣质量 kg (70kg炉, 仅存档)
    'G_CO0': 1.0,         # 基准CO生成速率 mol/(s·m²); 上界≈2×O2摩尔流/A≈1.7, 标定参数
    'A0': 10.0,           # 基准反应界面面积 m²
    'n_solid': 5,         # 固相率衰减指数
}

# 炉渣传质的工业尺度基准渣量 (替代 40/m_slag 无量纲因子; k_s' = k_s0·stir·(SLAG_REF_MASS/m_slag))
SLAG_REF_MASS = 10000.0   # kg, 标定参数

# 冲击区(火点)O₂分配份额 (氧账本: 其余全部 [Fe]+1/2O₂→(FeO) 供渣-金界面耦合反应)
# 注: 火点直接脱碳份额 C_direct 已移入 CALIBRATION_DEFAULTS (标定参数, 引擎从 self.p 读取)
O2_SHARE = {
    'Si': 0.08,           # Si 存在期间的火点直接氧化份额 (耗尽后份额转给FeO)
    'Mn': 0.02,           # Mn 同上
    'escape': 0.05,       # 二次燃烧耗氧/烟气逸散等未入熔池部分
}

# 点火爬坡时间 s: 开吹 O2 流量从 0 线性升到额定的点火期 (真实操作 ~30-90s)。
# 服务前期温度/成分曲线形态 (0-1min 废钢吸热微降得以显现); 0 = 关闭
O2_RAMP_TIME = 60.0

# 乳化液滴有效反应面积 m² (FeO还原/脱碳/脱磷的主要场所; 标定参数)
# 标定依据(非第一性原理): 转炉中期渣 FeO 存量台阶 8~15% 这一工业公认事实
# (张强图4c 谷底~13%)。准稳态 FeO = w_FeO_eq + 生成速率/(A_total·F_FeO),
# 生成~760 mol O/s, F_FeO~0.06-0.10 → A_total≈1050 m² → A_em≈300。
# 物理内涵: 有效面积受液滴携碳量限制(bloated droplet, 液滴飞行中自身碳耗尽),
# 远小于液滴几何总面积; 第四阶段用 Subagyo-Brooks 吹入数 N_B + 停留时间动态计算替代
A_EMULSION = 300.0

# 废钢热参数 (机理化热汇: 升温显热 + 熔化潜热, 替代1.3MJ/kg集总值)
SCRAP_CP = 700.0          # 固态废钢平均比热 J/(kg·K)
SCRAP_LATENT = 272000.0   # 熔化潜热 J/kg (纯铁~272 kJ/kg)
SCRAP_TAU_HEAT = 180.0    # 未熔废钢升温一阶弛豫时间常数 s (对应h·A量级, 标定参数)

# 标定参数默认值 —— V1.0 冻结值 (2026-07-09 冻结, 100 炉大废钢比数据标定结果)
# 一组参数共享于全部炉次 (设备/工艺级常数, 炉次差异由输入承载); Engine() 即生产模型。
# 标定框架=目标碳终止; 泛化验证=验证集≈训练集 (T MAE 24.6°C, Mn±0.02 命中 73%)。
# 原始默认值(标定前, alpha=1/gamma0.3/eta_post0.10/Q_rad2e6/A_em300/C_direct0.30/
#   eta_late1.0/tau180/scale1.0/mult2.0)见 git 历史 v6.x; 备份亦存 v1_frozen_params.json。
CALIBRATION_DEFAULTS = {
    'alpha_steel': 0.77182,   # 钢液传质系数乘子: k_m0_eff = alpha_steel × k_m0
    'alpha_slag': 1.16068,    # 炉渣传质系数乘子: k_s0_eff = alpha_slag × k_s0
    'G_CO0': 1.22544,         # CO生成速率基准 mol/(s·m²)
    'gamma': 0.44005,         # 搅拌指数
    'eta_post': 0.10035,      # 二次燃烧比例
    'Q_rad_W': 2224473.16,    # 辐射+炉衬综合散热 W
    'A_emulsion': 317.5931,   # 乳化液滴有效反应面积 m²
    'C_direct': 0.24952,      # 火点直接脱碳份额
    'eta_O2_late': 0.45,      # 末期"入铁氧保留率"下限 (由终渣FeO工业范围物理定值, 标定中固定)
    'scrap_tau_heat': 209.06472,  # 废钢升温弛豫时间常数 s
    'scrap_heat_scale': 0.60007,  # 废钢热汇缩放 (触下界, 吸收废钢比结构残差)
    'aO_cap_CO_mult': 2.04938,    # 碳氧缓冲本体氧系数 (Mn 平衡地板用, 不动全局氧平衡)
}

# 标定参数物理范围 (优化器边界约束)
CALIBRATION_BOUNDS = {
    'alpha_steel': (0.3, 3.0),
    'alpha_slag': (0.3, 3.0),
    'G_CO0': (0.1, 5.0),
    'gamma': (0.1, 0.8),
    'eta_post': (0.05, 0.30),
    'Q_rad_W': (0.5e6, 3.0e6),
    'A_emulsion': (100.0, 1500.0),
    'C_direct': (0.15, 0.45),
    'eta_O2_late': (0.4, 1.0),
    'scrap_tau_heat': (60.0, 600.0),
    'scrap_heat_scale': (0.6, 1.3),
    'aO_cap_CO_mult': (1.0, 4.0),
}

# 模拟控制
TIME_STEP_INTENSE = 1.0   # 脱碳剧烈期步长 s, w([C])>0.8%
TIME_STEP_NORMAL = 2.0    # 常规平稳期步长 s, 0.3%<w([C])≤0.8%
TIME_STEP_END = 5.0       # 吹炼末期步长 s, w([C])≤0.3%
SOLVER_TOL = 1e-10        # Brent寻根收敛容差
K_S_MIN = 1e-8            # ks数值下界 m/s
