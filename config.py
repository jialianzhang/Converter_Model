"""全局配置与物理常数"""
import os

# 项目根目录
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT_DIR, 'data')

# 钢液/炉渣物理常数
RHO_STEEL = 7000.0       # 钢液密度 kg/m³
RHO_SLAG = 3500.0        # 炉渣密度 kg/m³
CP_STEEL = 820.0         # 钢液比热容 J/(kg·K)
CP_SLAG = 1200.0         # 炉渣比热容 J/(kg·K)
R_GAS = 8.314            # 气体常数 J/(mol·K)
T_REF = 1873.0           # 参考温度 K (1600°C)

# 炉渣参数
SLAG_BASICITY_TARGET = 3.0   # 目标碱度 CaO/SiO₂
MGO_TARGET = 8.0             # 目标 MgO %

# Kitamura 标准条件参数 (Table 2)
KITAMURA_BASE = {
    'k_m0': 2.0e-4,       # 基准钢液传质系数 m/s
    'k_s0': 5.0e-4,       # 基准炉渣传质系数 m/s (工业转炉CO+O2搅拌)
    'q_Ar0': 0.00033,     # 基准氩气流量 Nm³/s (≈20 NL/min)
    'W_s0': 2.0,          # 基准炉渣质量 kg (渣钢比≈2%)
    'G_CO0': 0.01,        # 基准CO生成速率 mol/(s·m²)
    'A0': 10.0,           # 基准反应界面面积 m²
    'n_solid': 5,         # 固相率衰减指数
}

# 标定参数默认值
CALIBRATION_DEFAULTS = {
    'alpha_steel': 1.0,
    'alpha_slag': 1.0,
    'G_CO0': 0.01,
    'gamma': 0.3,
    'eta_post': 0.10,
    'alpha_rad': 0.03,
    'epsilon_dust': 0.01,
}

# 模拟控制
TIME_STEP_INTENSE = 1.0   # 脱碳剧烈期步长 s
TIME_STEP_NORMAL = 2.0    # 常规期步长 s
TIME_STEP_END = 5.0       # 末期步长 s
SOLVER_TOL = 1e-10        # 寻根收敛容差
K_S_MIN = 1e-8            # ks 数值下界
