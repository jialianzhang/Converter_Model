"""炉渣相图简化模型：CaO饱和浓度与固相率 (CaO-SiO2-FeO 系 C2S 区域)
以及 FeO/MnO 的 Raoultian 活度系数 (文献实测值分段插值)、Riboud 渣黏度模型"""
import numpy as np
from .constants import OXIDE_MOLAR_MASS


def calc_slag_viscosity(slag_wt: dict, T: float, f_solid: float = 0.0) -> float:
    """液态渣黏度 Riboud(1981) 模型 + Einstein-Roscoe 固相悬浮修正, 返回 Pa·s

    Riboud: η = A·T·exp(B/T)
      A = exp(−19.81 + 1.73·X_碱),  B = 31140 − 23896·X_碱
      X_碱 = X_CaO + X_MgO + X_FeO + X_MnO (摩尔分数; 本渣系无 CaF₂/Na₂O/Al₂O₃ 项)
    固相修正 (Roscoe): η_eff = η·(1 − 1.35·f_s)^−2.5 —— C₂S 析出使渣变稠(返干趋势)

    D 修复: 取代无出处的自造式 0.5·exp(8000/T)·(1+10f_s) —— 后者给出 35~74 Pa·s,
    高于真实转炉渣 (全液态 0.02~0.5 Pa·s) 约两个数量级, 曾把石灰溶解速率压低 ~5 倍。
    """
    mol = {}
    for ox, w in slag_wt.items():
        if ox in OXIDE_MOLAR_MASS and w > 0:
            mol[ox] = w / OXIDE_MOLAR_MASS[ox]
    tot = sum(mol.values())
    if tot <= 0 or T <= 0:
        return 0.1
    x_basic = sum(mol.get(ox, 0.0) for ox in ('CaO', 'MgO', 'FeO', 'MnO')) / tot
    A = np.exp(-19.81 + 1.73 * x_basic)
    B = 31140.0 - 23896.0 * x_basic
    eta = A * T * np.exp(B / T)
    fs = min(max(f_solid, 0.0), 0.7)
    eta *= (1.0 - 1.35 * fs) ** -2.5
    return float(np.clip(eta, 0.005, 50.0))


def calc_gamma_FeO(basicity: float) -> float:
    """碱性渣中 FeO 的 Raoultian 活度系数 gamma_FeO(x基准)

    文献实测锚点 (Timucin & Morris 1970; Suito & Inoue 1984 系列):
    B<=0.5: ~1.0; B=1: ~1.4; B=1.5: ~2.2; B=2~3: ~3.0; B>=4: ~2.5
    分段线性插值, 可标定。
    """
    pts_B = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
    pts_g = [1.0, 1.4, 2.2, 3.0, 3.0, 2.5]
    return float(np.clip(np.interp(basicity, pts_B, pts_g), 1.0, 3.5))


def calc_gamma_MnO(basicity: float) -> float:
    """碱性渣中 MnO 的 Raoultian 活度系数 (Suito-Inoue 量级 ~1-2, 随碱度升高)"""
    return float(np.clip(0.5 + 0.5 * basicity, 0.5, 2.0))


def calc_L_P(T: float, slag_wt: dict) -> float:
    """磷在渣-钢间的平衡分配比 L_P = (%P)_slag / [%P]_metal

    Suito & Inoue (1984) 经验式 (Kitamura 2009 脱磷模型采用的同源关系):
    log L_P = 0.072·[(%CaO) + 0.3(%MgO) + 0.6(%P2O5) + 0.2(%MnO)]
              + 2.5·log(%T.Fe) + 11570/T − 10.52
    %T.Fe: 渣中全铁 (本模型铁氧化物全部记作 FeO, T.Fe = 0.777·%FeO)

    物理内涵: 氧势由渣中铁氧化物承载(2.5·log T.Fe 项), 磷酸盐容量由
    CaO/MgO 等碱性组元承载 —— 正是"液滴/渣界面氧势接近渣侧"的乳化区脱磷图像。
    """
    w_TFe = max(slag_wt.get('FeO', 0.0) * 55.85 / 71.85, 0.1)
    log_LP = (0.072 * (slag_wt.get('CaO', 0.0) + 0.3 * slag_wt.get('MgO', 0.0)
                       + 0.6 * slag_wt.get('P2O5', 0.0) + 0.2 * slag_wt.get('MnO', 0.0))
              + 2.5 * np.log10(w_TFe) + 11570.0 / T - 10.52)
    return float(np.clip(10.0 ** log_LP, 1e-3, 1e6))


def calc_CaO_saturation(w_FeO: float, w_SiO2: float) -> float:
    if w_SiO2 < 0.01:
        return 55.0
    return 45.0 + 0.7 * w_SiO2 - 0.15 * w_FeO


def calc_solid_fraction(w_CaO: float, w_SiO2: float, w_FeO: float) -> float:
    w_CaO_sat = calc_CaO_saturation(w_FeO, w_SiO2)
    if w_CaO <= w_CaO_sat:
        return 0.0
    w_CaO_C2S = 65.0
    f_solid = (w_CaO - w_CaO_sat) / (w_CaO_C2S - w_CaO_sat)
    return np.clip(f_solid, 0.0, 0.8)
