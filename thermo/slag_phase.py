"""炉渣相图简化模型：CaO饱和浓度与固相率 (CaO-SiO2-FeO 系 C2S 区域)"""
import numpy as np


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
