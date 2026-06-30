"""炉渣相活度计算 - 离子结构理论 (张强 式4)"""
import numpy as np
from config import R_GAS

EXCHANGE_ENERGY = {
    ('SiO2', 'FeO'):  -50000.0, ('SiO2', 'MnO'):  -30000.0,
    ('SiO2', 'CaO'): -120000.0, ('SiO2', 'MgO'):  -80000.0,
    ('SiO2', 'P2O5'): -40000.0,
    ('FeO',  'CaO'):  -30000.0, ('FeO',  'MnO'):  10000.0,
    ('CaO',  'MnO'):  -20000.0, ('CaO',  'MgO'):  -15000.0,
    ('MnO',  'MgO'):   -5000.0, ('P2O5', 'CaO'): -180000.0,
}


def get_epsilon(ox_i: str, ox_j: str) -> float:
    key = tuple(sorted([ox_i, ox_j]))
    if key in EXCHANGE_ENERGY:
        return EXCHANGE_ENERGY[key]
    return 0.0


def calc_slag_activity_coeff(oxide: str, slag_comp: dict, T: float) -> float:
    denominator = 0.0
    for ox_j, w_j in slag_comp.items():
        if w_j <= 0:
            continue
        eps = get_epsilon(oxide, ox_j)
        denominator += w_j * np.exp(-eps / (R_GAS * T))
    if denominator < 1e-20:
        return 1.0
    return 1.0 / denominator


def calc_all_slag_activity_coeff(slag_comp: dict, T: float) -> dict:
    return {ox: calc_slag_activity_coeff(ox, slag_comp, T) for ox in slag_comp}
