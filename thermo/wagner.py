"""金属相活度计算 - Wagner 模型: lg f_i = Sum e_i^j * w([j])"""
import numpy as np
from config import T_REF, R_GAS

WAGNER_E_1E4 = {
    'C':  {'C': 14.0, 'Si': 8.0,  'Mn': -1.2, 'P': 5.1,  'S': 4.6,  'O': -34.0},
    'Si': {'C': 18.0, 'Si': 11.0, 'Mn': 0.2,  'P': 11.0, 'S': 5.6,  'O': -23.0},
    'Mn': {'C': -7.0, 'Si': 0.0,  'Mn': 0.0,  'P': -0.35,'S': -4.8, 'O': -8.3},
    'P':  {'C': 13.0, 'Si': 12.0, 'Mn': 0.0,  'P': 6.2,  'S': 2.8,  'O': 13.0},
    'S':  {'C': 11.0, 'Si': 6.3,  'Mn': -2.6, 'P': 2.9,  'S': -2.8, 'O': -27.0},
    'V':  {'C': -18.0,'Si': 4.2,  'Mn': 0.0,  'P': 0.0,  'S': -28.0,'O': -70.0},
    'Ti': {'C': -18.0,'Si': 5.5,  'Mn': 0.0,  'P': 0.0,  'S': -28.0,'O': -85.0},
    'Fe': {'C': 0.0,  'Si': 0.0,  'Mn': 0.0,  'P': 0.0,  'S': 0.0,  'O': 0.0},
    'O':  {'C': -45.0, 'Si': -13.1,'Mn': -2.1, 'P': 7.0,  'S': -13.3,'O': -20.0},
}


def get_e_ij(element_i: str, element_j: str, T: float) -> float:
    e0 = WAGNER_E_1E4.get(element_i, {}).get(element_j, 0.0) * 1e-4
    return e0 * T_REF / T


def calc_lgf_i(element: str, composition: dict, T: float) -> float:
    lgf = 0.0
    for elem_j, w_j in composition.items():
        if w_j > 0:
            lgf += get_e_ij(element, elem_j, T) * w_j
    return lgf


def calc_f_i(element: str, composition: dict, T: float) -> float:
    return 10.0 ** calc_lgf_i(element, composition, T)


def calc_all_f_i(composition: dict, T: float) -> dict:
    return {elem: calc_f_i(elem, composition, T) for elem in composition}
