"""反应平衡常数计算: lg K_M = A/T + B"""
import numpy as np
from .constants import OXIDATION_THERMO

def calc_lgK(element: str, T: float) -> float:
    """计算元素氧化反应的 lg K_M"""
    if element not in OXIDATION_THERMO:
        raise ValueError(f"Unknown element: {element}")
    d = OXIDATION_THERMO[element]
    return d['A'] / T + d['B']

def calc_K_M(element: str, T: float) -> float:
    """计算平衡常数 K_M"""
    return 10.0 ** calc_lgK(element, T)

def calc_all_lgK(T: float) -> dict:
    """计算所有元素的 lg K_M"""
    return {elem: calc_lgK(elem, T) for elem in OXIDATION_THERMO}

def calc_all_K_M(T: float) -> dict:
    """计算所有元素的 K_M"""
    return {elem: 10.0 ** lgK for elem, lgK in calc_all_lgK(T).items()}
