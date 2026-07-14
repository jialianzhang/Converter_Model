import pytest
import numpy as np
from thermo.wagner import calc_lgf_i, calc_f_i


def test_calc_lgf_C():
    # 铁水(4.5%C)中 f_C≈4~6 (lgf_C≈0.6~0.7), Sigworth-Elliott 系数
    composition = {'C': 4.5, 'Si': 0.34, 'Mn': 0.19, 'P': 0.117, 'S': 0.01}
    lgf_C = calc_lgf_i('C', composition, 1873.0)
    assert 0.3 < lgf_C < 1.0


def test_f_O_suppressed_in_high_C_melt():
    # 高碳铁液中溶解氧活度系数被强烈压低: f_O≈0.01~0.05
    composition = {'C': 4.5, 'Si': 0.34, 'Mn': 0.19, 'P': 0.117, 'S': 0.01}
    f_O = calc_f_i('O', composition, 1873.0)
    assert 0.005 < f_O < 0.1


def test_calc_f_i_positive():
    composition = {'C': 4.5, 'Si': 0.34, 'Mn': 0.19, 'P': 0.117}
    f_C = calc_f_i('C', composition, 1873.0)
    assert f_C > 0


def test_zero_composition():
    composition = {'C': 0.0, 'Si': 0.0, 'Mn': 0.0, 'P': 0.0}
    lgf = calc_lgf_i('Fe', composition, 1873.0)
    assert abs(lgf) < 1e-10
