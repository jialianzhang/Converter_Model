import pytest
import numpy as np
from thermo.wagner import calc_lgf_i, calc_f_i


def test_calc_lgf_C():
    composition = {'C': 4.5, 'Si': 0.34, 'Mn': 0.19, 'P': 0.117, 'S': 0.01}
    lgf_C = calc_lgf_i('C', composition, 1873.0)
    assert 0.001 < lgf_C < 0.05


def test_calc_f_i_positive():
    composition = {'C': 4.5, 'Si': 0.34, 'Mn': 0.19, 'P': 0.117}
    f_C = calc_f_i('C', composition, 1873.0)
    assert f_C > 0


def test_zero_composition():
    composition = {'C': 0.0, 'Si': 0.0, 'Mn': 0.0, 'P': 0.0}
    lgf = calc_lgf_i('Fe', composition, 1873.0)
    assert abs(lgf) < 1e-10
