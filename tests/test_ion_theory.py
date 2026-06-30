import pytest
import numpy as np
from thermo.ion_theory import calc_slag_activity_coeff


def test_calc_FeO_activity():
    slag_comp = {'FeO': 20.0, 'SiO2': 15.0, 'CaO': 45.0, 'MnO': 5.0, 'MgO': 8.0, 'P2O5': 3.0}
    f_FeO = calc_slag_activity_coeff('FeO', slag_comp, 1873.0)
    assert 0.1 < f_FeO < 10.0


def test_calc_SiO2_activity():
    slag_comp = {'FeO': 20.0, 'SiO2': 15.0, 'CaO': 45.0, 'MnO': 5.0, 'MgO': 8.0}
    f_SiO2 = calc_slag_activity_coeff('SiO2', slag_comp, 1873.0)
    assert 0.05 < f_SiO2 < 10.0
