import pytest
from model.oxygen_balance import OxygenBalance


def test_calc_a_O_b_from_slag():
    ob = OxygenBalance()
    aOb = ob.calc_a_O_b_from_slag(w_FeO_b=20.0, f_FeO=1.0, K_Fe=100.0)
    assert aOb > 0


def test_calc_J_M():
    ob = OxygenBalance()
    J = ob.calc_J_M(F_M=0.5, w_M_b=4.0, w_M_star=1.0)
    assert J == 1.5  # 0.5 * (4.0 - 1.0)


def test_residual():
    ob = OxygenBalance()
    J_fluxes = {'Si': 0.1, 'Mn': 0.2, 'P': 0.05, 'V': 0.0, 'Ti': 0.0, 'Fe': 0.3, 'C': 0.5}
    residual = ob.calc_residual(J_fluxes, F_O=5.0, a_O_b_slag=0.01, a_O_star=1e-4, f_O=1.0)
    assert isinstance(residual, float)


def test_calc_J_O_positive_driving_force():
    ob = OxygenBalance()
    JO = ob.calc_J_O(F_O=5.0, a_O_b=0.01, a_O_star=1e-6, f_O=1.0)
    assert JO > 0
