import pytest
from model.mixing_energy import MixingEnergy

def test_calc_epsilon_T():
    me = MixingEnergy()
    eps_T = me.calc_epsilon_T(Q_T=850.0, W=170000.0, X=1.5, D=0.05, U=300.0, theta=12.0)
    assert eps_T > 0

def test_calc_epsilon_B():
    me = MixingEnergy()
    eps_B = me.calc_epsilon_B(Q_B=11.5, T=1549.0, W=170000.0, H=1.8)
    assert eps_B > 0

def test_total_mixing_energy():
    me = MixingEnergy()
    eps_T = me.calc_epsilon_T(850.0, 170000.0, 1.5, 0.05, 300.0, 12.0)
    eps_B = me.calc_epsilon_B(11.5, 1549.0, 170000.0, 1.8)
    total = me.calc_total(850.0, 11.5, 1549.0, 170000.0, 1.5, 0.05, 300.0, 12.0, 1.8)
    assert abs(total - (eps_T + eps_B)) < 1e-6
