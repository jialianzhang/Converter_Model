import pytest
import numpy as np
from model.modified_mass_transfer import ModifiedMassTransfer

def test_calc_k_m():
    mt = ModifiedMassTransfer()
    km = mt.calc_k_m(km0=2e-4, eps_total=50.0, eps_total0=50.0, gamma=0.3)
    assert abs(km - 2e-4) < 1e-8

def test_calc_k_m_increases_with_energy():
    mt = ModifiedMassTransfer()
    km1 = mt.calc_k_m(km0=2e-4, eps_total=100.0, eps_total0=50.0, gamma=0.3)
    km2 = mt.calc_k_m(km0=2e-4, eps_total=50.0, eps_total0=50.0, gamma=0.3)
    assert km1 > km2

def test_calc_F_M():
    mt = ModifiedMassTransfer()
    FM = mt.calc_F_M(km=2e-4, rho_m=7000.0, M_M=0.02809)
    assert 0.01 < FM < 10.0

def test_ks_lower_bound():
    mt = ModifiedMassTransfer()
    ks = mt.calc_k_s(ks_prime=1e-8, f_solid=0.99, n=5)
    assert ks >= 1e-8
