import pytest
from model.lime_melting import LimeMelting


def test_calc_k_prime():
    lm = LimeMelting()
    kp = lm.calc_k_prime(rho_s=3500.0, eta=0.5, D=1e-9, d=0.01, U=0.5)
    assert kp > 0


def test_calc_dr_dt():
    lm = LimeMelting()
    drdt = lm.calc_dr_dt(k_S=1e-5, rho_S=3500.0, rho_CaO=3300.0, w_CaO_sat=55.0, w_CaO_slag=40.0)
    assert drdt > 0


def test_calc_delta_m_lime():
    lm = LimeMelting()
    dm = lm.calc_delta_m(m_BD=10.0, r=0.01, delta_r=0.001)
    assert dm > 0
    assert dm < 10.0


def test_c2s_treatment():
    lm = LimeMelting()
    kp = lm.calc_k_prime_C2S(k_prime_normal=1e-4)
    assert kp == 1e-5
