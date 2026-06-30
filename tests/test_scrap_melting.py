import pytest
from model.scrap_melting import ScrapMelting


def test_calc_T_MP():
    sm = ScrapMelting()
    T_mp = sm.calc_T_MP(w_C=4.34)
    assert 1400 < T_mp < 1440


def test_calc_T_MP_low_C():
    sm = ScrapMelting()
    T_mp = sm.calc_T_MP(w_C=0.1)
    assert T_mp > 1700


def test_melting_mechanism_diffusion():
    sm = ScrapMelting()
    mech = sm.determine_mechanism(T_bath=1400.0, T_MP=1500.0)
    assert mech == 'diffusion'


def test_melting_mechanism_forced():
    sm = ScrapMelting()
    mech = sm.determine_mechanism(T_bath=1600.0, T_MP=1500.0)
    assert mech == 'forced'
