import pytest
from thermo.slag_phase import calc_CaO_saturation, calc_solid_fraction


def test_calc_CaO_saturation():
    w_CaO_sat = calc_CaO_saturation(w_FeO=15.0, w_SiO2=15.0)
    assert 35.0 < w_CaO_sat < 65.0


def test_calc_CaO_saturation_high_FeO():
    w_CaO_sat_high = calc_CaO_saturation(w_FeO=30.0, w_SiO2=10.0)
    w_CaO_sat_low = calc_CaO_saturation(w_FeO=5.0, w_SiO2=20.0)
    assert w_CaO_sat_high < w_CaO_sat_low


def test_calc_solid_fraction_zero_when_unsaturated():
    fs = calc_solid_fraction(w_CaO=30.0, w_SiO2=20.0, w_FeO=15.0)
    assert fs == 0.0


def test_calc_solid_fraction_nonzero_when_saturated():
    w_CaO_sat = calc_CaO_saturation(w_FeO=15.0, w_SiO2=15.0)
    fs = calc_solid_fraction(w_CaO=w_CaO_sat + 10.0, w_SiO2=15.0, w_FeO=15.0)
    assert fs > 0.0
