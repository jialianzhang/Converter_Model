# tests/test_constants.py
import pytest
from thermo.constants import OXIDATION_THERMO, MOLAR_MASS, OXIDE_MOLAR_MASS

def test_oxidation_thermo_has_seven_elements():
    assert 'Si' in OXIDATION_THERMO
    assert 'Mn' in OXIDATION_THERMO
    assert 'P' in OXIDATION_THERMO
    assert 'V' in OXIDATION_THERMO
    assert 'Ti' in OXIDATION_THERMO
    assert 'Fe' in OXIDATION_THERMO
    assert 'C' in OXIDATION_THERMO

def test_oxidation_thermo_format():
    for elem, data in OXIDATION_THERMO.items():
        assert 'A' in data
        assert 'B' in data
        assert 'n' in data
        assert 'delta_H' in data

def test_C_oxidation_is_CO():
    assert OXIDATION_THERMO['C']['n'] == 1

def test_P_oxidation_is_P2O5():
    assert OXIDATION_THERMO['P']['n'] == 2.5

def test_molar_mass():
    assert MOLAR_MASS['Fe'] == 0.05585
    assert MOLAR_MASS['C'] == 0.01201
    assert MOLAR_MASS['O'] == 0.01600
