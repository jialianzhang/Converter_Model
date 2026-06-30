import pytest
import numpy as np
from thermo.equilibrium import calc_lgK, calc_K_M

def test_calc_lgK_C_at_1873K():
    lgK = calc_lgK('C', 1873.0)
    # lgK = -11680/1873 + 9.04 ≈ 2.80
    assert 2.0 < lgK < 4.0

def test_calc_lgK_Si_at_1873K():
    lgK = calc_lgK('Si', 1873.0)
    assert -7.0 < lgK < -3.0

def test_calc_K_M_positive():
    K = calc_K_M('Fe', 1873.0)
    assert K > 0

def test_unknown_element_raises():
    with pytest.raises(ValueError):
        calc_lgK('Xx', 1873.0)
