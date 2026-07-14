import pytest
import numpy as np
from thermo.equilibrium import calc_lgK, calc_K_M

def test_calc_lgK_C_at_1873K():
    lgK = calc_lgK('C', 1873.0)
    # JSPS标准值: lgK = 1160/1873 + 2.00 ≈ 2.62
    assert 2.0 < lgK < 3.5

def test_calc_lgK_Si_at_1873K():
    lgK = calc_lgK('Si', 1873.0)
    # JSPS标准值: lgK = 30110/1873 - 11.40 ≈ 4.68 (Si强烈氧化, K_Si≈4.8e4)
    # 旧版A值为负导致 lgK<0(Si无法氧化)是错误的, 已在V3修正
    assert 4.0 < lgK < 5.5

def test_calc_K_M_positive():
    K = calc_K_M('Fe', 1873.0)
    assert K > 0

def test_unknown_element_raises():
    with pytest.raises(ValueError):
        calc_lgK('Xx', 1873.0)
