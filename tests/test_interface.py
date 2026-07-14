import pytest
import numpy as np
from model.interface import InterfaceConcentration


def test_calc_interface_wM():
    ic = InterfaceConcentration()
    wM_star = ic.calc_w_M_star(F_M=0.5, F_MOn=0.3, w_M_b=4.5, w_MOn_b=1.0, E_M=10000.0, a_O_star=0.01, n=2.0)
    assert wM_star > 0
    assert wM_star < 4.5


def test_calc_all_includes_carbon_special_case():
    ic = InterfaceConcentration()
    F_dict = {'F_C': 1.0, 'F_Si': 0.5, 'F_Mn': 0.3, 'F_P': 0.2, 'F_Fe': 2.0,
              'F_SiO2': 0.3, 'F_MnO': 0.2, 'F_P2O5': 0.15, 'F_FeO': 1.0,
              'F_O': 5.0, 'F_V': 0.1, 'F_Ti': 0.1, 'F_V2O3': 0.05, 'F_TiO2': 0.05}
    steel_b = {'w([C])b': 4.5, 'w([Si])b': 0.3, 'w([Mn])b': 0.19,
               'w([P])b': 0.1, 'w([V])b': 0.0, 'w([Ti])b': 0.0, 'w([Fe])b': 95.0}
    slag_b = {'w((SiO2))b': 15.0, 'w((MnO))b': 5.0, 'w((P2O5))b': 3.0,
              'w((FeO))b': 20.0, 'w((V2O3))b': 0.0, 'w((TiO2))b': 0.0}
    E_M = {'C': 50, 'Si': 100, 'Mn': 80, 'P': 60, 'V': 40, 'Ti': 50, 'Fe': 20}
    result = ic.calc_all(F_dict, steel_b, slag_b, E_M, a_O_star=1e-4, G_CO=0.01, K_C=50.0)
    assert 'w([C])*' in result
    assert len(result) == 6  # C, Si, Mn, P, V, Ti (Fe为溶剂, FeO通量由渣侧在engine中单独控制)
    for v in result.values():
        assert v >= 0
    # CO不会重新渗碳: 界面碳浓度不高于本体
    assert result['w([C])*'] <= steel_b['w([C])b'] + 1e-12


def test_carbon_formula_uses_G_CO():
    ic = InterfaceConcentration()
    F_dict = {'F_C': 1.0, 'F_Si': 0.5, 'F_Mn': 0.3, 'F_P': 0.2, 'F_Fe': 2.0,
              'F_SiO2': 0.3, 'F_MnO': 0.2, 'F_P2O5': 0.15, 'F_FeO': 1.0,
              'F_O': 5.0, 'F_V': 0.1, 'F_Ti': 0.1, 'F_V2O3': 0.05, 'F_TiO2': 0.05}
    steel_b = {'w([C])b': 4.5, 'w([Si])b': 0.3, 'w([Mn])b': 0.19,
               'w([P])b': 0.1, 'w([V])b': 0.0, 'w([Ti])b': 0.0, 'w([Fe])b': 95.0}
    slag_b = {'w((SiO2))b': 15.0, 'w((MnO))b': 5.0, 'w((P2O5))b': 3.0,
              'w((FeO))b': 20.0, 'w((V2O3))b': 0.0, 'w((TiO2))b': 0.0}
    E_M = {'C': 50, 'Si': 100, 'Mn': 80, 'P': 60, 'V': 40, 'Ti': 50, 'Fe': 20}
    # a_O*=0.01 高于该条件下 C/CO 平衡值, 界面浓度未触发钳位, G_CO 才有区分度
    r1 = ic.calc_all(F_dict, steel_b, slag_b, E_M, a_O_star=0.01, G_CO=0.01, K_C=50.0)
    r2 = ic.calc_all(F_dict, steel_b, slag_b, E_M, a_O_star=0.01, G_CO=10.0, K_C=50.0)
    assert r1['w([C])*'] != r2['w([C])*']
