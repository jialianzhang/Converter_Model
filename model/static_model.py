"""转炉静态模型：物料平衡 + 热量平衡"""
import numpy as np
from thermo.constants import MOLAR_MASS, OXIDATION_THERMO

LIME_COMP = {'CaO': 90.0, 'MgO': 3.0, 'SiO2': 2.0, 'Al2O3': 1.0, 'S': 0.05, '烧损': 3.95}
DOLOMITE_COMP = {'CaO': 30.0, 'MgO': 21.0, 'SiO2': 2.0, 'Al2O3': 1.0, '烧损': 46.0}
SINTER_COMP = {'FeO': 9.0, 'Fe2O3': 64.0, 'CaO': 13.0, 'SiO2': 7.0, 'MgO': 2.0, 'Al2O3': 2.0, 'MnO': 1.5, '其他': 1.5}

class StaticModel:
    def __init__(self):
        self.cp_steel = 820.0
        self.cp_slag = 1200.0

    def calculate(self, inputs: dict) -> dict:
        hm = {k: inputs[k] for k in ['w([C])', 'w([Si])', 'w([Mn])', 'w([P])',
                                       'w([S])', 'w([V])', 'w([Ti])']}
        m_hm = inputs['铁水质量'] * 1000
        T_hm = inputs['铁水温度'] + 273.15
        m_scrap = inputs['废钢加入量'] * 1000
        R = inputs['炉渣碱度']

        oxided = {}
        for elem in ['Si', 'Mn', 'P', 'V', 'Ti']:
            oxided[elem] = hm[f'w([{elem}])'] / 100.0 * m_hm

        w_C_target = inputs['终点碳目标']
        dC = hm['w([C])'] - w_C_target
        oxided['C'] = max(0.0, dC / 100.0 * m_hm)
        oxided['Fe'] = 0.01 * m_hm

        O2_kg = 0.0
        for elem, mass_kg in oxided.items():
            n = OXIDATION_THERMO[elem]['n']
            M_elem = MOLAR_MASS[elem]
            moles = mass_kg / M_elem
            O2_kg += moles * n * 0.016

        O2_actual = O2_kg / 0.95 / 1.429

        slag_oxides = {'SiO2': 0.0, 'MnO': 0.0, 'P2O5': 0.0, 'FeO': 0.0, 'CaO': 0.0, 'MgO': 0.0}
        for elem in ['Si', 'Mn', 'P', 'V', 'Ti', 'Fe']:
            mass_kg = oxided.get(elem, 0.0)
            if mass_kg <= 0: continue
            M_elem = MOLAR_MASS[elem]
            if elem == 'Si': slag_oxides['SiO2'] += mass_kg / M_elem * 0.06009
            elif elem == 'Mn': slag_oxides['MnO'] += mass_kg / M_elem * 0.07094
            elif elem == 'P': slag_oxides['P2O5'] += mass_kg / M_elem * 0.14194 / 2.0
            elif elem == 'Fe': slag_oxides['FeO'] += mass_kg / M_elem * 0.07185

        SiO2_in_slag = slag_oxides['SiO2']
        CaO_need = R * SiO2_in_slag
        lime_kg = max(0.0, CaO_need / (LIME_COMP['CaO'] / 100.0))

        MgO_target_pct = inputs['MgO目标含量']
        slag_mass_est = sum(slag_oxides.values()) + lime_kg * (LIME_COMP['CaO'] + LIME_COMP['MgO']) / 100.0
        MgO_need = MgO_target_pct / 100.0 * slag_mass_est
        dolomite_kg = max(0.0, (MgO_need - lime_kg * LIME_COMP['MgO'] / 100.0) / (DOLOMITE_COMP['MgO'] / 100.0))

        m_total_metal = m_hm + m_scrap
        m_oxidized = sum(oxided.values())
        steel_output_t = (m_total_metal - m_oxidized) / 1000.0
        iron_loss = (inputs['铁水质量'] + inputs['废钢加入量']) / steel_output_t * 1000.0

        return {
            '总耗氧量': round(O2_actual, 0),
            '石灰加入量': round(lime_kg, 0),
            '白云石加入量': round(dolomite_kg, 0),
            '合理废钢量': inputs['废钢加入量'],
            '出钢量': round(steel_output_t, 2),
            '铁耗': round(iron_loss, 2),
        }
