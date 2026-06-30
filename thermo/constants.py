"""热力学常数：元素氧化反应数据、摩尔质量、氧化放热量"""

# 摩尔质量 kg/mol
MOLAR_MASS = {
    'C':  0.01201, 'Si': 0.02809, 'Mn': 0.05494,
    'P':  0.03097, 'S':  0.03207, 'V':  0.05094,
    'Ti': 0.04787, 'Fe': 0.05585, 'O':  0.01600,
    'Ca': 0.04008, 'Mg': 0.02431, 'Al': 0.02698,
}

# 氧化物摩尔质量 kg/mol
OXIDE_MOLAR_MASS = {
    'SiO2':  0.06009, 'MnO':   0.07094, 'P2O5': 0.14194,
    'V2O3':  0.14988, 'TiO2':  0.07987, 'FeO':  0.07185,
    'CaO':   0.05608, 'MgO':   0.04031, 'Al2O3': 0.10196,
    'Fe2O3': 0.15969,
}

# 各元素氧化反应热力学数据
# n: [M] + n[O] = (MO_n) 中的化学计量系数
# delta_H: 氧化反应放热量 J/mol (放热为正)
# lgK = A/T + B (平衡常数)
OXIDATION_THERMO = {
    'C':   {'n': 1,   'A': -11680.0, 'B':  9.04,  'delta_H':  140000.0},
    'Si':  {'n': 2,   'A': -30150.0, 'B': 11.52,  'delta_H':  580000.0},
    'Mn':  {'n': 1,   'A': -12720.0, 'B':  5.65,  'delta_H':  280000.0},
    'P':   {'n': 2.5, 'A':  -8420.0, 'B':  0.94,  'delta_H':  750000.0},
    'V':   {'n': 1.5, 'A': -17800.0, 'B':  8.21,  'delta_H':  450000.0},
    'Ti':  {'n': 2,   'A': -32200.0, 'B': 12.86,  'delta_H':  680000.0},
    'Fe':  {'n': 1,   'A':  -6340.0, 'B':  2.75,  'delta_H':  120000.0},
}

# C₂S 和 C₃P 形成热 J/mol
HEAT_C2S_FORMATION = 120000.0   # 2CaO + SiO₂ → 2CaO·SiO₂
HEAT_C3P_FORMATION = 200000.0   # 3CaO + P₂O₅ → 3CaO·P₂O₅

# 二次燃烧热 CO + 1/2 O₂ → CO₂  J/mol CO
HEAT_POST_COMBUSTION = 283000.0

# 辅料分解热 J/kg
HEAT_LIME_DECOMPOSITION = 1780000.0      # CaCO₃ → CaO + CO₂
HEAT_DOLOMITE_DECOMPOSITION = 3020000.0   # CaMg(CO₃)₂ → CaO + MgO + 2CO₂
