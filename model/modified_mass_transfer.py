"""动态修正传质系数: k_m, k_s -> F_M, F_MOn"""
import numpy as np
from config import K_S_MIN

class ModifiedMassTransfer:
    def calc_k_m(self, km0: float, eps_total: float, eps_total0: float, gamma: float) -> float:
        return km0 * (eps_total / eps_total0) ** gamma

    def calc_k_s_prime(self, ks0: float, q_Ar: float, q_CO: float, q_Ar0: float, W_s: float, W_s0: float) -> float:
        return ks0 * ((q_Ar + 3.0 * q_CO) / q_Ar0) * (W_s0 / max(W_s, 1e-6))

    def calc_k_s(self, ks_prime: float, f_solid: float, n: float = 5.0) -> float:
        ks = ks_prime * 10.0 ** (-n * f_solid)
        return max(ks, K_S_MIN)

    def calc_F_M(self, km: float, rho_m: float, M_M: float) -> float:
        return km * rho_m / (100.0 * M_M)

    def calc_F_MOn(self, ks: float, rho_s: float, M_MOn: float) -> float:
        return ks * rho_s / (100.0 * M_MOn)

    def calc_all_F(self, km: float, ks: float, rho_m: float, rho_s: float,
                   molar_masses: dict, oxide_molar_masses: dict) -> dict:
        elements = ['C', 'Si', 'Mn', 'P', 'V', 'Ti', 'Fe']
        oxides = ['SiO2', 'MnO', 'P2O5', 'V2O3', 'TiO2', 'FeO']
        result = {'F_O': self.calc_F_M(km, rho_m, 0.016)}
        for elem in elements:
            if elem in molar_masses:
                result[f'F_{elem}'] = self.calc_F_M(km, rho_m, molar_masses[elem])
        for ox in oxides:
            if ox in oxide_molar_masses:
                result[f'F_{ox}'] = self.calc_F_MOn(ks, rho_s, oxide_molar_masses[ox])
        return result
