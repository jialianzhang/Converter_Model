"""石灰熔化模型: Kitamura完整模型 + 张强简化模型"""
import numpy as np


class LimeMelting:
    def calc_k_prime(self, rho_s, eta, D, d, U):
        """Kitamura式11: k' = 0.35 * rho_s^0.36 * eta^(-0.36) * D^0.67 * d^(-0.31) * U^0.67"""
        return 0.35 * (rho_s ** 0.36) * (eta ** (-0.36)) * (D ** 0.67) * (d ** (-0.31)) * (U ** 0.67)

    def calc_dr_dt(self, k_S, rho_S, rho_CaO, w_CaO_sat, w_CaO_slag):
        """张强式9: -dr/dt = k_S * (rho_S/rho_CaO) * (w(CaO)_sat - w(CaO)_slag)/100"""
        return k_S * (rho_S / rho_CaO) * (w_CaO_sat - w_CaO_slag) / 100.0

    def calc_delta_m(self, m_BD, r, delta_r):
        """张强式10: delta_m = m_BD - m_BD * ((r - delta_r)^3 / r^3)"""
        if r < 1e-10:
            return m_BD
        return m_BD - m_BD * ((r - delta_r) ** 3) / (r ** 3)

    def calc_k_prime_C2S(self, k_prime_normal):
        """C2S饱和区: 传质系数降至1/10"""
        return k_prime_normal / 10.0

    def is_fully_melted(self, r):
        """r < 1mm -> 完全熔解"""
        return r < 0.001

    def calc_mass_change_rate(self, n_CaO, A_CaO, k_prime, rho_s, r_prime, delta_CaO_pct, x):
        """Kitamura式12: -dW_CaO/dt = n_CaO * A_CaO * {k'*rho_s/(100*r')} * delta(%CaO) * r' * x"""
        return n_CaO * A_CaO * (k_prime * rho_s / (100.0 * r_prime)) * delta_CaO_pct * r_prime * x
