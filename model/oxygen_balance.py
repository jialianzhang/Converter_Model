"""界面氧平衡残差 (式6修正): Residual = Sum(n_i*J_i) - J_O"""
import numpy as np


class OxygenBalance:
    def calc_J_M(self, F_M, w_M_b, w_M_star):
        return F_M * (w_M_b - w_M_star)

    def calc_a_O_b_from_slag(self, w_FeO_b, f_FeO, K_Fe):
        if K_Fe < 1e-30:
            return 1e-6
        return max(w_FeO_b * f_FeO / K_Fe, 1e-6)

    def calc_J_O(self, F_O, a_O_b, a_O_star, f_O=1.0):
        w_O_b = a_O_b / f_O if f_O > 0 else a_O_b
        w_O_star = a_O_star / f_O if f_O > 0 else a_O_star
        return F_O * (w_O_b - w_O_star)

    def calc_residual(self, J_fluxes, F_O, a_O_b_slag, a_O_star, f_O=1.0):
        residual = (
            2.0 * J_fluxes.get('Si', 0.0) +
            1.0 * J_fluxes.get('Mn', 0.0) +
            2.5 * J_fluxes.get('P', 0.0) +
            1.5 * J_fluxes.get('V', 0.0) +
            2.0 * J_fluxes.get('Ti', 0.0) +
            1.0 * J_fluxes.get('Fe', 0.0) +
            1.0 * J_fluxes.get('C', 0.0)
        )
        J_O = self.calc_J_O(F_O, a_O_b_slag, a_O_star, f_O)
        return residual - J_O
