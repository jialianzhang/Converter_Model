"""混合能计算: Nakanishi 公式 (式14, 式15)"""
import numpy as np

class MixingEnergy:
    def calc_epsilon_T(self, Q_T: float, W: float, X: float, D: float, U: float, theta: float) -> float:
        return 0.0403 * Q_T * D * U**2 * (np.cos(np.radians(theta)))**2 / (W * X)

    def calc_epsilon_B(self, Q_B: float, T: float, W: float, H: float) -> float:
        return (28.5 * Q_B * T / W) * np.log10(1.0 + H / 148.0)

    def calc_total(self, Q_T, Q_B, T, W, X, D, U, theta, H):
        return self.calc_epsilon_T(Q_T, W, X, D, U, theta) + self.calc_epsilon_B(Q_B, T, W, H)
