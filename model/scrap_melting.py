"""废钢熔化模型 (张强 式11)"""
import numpy as np


class ScrapMelting:
    def calc_T_MP(self, w_C):
        """废钢熔点 K: T_MP = 1809 - 8.13*w(C)^2 - 54*w(C)"""
        return 1809.0 - 8.13 * w_C**2 - 54.0 * w_C

    def determine_mechanism(self, T_bath, T_MP):
        """判定熔化机制"""
        if T_bath < T_MP:
            return 'diffusion'
        return 'forced'

    def calc_diffusion_rate(self, delta_C, A_scrap, k_mass):
        """扩散熔化速率 (碳浓度梯度驱动, 简化)"""
        return k_mass * A_scrap * delta_C

    def calc_forced_rate(self, delta_T, A_scrap, h_heat):
        """强制熔化速率 (温差驱动, 简化)"""
        return h_heat * A_scrap * delta_T
