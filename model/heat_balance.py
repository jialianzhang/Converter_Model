"""热量平衡与温度更新 (式13): delta_T/delta_t = 1/(c_M*m_M + c_S*m_S) * (Q_HG - Q_HC - Q_CM)/delta_t"""
import numpy as np


class HeatBalance:
    def calc_Q_HG(self, oxidation_heats, post_comb_ratio=0.1):
        """热收入: sum(element oxidation heat) + post combustion"""
        Q = sum(oxidation_heats.values())
        Q_post = oxidation_heats.get('C', 0.0) * post_comb_ratio
        return Q + Q_post

    def calc_Q_HC(self, heat_losses):
        """热支出: sum of all heat losses"""
        return sum(heat_losses.values())

    def calc_delta_T(self, Q_HG, Q_HC, Q_CM, c_M, m_M, c_S, m_S):
        """温度变化 K: delta_T = (Q_HG - Q_HC - Q_CM) / (c_M*m_M + c_S*m_S)"""
        denominator = c_M * m_M + c_S * m_S
        if denominator < 1e-10:
            return 0.0
        return (Q_HG - Q_HC - Q_CM) / denominator

    def calc_temp_change_rate(self, Q_HG, Q_HC, Q_CM, c_M, m_M, c_S, m_S, dt):
        """温度变化率 K/s (式13完整形式)"""
        denominator = c_M * m_M + c_S * m_S
        if denominator < 1e-10:
            return 0.0
        return (Q_HG - Q_HC - Q_CM) / denominator / dt
