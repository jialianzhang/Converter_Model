"""物质平衡更新 (式12): m_i^(l) = m_i^(l-1) - delta_m_i + sum(delta_m_i^CM)"""
import numpy as np


class MassBalance:
    def update(self, m_prev, delta_m_reaction, delta_m_melting):
        """更新各组分质量 (kg)"""
        m_new = {}
        all_keys = set(list(m_prev.keys()) + list(delta_m_reaction.keys()) + list(delta_m_melting.keys()))
        for k in all_keys:
            prev = m_prev.get(k, 0.0)
            react = delta_m_reaction.get(k, 0.0)
            melt = delta_m_melting.get(k, 0.0)
            m_new[k] = max(0.0, prev - react + melt)
        return m_new

    def update_steel_and_slag(self, steel_mass, slag_mass, steel_reactions, slag_reactions,
                               steel_melting, slag_melting):
        """同时更新钢液和炉渣质量"""
        steel = self.update(steel_mass, steel_reactions, steel_melting)
        slag = self.update(slag_mass, slag_reactions, slag_melting)
        return steel, slag
