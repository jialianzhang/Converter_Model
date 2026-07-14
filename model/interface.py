"""界面浓度计算: 式7通用 + 碳专属气体界面公式"""
import numpy as np


class InterfaceConcentration:
    def calc_w_M_star(self, F_M, F_MOn, w_M_b, w_MOn_b, E_M, a_O_star, n):
        """w([M])* = (F_MOn*w((MOn))b + F_M*w([M])b) / (F_M + F_MOn*E_M*(a*_O)^n)"""
        numerator = F_MOn * w_MOn_b + F_M * w_M_b
        denominator = F_M + F_MOn * E_M * (a_O_star ** n)
        if denominator < 1e-30:
            return 0.0
        return numerator / denominator

    def calc_all(self, F_dict, steel_b, slag_b, E_M, a_O_star, G_CO, K_C, f_C=1.0):
        """All elements interface concentrations. C uses gas formula, others use slag-metal (式7).

        Fe 不在此处理：Fe 为溶剂(a_Fe=1), 其 FeO 通量由渣侧传质单独控制 (engine Step 4b)。
        碳的 P*_CO = K_C·f_C·w([C])*·a_O*, 并内置 w([C])*<=w([C])b 钳位
        (CO 为气泡逸出, 不会重新渗碳; 钳位内置保证求解残差与通量施加一致)。
        """
        elements_map = {
            'Si': ('F_Si', 'F_SiO2', 'SiO2', 2.0),
            'Mn': ('F_Mn', 'F_MnO', 'MnO', 1.0),
            'P':  ('F_P', 'F_P2O5', 'P2O5', 2.5),
            'V':  ('F_V', 'F_V2O3', 'V2O3', 1.5),
            'Ti': ('F_Ti', 'F_TiO2', 'TiO2', 2.0),
        }
        result = {}

        # Carbon: w([C])* = (F_C*w([C])b + G_CO) / (F_C + G_CO*K_C*f_C*a_O*)
        F_C = F_dict.get('F_C', 0.0)
        w_C_b = steel_b.get('w([C])b', 0.0)
        num = F_C * w_C_b + G_CO
        den = F_C + G_CO * K_C * f_C * a_O_star
        w_C_star = num / den if den > 1e-30 else 0.0
        result['w([C])*'] = min(w_C_star, w_C_b)

        for elem, (fk, fo, ox, n) in elements_map.items():
            FM = F_dict.get(fk, 0.0)
            wMb = steel_b.get(f'w([{elem}])b', 0.0)
            FMO = F_dict.get(fo, 0.0)
            wMOb = slag_b.get(f'w(({ox}))b', 0.0)
            E = E_M.get(elem, 1.0)
            result[f'w([{elem}])*'] = self.calc_w_M_star(FM, FMO, wMb, wMOb, E, a_O_star, n)
        return result
