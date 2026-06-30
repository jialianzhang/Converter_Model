"""主模拟引擎: 13步 Euler-forward 时间步循环"""
import numpy as np
from config import (RHO_STEEL, RHO_SLAG, CP_STEEL, CP_SLAG, KITAMURA_BASE,
                     TIME_STEP_INTENSE, TIME_STEP_NORMAL, TIME_STEP_END, SOLVER_TOL)
from thermo.constants import OXIDATION_THERMO, MOLAR_MASS, OXIDE_MOLAR_MASS
from thermo.equilibrium import calc_K_M, calc_lgK
from thermo.wagner import calc_all_f_i
from thermo.ion_theory import calc_slag_activity_coeff
from thermo.slag_phase import calc_CaO_saturation, calc_solid_fraction
from model.static_model import StaticModel
from model.mixing_energy import MixingEnergy
from model.modified_mass_transfer import ModifiedMassTransfer
from model.interface import InterfaceConcentration
from model.oxygen_balance import OxygenBalance
from model.root_solver import RootSolver
from model.lime_melting import LimeMelting
from model.scrap_melting import ScrapMelting
from model.mass_balance import MassBalance
from model.heat_balance import HeatBalance


class Engine:
    def __init__(self):
        self.static = StaticModel()
        self.mixing = MixingEnergy()
        self.transfer = ModifiedMassTransfer()
        self.interface = InterfaceConcentration()
        self.oxy = OxygenBalance()
        self.solver = RootSolver(tol=SOLVER_TOL)
        self.lime = LimeMelting()
        self.scrap = ScrapMelting()
        self.mass_bal = MassBalance()
        self.heat_bal = HeatBalance()

    def select_time_step(self, w_C):
        if w_C > 0.8:
            return TIME_STEP_INTENSE
        elif w_C > 0.3:
            return TIME_STEP_NORMAL
        else:
            return TIME_STEP_END

    def run(self, inputs, max_time=1500.0):
        # Step 0: Static model
        static_result = self.static.calculate(inputs)

        # Initialize state
        m_hm = inputs['铁水质量'] * 1000.0  # kg
        m_scrap = inputs['废钢加入量'] * 1000.0  # kg
        m_steel = m_hm + m_scrap
        T = inputs['铁水温度'] + 273.15
        lime_total = static_result['石灰加入量']
        dolomite_total = static_result['白云石加入量']

        # Steel composition (mass in kg, not fraction)
        steel_comp = {}
        for elem in ['C', 'Si', 'Mn', 'P', 'S', 'V', 'Ti']:
            steel_comp[elem] = inputs[f'w([{elem}])'] / 100.0 * m_hm
        steel_comp['Fe'] = m_steel - sum(steel_comp.values())

        # Slag composition (initialize from static, wt%)
        w_SiO2_b = 15.0
        w_CaO_b = 45.0
        w_FeO_b = 20.0
        w_MnO_b = 5.0
        w_MgO_b = 8.0
        w_P2O5_b = 3.0
        m_slag = lime_total * 0.9 + dolomite_total * 0.54
        # Track CaO mass separately for additive lime-melting updates
        m_CaO_initial = m_slag * w_CaO_b / 100.0  # kg CaO at t=0

        # Time series
        ts = {'time': [], 'T': [], 'w_C': [], 'w_Si': [], 'w_Mn': [], 'w_P': [],
              'w_FeO': [], 'w_SiO2': [], 'w_CaO': []}

        # Previous step variables
        q_CO_prev = 0.0
        O2_consumed = 0.0
        t = 0.0
        O2_total_target = static_result['总耗氧量']

        while t < max_time:
            # Adaptive time step
            w_C_current = steel_comp.get('C', 0.0) / m_steel * 100.0 if m_steel > 0 else 0.0
            dt = self.select_time_step(w_C_current)

            # Step 1: Mixing energy
            eps_T = self.mixing.calc_epsilon_T(inputs['供氧流量'], m_steel, inputs['枪位高度'],
                                                inputs['喷嘴直径'], 300.0, inputs['喷嘴夹角'])
            eps_B = self.mixing.calc_epsilon_B(inputs['底吹氩气流量'], T, m_steel, inputs['熔池深度'])
            eps_total = eps_T + eps_B

            # Step 2: Mass transfer coefficients (use q_CO_prev)
            eps_total0 = 50.0  # baseline
            km0 = KITAMURA_BASE['k_m0']
            ks0 = KITAMURA_BASE['k_s0']
            km = self.transfer.calc_k_m(km0, eps_total, eps_total0, gamma=0.3)
            q_Ar = inputs['底吹氩气流量'] / 60000.0  # Nl/min -> Nm3/s
            # Cap the gas-flow ratio in ks_prime to prevent lab-scale q_Ar0
            # from blowing up ks when industrial q_CO reaches Nm3/s levels.
            q_gas_ratio = min((q_Ar + 3.0 * q_CO_prev) / KITAMURA_BASE['q_Ar0'], 10.0)
            ks_prime = ks0 * q_gas_ratio * (KITAMURA_BASE['W_s0'] / max(m_slag, 1e-6))
            f_solid = calc_solid_fraction(w_CaO_b, w_SiO2_b, w_FeO_b)
            ks = self.transfer.calc_k_s(ks_prime, f_solid, KITAMURA_BASE['n_solid'])

            # Convert to F coefficients
            F_dict = self.transfer.calc_all_F(km, ks, RHO_STEEL, RHO_SLAG, MOLAR_MASS, OXIDE_MOLAR_MASS)

            # Step 3: Equilibrium constants + E_M
            K_M = {elem: calc_K_M(elem, T) for elem in OXIDATION_THERMO}

            # Build steel composition in wt% for Wagner model
            steel_wt_pct = {}
            wagner_elements = ['C', 'Si', 'Mn', 'P', 'S', 'V', 'Ti', 'O', 'Fe']
            for elem in wagner_elements:
                steel_wt_pct[elem] = steel_comp.get(elem, 0.0) / m_steel * 100.0 if m_steel > 0 else 0.0
            f_steel = calc_all_f_i(steel_wt_pct, T)

            slag_comp_pct = {'FeO': w_FeO_b, 'SiO2': w_SiO2_b, 'CaO': w_CaO_b,
                             'MnO': w_MnO_b, 'MgO': w_MgO_b, 'P2O5': w_P2O5_b}
            f_slag = {}
            for ox in slag_comp_pct:
                try:
                    f_slag[ox] = calc_slag_activity_coeff(ox, slag_comp_pct, T)
                except Exception:
                    f_slag[ox] = 1.0

            C_slag = RHO_SLAG / 0.060  # approximate molar concentration mol/m3

            E_M = {}
            for elem in ['Si', 'Mn', 'P', 'V', 'Ti', 'Fe']:
                d = OXIDATION_THERMO[elem]
                f_M = f_steel.get(elem, 1.0)
                ox_map = {'Si': 'SiO2', 'Mn': 'MnO', 'P': 'P2O5', 'V': 'V2O3', 'Ti': 'TiO2', 'Fe': 'FeO'}
                ox = ox_map[elem]
                f_MOn = f_slag.get(ox, 1.0)
                N_MOn = OXIDE_MOLAR_MASS.get(ox, 0.1)
                E_M[elem] = 100.0 * C_slag * N_MOn * f_M * K_M[elem] / (RHO_SLAG * f_MOn)

            # Step 4: Solve a*_O via Brent
            steel_b = {}
            for elem in ['C', 'Si', 'Mn', 'P', 'V', 'Ti', 'Fe']:
                steel_b[f'w([{elem}])b'] = steel_comp.get(elem, 0.0) / m_steel * 100.0 if m_steel > 0 else 0.0
            slag_b = {'w((SiO2))b': w_SiO2_b, 'w((MnO))b': w_MnO_b, 'w((P2O5))b': w_P2O5_b,
                       'w((FeO))b': w_FeO_b, 'w((V2O3))b': 0.0, 'w((TiO2))b': 0.0}

            K_C = K_M.get('C', 50.0)
            # Industrial-scale G_CO:
            # The Kitamura formula G_CO = G_CO0*(q_gas/q_Ar0)*(W_s0/W_s) was
            # calibrated at 100-kg lab scale.  At 170-t industrial scale the
            # W_s0/W_s factor (~0.0004) suppresses G_CO by ~2500x.  We remove
            # that factor and also replace the lab-scale q_Ar0 normalisation
            # (which causes blow-up when q_CO reaches industrial Nm3/s) with a
            # physically-grounded reference q_CO_ref.
            G_CO0_eff = KITAMURA_BASE['G_CO0'] * (max(m_slag, 1e-6) / KITAMURA_BASE['W_s0'])
            # Reference CO flow rate ~ O2 injection rate [Nm3/s].
            q_CO_ref = max(inputs['供氧流量'] / 60.0, 1e-6)
            G_CO = G_CO0_eff * (q_Ar / KITAMURA_BASE['q_Ar0'] + q_CO_prev / q_CO_ref)
            G_CO = max(G_CO, 1e-6)  # numerical floor
            # Industrial-scale calibration: A0=10 [m2] was calibrated at
            # lab scale (100 kg).  For 170-t BOF the effective reaction
            # interface area is ~100x larger due to larger hot-spot,
            # emulsification, and bath surface.
            A0_eff = KITAMURA_BASE['A0'] * 100.0
            A = A0_eff * (1.0 - f_solid) ** (2.0 / 3.0)

            K_Fe = K_M.get('Fe', 10.0)
            f_FeO_slag = f_slag.get('FeO', 1.0)
            a_O_b_slag = self.oxy.calc_a_O_b_from_slag(w_FeO_b, f_FeO_slag, K_Fe)
            # Cap at physically-reasonable range for dynamic BOF conditions.
            a_O_b_slag = min(a_O_b_slag, 5.0)

            # O2 injection contributes to total O supply (mol O / m2 / s)
            O2_mol_per_s = inputs['供氧流量'] / (60.0 * 0.0224)  # Nm3/min -> mol O2/s
            J_O2_inj = O2_mol_per_s * 2.0 / max(A, 1e-6) * 0.015  # mol O/(m2·s), 1.5% efficiency

            def residual_func(a_O_star):
                w_star = self.interface.calc_all(F_dict, steel_b, slag_b, E_M, a_O_star, G_CO, K_C)
                J_fluxes = {}
                for elem in ['C', 'Si', 'Mn', 'P', 'V', 'Ti', 'Fe']:
                    w_b = steel_b.get(f'w([{elem}])b', 0.0)
                    w_s = w_star.get(f'w([{elem}])*', 0.0)
                    F_M = F_dict.get(f'F_{elem}', 0.0)
                    J_fluxes[elem] = self.oxy.calc_J_M(F_M, w_b, w_s)
                slag_residual = self.oxy.calc_residual(J_fluxes, F_dict.get('F_O', 1.0), a_O_b_slag, a_O_star, f_steel.get('O', 1.0))
                # Subtract O2 injection contribution from the residual
                # (positive slag_residual = O demand > slag supply;
                #  subtract J_O2_inj so both sources together balance demand)
                return slag_residual - J_O2_inj

            try:
                a_O_star = self.solver.solve(residual_func, bracket=(1e-8, 1.0))
            except Exception:
                a_O_star = 1e-6  # fallback

            # Steps 5-6: Interface concentrations and reaction rates
            w_star = self.interface.calc_all(F_dict, steel_b, slag_b, E_M, a_O_star, G_CO, K_C)

            delta_m = {}
            for elem in ['C', 'Si', 'Mn', 'P', 'V', 'Ti', 'Fe']:
                w_b = steel_b.get(f'w([{elem}])b', 0.0)
                w_s = w_star.get(f'w([{elem}])*', 0.0)
                F_M = F_dict.get(f'F_{elem}', 0.0)
                M_i = MOLAR_MASS.get(elem, 0.05)
                delta_m[elem] = A * F_M * M_i * (w_b - w_s) * dt

            # Track q_CO for next step
            if 'C' in delta_m and dt > 0:
                q_CO_prev = delta_m['C'] / MOLAR_MASS['C'] / dt * 0.0224  # mol/s -> Nm3/s

            # Steps 7-8: Mass balance
            for elem in ['C', 'Si', 'Mn', 'P', 'V', 'Ti']:
                if elem in delta_m and elem in steel_comp:
                    steel_comp[elem] = max(0.0, steel_comp[elem] - delta_m[elem])
            steel_comp['Fe'] = max(0.0, steel_comp.get('Fe', 0.0) - delta_m.get('Fe', 0.0))
            m_steel = sum(steel_comp.values())

            # Simplified slag update
            if 'Si' in delta_m:
                w_SiO2_b += delta_m['Si'] / MOLAR_MASS['Si'] * 0.06009 / max(m_slag, 1e-6) * 100.0
            if 'Fe' in delta_m:
                w_FeO_b += delta_m['Fe'] / MOLAR_MASS['Fe'] * 0.07185 / max(m_slag, 1e-6) * 100.0
            m_slag += sum(delta_m.values()) * 1.5

            # Steps 9-10: Lime/scrap melting (simplified)
            # Track incremental lime melted this step for heat balance
            lime_melted_prev = min(lime_total, lime_total * max(0.0, t - dt) / 900.0)
            lime_melted = min(lime_total, lime_total * t / 900.0)
            delta_lime = lime_melted - lime_melted_prev  # kg melted this step
            # CaO from lime (90 %) + dolomite (30 %) enters slag additively
            m_CaO_from_add = lime_melted * 0.9 + dolomite_total * 0.3 * min(1.0, t / 900.0)
            w_CaO_b = (m_CaO_initial + m_CaO_from_add) / max(m_slag, 1e-6) * 100.0
            w_CaO_b = min(w_CaO_b, 65.0)

            # Scrap melts over time (track incremental)
            T_MP_scrap = self.scrap.calc_T_MP(0.2)
            if T > T_MP_scrap:
                m_scrap_melted_prev = min(m_scrap, m_scrap * max(0.0, t - dt) / 600.0)
                m_scrap_melted = min(m_scrap, m_scrap * t / 600.0)
            else:
                m_scrap_melted_prev = m_scrap * max(0.0, t - dt) / 1200.0 * 0.5
                m_scrap_melted = m_scrap * t / 1200.0 * 0.5
            delta_scrap = m_scrap_melted - m_scrap_melted_prev  # kg melted this step
            steel_comp['Fe'] += delta_scrap  # incremental Fe from scrap
            m_steel = sum(steel_comp.values())

            # Steps 11-12: Heat balance
            oxidation_heat = {}
            for elem in ['C', 'Si', 'Mn', 'P', 'V', 'Ti', 'Fe']:
                if elem in delta_m:
                    moles = delta_m[elem] / MOLAR_MASS.get(elem, 0.05)
                    oxidation_heat[elem] = moles * OXIDATION_THERMO[elem]['delta_H']

            Q_HG = self.heat_bal.calc_Q_HG(oxidation_heat, post_comb_ratio=0.1)
            heat_losses = {'radiation': 0.03 * Q_HG, 'offgas': 0.10 * Q_HG,
                           'lining': 0.02 * Q_HG, 'dust': 0.01 * Q_HG}
            Q_HC = self.heat_bal.calc_Q_HC(heat_losses)
            # Lime decomposition: endothermic, proportional to incremental lime melted
            Q_CM = delta_lime * 1780000.0
            # Scrap melting: endothermic, ~250 kJ/kg
            Q_scrap = delta_scrap * 250000.0

            delta_T = self.heat_bal.calc_delta_T(Q_HG, Q_HC, Q_CM + Q_scrap, CP_STEEL, m_steel, CP_SLAG, m_slag)
            T += max(-50.0, min(delta_T, 100.0))  # clamp extreme changes

            # Record
            ts['time'].append(t)
            ts['T'].append(T)
            w_C_pct = steel_comp.get('C', 0.0) / m_steel * 100.0 if m_steel > 0 else 0.0
            ts['w_C'].append(w_C_pct)
            ts['w_Si'].append(steel_comp.get('Si', 0.0) / m_steel * 100.0 if m_steel > 0 else 0.0)
            ts['w_Mn'].append(steel_comp.get('Mn', 0.0) / m_steel * 100.0 if m_steel > 0 else 0.0)
            ts['w_P'].append(steel_comp.get('P', 0.0) / m_steel * 100.0 if m_steel > 0 else 0.0)
            ts['w_FeO'].append(w_FeO_b)
            ts['w_SiO2'].append(w_SiO2_b)
            ts['w_CaO'].append(w_CaO_b)

            t += dt
            O2_consumed += inputs['供氧流量'] * dt / 60.0

            # Composite termination: min O2 consumed AND C <= target
            if O2_consumed >= O2_total_target * 0.5 and w_C_pct <= inputs['终点碳目标']:
                break

        return {
            '终点C': round(ts['w_C'][-1], 4) if ts['w_C'] else 0.0,
            '终点T': round(ts['T'][-1] - 273.15, 1) if ts['T'] else 0.0,
            '时序数据': ts,
            '步数': len(ts['time']),
        }
