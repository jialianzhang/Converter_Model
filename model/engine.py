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
        if w_C > 0.8: return TIME_STEP_INTENSE
        elif w_C > 0.3: return TIME_STEP_NORMAL
        else: return TIME_STEP_END

    def run(self, inputs, max_time=1500.0):
        static_result = self.static.calculate(inputs)

        # --- State initialization ---
        m_hm = inputs['铁水质量'] * 1000.0
        m_scrap_init = inputs['废钢加入量'] * 1000.0
        m_steel = m_hm + m_scrap_init
        T = inputs['铁水温度'] + 273.15

        steel_kg = {}
        for elem in ['C','Si','Mn','P','S','V','Ti']:
            steel_kg[elem] = inputs[f'w([{elem}])'] / 100.0 * m_hm
        steel_kg['Fe'] = m_steel - sum(steel_kg.values())
        m_scrap_remaining = m_scrap_init

        lime_total = static_result['石灰加入量']
        dolomite_total = static_result['白云石加入量']
        # Initial slag: hot metal carry-over (~1% of HM mass)
        m_slag_initial = m_hm * 0.008
        slag_kg = {
            'SiO2': m_slag_initial*0.15, 'CaO': m_slag_initial*0.40,
            'FeO':  m_slag_initial*0.20, 'MnO': m_slag_initial*0.08,
            'MgO':  m_slag_initial*0.10, 'P2O5':m_slag_initial*0.07,
        }
        lime_remaining = lime_total
        dolo_remaining = dolomite_total
        sinter_total = inputs['铁水质量'] * 18.0
        sinter_remaining = sinter_total

        ts = {'time_s':[],'T_C':[],'T_K':[],
              'w_C':[],'w_Si':[],'w_Mn':[],'w_P':[],'w_S':[],'w_Fe':[],
              'w_FeO':[],'w_SiO2':[],'w_CaO':[],'w_MgO':[],'w_MnO':[],'w_P2O5':[],
              'm_steel_kg':[],'m_slag_kg':[]}

        q_CO_prev = 0.0;  O2_consumed = 0.0;  t = 0.0
        self._sinter_deferred = 0.0
        O2_total_target = static_result['总耗氧量']

        while t < max_time:
            w_C_pct = steel_kg.get('C',0)/m_steel*100 if m_steel>0 else 0
            dt = self.select_time_step(w_C_pct)

            m_slag = max(sum(slag_kg.values()), 1.0)
            slag_wt = {k:v/m_slag*100 for k,v in slag_kg.items()}
            f_solid = calc_solid_fraction(slag_wt.get('CaO',45),slag_wt.get('SiO2',15),slag_wt.get('FeO',20))
            steel_wt_pct = {e:steel_kg.get(e,0)/m_steel*100 if m_steel>0 else 0
                           for e in ['C','Si','Mn','P','S','V','Ti','Fe']}

            # === Steps 1-2: Mixing energy + dynamic mass transfer ===
            eps_T = self.mixing.calc_epsilon_T(inputs['供氧流量'],m_steel,inputs['枪位高度'],
                                                inputs['喷嘴直径'],300.0,inputs['喷嘴夹角'])
            eps_B = self.mixing.calc_epsilon_B(inputs['底吹氩气流量'],T,m_steel,inputs['熔池深度'])
            eps_total = eps_T + eps_B
            if t == 0.0: self._eps_total0 = max(eps_total, 0.01)

            km = self.transfer.calc_k_m(KITAMURA_BASE['k_m0'], eps_total, self._eps_total0, gamma=0.3)
            q_Ar = inputs['底吹氩气流量'] / 60000.0
            ks_stir = min(1.0 + 3.0*q_CO_prev/max(q_Ar,1e-10), 10.0)
            ks_prime = KITAMURA_BASE['k_s0'] * ks_stir * (40.0 / max(m_slag,1e-6))
            ks = self.transfer.calc_k_s(ks_prime, f_solid, KITAMURA_BASE['n_solid'])
            F_dict = self.transfer.calc_all_F(km, ks, RHO_STEEL, RHO_SLAG, MOLAR_MASS, OXIDE_MOLAR_MASS)

            # === Step 3: Thermodynamics (E_M + activity) ===
            K_M = {elem: calc_K_M(elem, T) for elem in OXIDATION_THERMO}
            f_steel = calc_all_f_i(steel_wt_pct, T)
            slag_pct = {k:slag_wt[k] for k in ['FeO','SiO2','CaO','MnO','MgO','P2O5'] if k in slag_wt}
            f_slag = {}
            for ox in slag_pct:
                try: f_slag[ox] = calc_slag_activity_coeff(ox, slag_pct, T)
                except: f_slag[ox] = 1.0

            C_slag = RHO_SLAG / 0.060
            E_M = {}
            ox_map = {'Si':'SiO2','Mn':'MnO','P':'P2O5','V':'V2O3','Ti':'TiO2','Fe':'FeO'}
            for elem in ['Si','Mn','P','V','Ti','Fe']:
                ox = ox_map[elem]
                fM = f_steel.get(elem,1.0);  fMOn = f_slag.get(ox,1.0)
                N_MOn = OXIDE_MOLAR_MASS.get(ox,0.1)
                E_M[elem] = 100.0 * C_slag * N_MOn * fM * K_M[elem] / (RHO_SLAG * max(fMOn,1e-6))

            O2_mol_per_s = inputs['供氧流量'] / (60.0 * 0.0224)
            O2_total_mol = O2_mol_per_s * 2.0 * dt

            # === Step 4a: Impact zone — direct oxidation with dynamic C share ===
            # Si and Mn get fixed baseline O2 while present (high O affinity).
            # As Si/Mn deplete, their O2 share transfers to C → natural 3-stage decarb.
            # Fe gets fixed baseline for FeO generation (feeds slag-metal reactions).
            direct_oxidized = {}

            si_base = 0.08; mn_base = 0.02; c_base = 0.48; fe_base = 0.03

            # Si: full allocation while Si > 0.01 kg, else zero → complete oxidation
            si_frac = si_base if steel_kg.get('Si', 0) > 0.01 else 0.0
            o2_Si = O2_total_mol * si_frac
            if steel_kg.get('Si', 0) > 0.001:
                max_Si_mol = steel_kg['Si'] / MOLAR_MASS['Si']
                si_mol = min(max_Si_mol, o2_Si / 2.0)
                dm = si_mol * MOLAR_MASS['Si']
                steel_kg['Si'] -= dm
                slag_kg['SiO2'] += dm * 0.06009 / 0.02809
                direct_oxidized['Si'] = dm

            # Mn: full allocation while Mn > 0.01 kg, else zero
            mn_frac = mn_base if steel_kg.get('Mn', 0) > 0.01 else 0.0
            o2_Mn = O2_total_mol * mn_frac
            if steel_kg.get('Mn', 0) > 0.001:
                max_Mn_mol = steel_kg['Mn'] / MOLAR_MASS['Mn']
                mn_mol = min(max_Mn_mol, o2_Mn)
                dm = mn_mol * MOLAR_MASS['Mn']
                steel_kg['Mn'] -= dm
                slag_kg['MnO'] += dm * 0.07094 / 0.05494
                direct_oxidized['Mn'] = dm

            # C: baseline + O2 released from depleting Si and Mn → 3-stage pattern
            c_frac = c_base + (si_base - si_frac) + (mn_base - mn_frac)
            o2_C = O2_total_mol * c_frac
            if steel_kg.get('C', 0) > 0.001:
                max_C_mol = steel_kg['C'] / MOLAR_MASS['C']
                c_mol = min(max_C_mol, o2_C)
                dm = c_mol * MOLAR_MASS['C']
                steel_kg['C'] -= dm
                direct_oxidized['C'] = dm

            # Fe: fixed FeO generation
            o2_Fe = O2_total_mol * fe_base
            if steel_kg.get('Fe', 0) > 0.001:
                max_Fe_mol = steel_kg['Fe'] / MOLAR_MASS['Fe']
                fe_mol = min(max_Fe_mol, o2_Fe)
                dm = fe_mol * MOLAR_MASS['Fe']
                steel_kg['Fe'] -= dm
                slag_kg['FeO'] += dm * 0.07185 / 0.05585
                direct_oxidized['Fe'] = dm

            m_steel = sum(steel_kg.values())
            steel_wt_pct = {e:steel_kg.get(e,0)/m_steel*100 if m_steel>0 else 0
                           for e in ['C','Si','Mn','P','S','V','Ti','Fe']}
            m_slag = max(sum(slag_kg.values()), 1.0)
            slag_wt = {k:v/m_slag*100 for k,v in slag_kg.items()}

            # === Step 4b: Coupled reaction solver (slag-metal interface, all elements) ===
            # C, Si, Mn, P, Fe all compete for oxygen at the slag-metal interface.
            # Oxygen sources: (a) FeO dissociation in slag → J_O, (b) dissolved [O]
            # diffusing from impact zone → J_O2_inj.
            # The oxygen balance Σ(n_i·J_i) = J_O + J_O2_inj determines a_O*.
            steel_b = {f'w([{e}])b':steel_wt_pct.get(e,0) for e in ['C','Si','Mn','P','V','Ti','Fe']}
            slag_b = {f'w(({o}))b':slag_wt.get(o,0) for o in ['SiO2','MnO','P2O5','FeO','V2O3','TiO2']}

            K_C = K_M.get('C', 50.0)
            G_CO = KITAMURA_BASE['G_CO0'] * (1.0 + 3.0*q_CO_prev/max(q_Ar,1e-10))
            G_CO = max(G_CO, 1e-6)

            A = KITAMURA_BASE['A0'] * 75.0 * (1.0 - f_solid)**(2.0/3.0)

            # a_O_b from slag FeO (equilibrium, mole fraction basis)
            x_FeO_slag = slag_wt.get('FeO', 20) / 100.0 * 0.060 / OXIDE_MOLAR_MASS['FeO']
            a_O_b_slag = x_FeO_slag * f_slag.get('FeO', 1.0) / max(K_M.get('Fe', 4.0), 1e-30)
            a_O_b_slag = max(a_O_b_slag, 1e-10)

            # Dynamic dissolved [O] flux from impact zone to slag-metal interface.
            # More dissolved O reaches the interface when C is low (less O consumed
            # by direct decarburization along the diffusion path).
            o2_interface_frac = 0.03 * (1.0 - min(w_C_pct / 5.0, 0.9))
            J_O2_inj = O2_mol_per_s * 2.0 / max(A, 1e-6) * o2_interface_frac

            def residual_func(a_O_star):
                w_star = self.interface.calc_all(F_dict, steel_b, slag_b, E_M, a_O_star, G_CO, K_C)
                J_fluxes = {}
                for e in ['C','Si','Mn','P','Fe']:
                    F_M = F_dict.get(f'F_{e}', 0)
                    J_fluxes[e] = self.oxy.calc_J_M(F_M, steel_b.get(f'w([{e}])b', 0), w_star.get(f'w([{e}])*', 0))
                J_fluxes['V'] = 0.0
                J_fluxes['Ti'] = 0.0
                return self.oxy.calc_residual(J_fluxes, F_dict.get('F_O', 1.0), a_O_b_slag, a_O_star, f_steel.get('O', 1.0)) - J_O2_inj

            try: a_O_star = self.solver.solve(residual_func, bracket=(1e-10, 0.5))
            except: a_O_star = 1e-8

            # === Steps 5-8: Reaction rates + mass balance ===
            w_star_all = self.interface.calc_all(F_dict, steel_b, slag_b, E_M, a_O_star, G_CO, K_C)

            delta_kg = {}
            for elem in ['C','Si','Mn','P','Fe']:
                w_b = steel_wt_pct.get(elem, 0)
                w_s = w_star_all.get(f'w([{elem}])*', 0)
                if elem == 'C': w_s = min(w_s, w_b)  # CO→[C] impossible in BOF
                F_M = F_dict.get(f'F_{elem}', 0)
                M_i = MOLAR_MASS.get(elem, 0.05)
                delta_kg[elem] = A * F_M * M_i * (w_b - w_s) * dt

            # Steel mass balance (all elements from coupled model)
            for elem in ['C','Si','Mn','P']:
                steel_kg[elem] = max(0.0, steel_kg.get(elem, 0) - delta_kg.get(elem, 0))
            steel_kg['Fe'] = max(0.0, steel_kg.get('Fe', 0) - delta_kg.get('Fe', 0))
            m_steel = sum(steel_kg.values())

            # Slag mass balance (bidirectional: allows FeO/MnO reduction when J<0)
            elem_to_oxide = {
                'Si':('SiO2',0.06009/0.02809), 'Mn':('MnO',0.07094/0.05494),
                'P':('P2O5',0.14194/0.03097/2.0), 'Fe':('FeO',0.07185/0.05585),
            }
            for elem,(ox,factor) in elem_to_oxide.items():
                dm = delta_kg.get(elem,0)
                slag_kg[ox] = max(0.0, slag_kg.get(ox,0) + dm*factor)

            C_mol_coupled = delta_kg.get('C',0)/MOLAR_MASS['C']
            C_mol_direct = direct_oxidized.get('C',0)/MOLAR_MASS['C']
            q_CO_prev = (C_mol_coupled + C_mol_direct)/max(dt,1e-6)*0.0224 if dt>0 else 0

            # === Steps 9-10: Flux + scrap melting ===
            O2_progress = O2_consumed/max(O2_total_target,1.0)
            flux_frac = max(0.0, min(1.0, (O2_progress-0.10)/0.80))

            # Lime (Kitamura dissolution)
            lime_target = flux_frac * lime_total
            dl_target = max(0.0, lime_target - (lime_total-lime_remaining))
            eta_slag = 0.5*np.exp(8000.0/T)*(1.0+10.0*f_solid)
            D_CaO = 1e-9*np.exp(-5000.0/T)
            w_CaO_sat_lim = calc_CaO_saturation(slag_wt.get('FeO',20), slag_wt.get('SiO2',15))
            if lime_remaining>0.01 and slag_wt.get('CaO',45)<w_CaO_sat_lim and dl_target>0:
                kp = self.lime.calc_k_prime(RHO_SLAG, eta_slag, D_CaO, d=0.01, U=0.5)
                dc = w_CaO_sat_lim-slag_wt.get('CaO',45)
                nl = lime_remaining/(3300.0*4/3*np.pi*0.005**3)
                Al = nl*4*np.pi*0.005**2
                dl = self.lime.calc_mass_change_rate(nl,Al,kp,RHO_SLAG,3300.0,dc,x=1.0)*dt
                dl = min(dl, lime_remaining, dl_target)
            else: dl = 0.0
            lime_remaining -= dl
            slag_kg['CaO']+=dl*0.90; slag_kg['MgO']+=dl*0.03; slag_kg['SiO2']+=dl*0.02

            # Dolomite (O2-driven)
            dolo_target = flux_frac*dolomite_total
            dd_target = max(0.0, dolo_target-(dolomite_total-dolo_remaining))
            dd = min(dd_target*dt/60.0, dolo_remaining) if dolo_remaining>0.01 and dd_target>0 else 0.0
            dolo_remaining -= dd
            slag_kg['CaO']+=dd*0.30; slag_kg['MgO']+=dd*0.21; slag_kg['SiO2']+=dd*0.02

            # Sinter (O2-driven, mass enters early for slag chemistry)
            sinter_target = flux_frac*sinter_total
            ds_target = max(0.0, sinter_target-(sinter_total-sinter_remaining))
            dsinter = min(ds_target*dt/60.0, sinter_remaining) if sinter_remaining>0.01 and ds_target>0 else 0.0
            sinter_remaining -= dsinter
            slag_kg['FeO']+=dsinter*(0.09+0.64*0.9); slag_kg['CaO']+=dsinter*0.13
            slag_kg['SiO2']+=dsinter*0.07; slag_kg['MgO']+=dsinter*0.02; slag_kg['MnO']+=dsinter*0.015

            # Scrap melting (parabolic, Zhang Eq.7-2)
            t_end_scrap = O2_total_target/max(inputs['供氧流量'],1.0)*60.0*0.80
            if t<t_end_scrap and m_scrap_remaining>0.01:
                frac=t/t_end_scrap; ds=m_scrap_init*6.0*frac*(1.0-frac)/t_end_scrap*dt
                ds=min(ds, m_scrap_remaining)
            else: ds=min(m_scrap_remaining, m_scrap_init/t_end_scrap*dt)
            m_scrap_remaining-=ds; steel_kg['Fe']+=ds; m_steel=sum(steel_kg.values())

            # === Steps 11-12: Heat balance ===
            ox_heat = {}
            for elem in ['C','Si','Mn','P','Fe']:
                dm_coupled = delta_kg.get(elem, 0)
                dm_direct = direct_oxidized.get(elem, 0)
                dm_total = abs(dm_coupled) + abs(dm_direct)
                if dm_total > 0:
                    ox_heat[elem] = dm_total / MOLAR_MASS[elem] * OXIDATION_THERMO[elem]['delta_H']
            Q_oxid = sum(ox_heat.values())
            Q_post = ox_heat.get('C',0)*0.10
            dm_Si = delta_kg.get('Si',0) + direct_oxidized.get('Si', 0)
            Q_C2S = dm_Si/MOLAR_MASS['Si']*120000.0 if dm_Si>0 else 0.0
            Q_HG = Q_oxid+Q_post+Q_C2S

            Q_scrap = ds*1300000.0
            Q_CM = dl*1780000.0+dd*3020000.0
            Q_slag_heat = 0.0
            if T>1100.0: Q_slag_heat += (dl*0.95+dd*0.53)*CP_SLAG*(T-1100.0)
            if T>298.0: Q_slag_heat += dsinter*1000.0*(T-298.0)
            Q_offgas = q_CO_prev/0.0224*33.0*(T-298.0)*dt if q_CO_prev>0 and T>298 else 0.0
            Q_rad_lining = 2000000.0*dt
            # Sinter heat: deferred to mid-blow when C oxidation provides enough heat
            Q_sinter_now = dsinter*1800000.0 if w_C_pct<2.5 else 0.0
            self._sinter_deferred += dsinter*1800000.0 - Q_sinter_now
            if self._sinter_deferred>0 and w_C_pct<1.0:
                release = min(self._sinter_deferred, self._sinter_deferred*dt/90.0)
                Q_sinter_now += release; self._sinter_deferred -= release

            Q_consumed = Q_scrap+Q_CM+Q_slag_heat+Q_offgas+Q_rad_lining+Q_sinter_now
            delta_T = self.heat_bal.calc_delta_T(max(0,Q_HG),max(0,Q_consumed),0,CP_STEEL,m_steel,CP_SLAG,m_slag)
            T += max(-50.0, min(delta_T, 100.0))

            # === Record ===
            m_slag = sum(slag_kg.values())
            ts['time_s'].append(t); ts['T_C'].append(T-273.15); ts['T_K'].append(T)
            ts['w_C'].append(steel_kg.get('C',0)/m_steel*100 if m_steel>0 else 0)
            ts['w_Si'].append(steel_kg.get('Si',0)/m_steel*100 if m_steel>0 else 0)
            ts['w_Mn'].append(steel_kg.get('Mn',0)/m_steel*100 if m_steel>0 else 0)
            ts['w_P'].append(steel_kg.get('P',0)/m_steel*100 if m_steel>0 else 0)
            ts['w_S'].append(steel_kg.get('S',0)/m_steel*100 if m_steel>0 else 0)
            ts['w_Fe'].append(steel_kg.get('Fe',0)/m_steel*100 if m_steel>0 else 0)
            ts['w_FeO'].append(slag_kg.get('FeO',0)/m_slag*100 if m_slag>0 else 0)
            ts['w_SiO2'].append(slag_kg.get('SiO2',0)/m_slag*100 if m_slag>0 else 0)
            ts['w_CaO'].append(slag_kg.get('CaO',0)/m_slag*100 if m_slag>0 else 0)
            ts['w_MgO'].append(slag_kg.get('MgO',0)/m_slag*100 if m_slag>0 else 0)
            ts['w_MnO'].append(slag_kg.get('MnO',0)/m_slag*100 if m_slag>0 else 0)
            ts['w_P2O5'].append(slag_kg.get('P2O5',0)/m_slag*100 if m_slag>0 else 0)
            ts['m_steel_kg'].append(m_steel); ts['m_slag_kg'].append(m_slag)

            t += dt
            O2_consumed += inputs['供氧流量']*dt/60.0
            w_C_pct = steel_kg.get('C',0)/m_steel*100 if m_steel>0 else 0

            if O2_consumed>=O2_total_target*0.5 and w_C_pct<=inputs['终点碳目标']: break

        return {
            '终点C': round(ts['w_C'][-1],4) if ts['w_C'] else 0.0,
            '终点T': round(ts['T_C'][-1],1) if ts['T_C'] else 0.0,
            '时序数据': ts,
            '步数': len(ts['time_s']),
        }
