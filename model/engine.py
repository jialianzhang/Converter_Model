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
        """Execute the 13-step Euler-forward simulation."""
        # Step 0: Static model
        static_result = self.static.calculate(inputs)

        # --- Initialize state ---
        m_hm = inputs['铁水质量'] * 1000.0
        m_scrap_init = inputs['废钢加入量'] * 1000.0
        m_steel = m_hm + m_scrap_init
        T = inputs['铁水温度'] + 273.15

        # Steel composition (kg per element)
        steel_kg = {}
        for elem in ['C','Si','Mn','P','S','V','Ti']:
            steel_kg[elem] = inputs[f'w([{elem}])'] / 100.0 * m_hm
        steel_kg['Fe'] = m_steel - sum(steel_kg.values())
        m_scrap_remaining = m_scrap_init  # unmelted scrap

        # Slag composition tracked by MASS (kg per oxide)
        lime_total = static_result['石灰加入量']
        dolomite_total = static_result['白云石加入量']
        m_slag_initial = lime_total * 0.9 + dolomite_total * 0.54
        slag_kg = {
            'SiO2':  m_slag_initial * 0.15,
            'CaO':   m_slag_initial * 0.45,
            'FeO':   m_slag_initial * 0.20,
            'MnO':   m_slag_initial * 0.05,
            'MgO':   m_slag_initial * 0.08,
            'P2O5':  m_slag_initial * 0.03,
        }
        lime_melted_cum = 0.0
        dolomite_melted_cum = 0.0
        sinter_melted_cum = 0.0
        # 烧结矿冷却剂用量: 根据热平衡估算 (静态模型反算)
        sinter_total = inputs['铁水质量']*18.0  # ~18kg/t铁水 = ~3t for 170t

        # Time series: full steel + slag composition
        ts = {'time_s':[], 'T_C':[], 'T_K':[],
              # Steel wt%
              'w_C':[], 'w_Si':[], 'w_Mn':[], 'w_P':[], 'w_S':[], 'w_Fe':[],
              # Slag wt%
              'w_FeO':[], 'w_SiO2':[], 'w_CaO':[], 'w_MgO':[], 'w_MnO':[], 'w_P2O5':[],
              # Masses
              'm_steel_kg':[], 'm_slag_kg':[]}

        q_CO_prev = 0.0
        O2_consumed = 0.0
        t = 0.0
        O2_total_target = static_result['总耗氧量']

        while t < max_time:
            # --- Adaptive time step ---
            w_C_pct = steel_kg.get('C',0.0)/m_steel*100.0 if m_steel>0 else 0.0
            dt = self.select_time_step(w_C_pct)

            # --- Slag mass and wt% ---
            m_slag = sum(slag_kg.values())
            m_slag = max(m_slag, 1.0)
            slag_wt = {k: v/m_slag*100.0 for k,v in slag_kg.items()}
            f_solid = calc_solid_fraction(slag_wt.get('CaO',45), slag_wt.get('SiO2',15), slag_wt.get('FeO',20))

            # --- Steel wt% ---
            steel_wt_pct = {}
            for elem in ['C','Si','Mn','P','S','V','Ti','Fe']:
                steel_wt_pct[elem] = steel_kg.get(elem,0.0)/m_steel*100.0 if m_steel>0 else 0.0

            # === Steps 1-2: Mixing energy and mass transfer ===
            eps_T = self.mixing.calc_epsilon_T(inputs['供氧流量'],m_steel,inputs['枪位高度'],
                                                inputs['喷嘴直径'],300.0,inputs['喷嘴夹角'])
            eps_B = self.mixing.calc_epsilon_B(inputs['底吹氩气流量'],T,m_steel,inputs['熔池深度'])
            eps_total = eps_T+eps_B
            # Calibrate baseline to initial operating condition so km≈km0 at start
            if t == 0.0:
                self._eps_total0 = max(eps_total, 0.01)

            km = self.transfer.calc_k_m(KITAMURA_BASE['k_m0'],eps_total,self._eps_total0,gamma=0.3)
            q_Ar = inputs['底吹氩气流量']/60000.0
            q_gas_ratio = min((q_Ar+3.0*q_CO_prev)/KITAMURA_BASE['q_Ar0'], 10.0)
            ks_prime = KITAMURA_BASE['k_s0']*q_gas_ratio*(KITAMURA_BASE['W_s0']/max(m_slag,1e-6))
            ks = self.transfer.calc_k_s(ks_prime,f_solid,KITAMURA_BASE['n_solid'])
            F_dict = self.transfer.calc_all_F(km,ks,RHO_STEEL,RHO_SLAG,MOLAR_MASS,OXIDE_MOLAR_MASS)

            # === Step 3: Equilibrium constants + E_M ===
            K_M = {elem: calc_K_M(elem,T) for elem in OXIDATION_THERMO}
            f_steel = calc_all_f_i(steel_wt_pct,T)
            slag_pct = {k:slag_wt[k] for k in ['FeO','SiO2','CaO','MnO','MgO','P2O5'] if k in slag_wt}
            f_slag = {}
            for ox in slag_pct:
                try: f_slag[ox]=calc_slag_activity_coeff(ox,slag_pct,T)
                except: f_slag[ox]=1.0

            C_slag = RHO_SLAG/0.060
            E_M = {}
            ox_map = {'Si':'SiO2','Mn':'MnO','P':'P2O5','V':'V2O3','Ti':'TiO2','Fe':'FeO'}
            for elem in ['Si','Mn','P','V','Ti','Fe']:
                ox = ox_map[elem]
                fM = f_steel.get(elem,1.0)
                fMOn = f_slag.get(ox,1.0)
                N_MOn = OXIDE_MOLAR_MASS.get(ox,0.1)
                E_M[elem] = 100.0*C_slag*N_MOn*fM*K_M[elem]/(RHO_SLAG*max(fMOn,1e-6))

            # === Step 4: Solve a*_O via Brent ===
            steel_b = {f'w([{e}])b':steel_wt_pct.get(e,0.0) for e in ['C','Si','Mn','P','V','Ti','Fe']}
            slag_b = {f'w(({o}))b':slag_wt.get(o,0.0) for o in ['SiO2','MnO','P2O5','FeO','V2O3','TiO2']}

            K_C = K_M.get('C',50.0)
            G_CO = KITAMURA_BASE['G_CO0']*(max(m_slag,1e-6)/KITAMURA_BASE['W_s0'])
            G_CO = G_CO*max(q_Ar/KITAMURA_BASE['q_Ar0']+q_CO_prev/max(inputs['供氧流量']/60.0,1e-6), 0.1)
            G_CO = max(G_CO,1e-6)

            A0_eff = KITAMURA_BASE['A0']*75.0
            A = A0_eff*(1.0-f_solid)**(2.0/3.0)

            K_Fe = K_M.get('Fe',10.0)
            a_O_b_slag = self.oxy.calc_a_O_b_from_slag(slag_wt.get('FeO',20),f_slag.get('FeO',1.0),K_Fe)
            a_O_b_slag = min(a_O_b_slag,5.0)

            O2_mol_per_s = inputs['供氧流量']/(60.0*0.0224)
            J_O2_inj = O2_mol_per_s*2.0/max(A,1e-6)*0.05

            def residual_func(a_O_star):
                w_star = self.interface.calc_all(F_dict,steel_b,slag_b,E_M,a_O_star,G_CO,K_C)
                J_fluxes = {}
                for e in ['C','Si','Mn','P','V','Ti','Fe']:
                    F_M = F_dict.get(f'F_{e}',0.0)
                    J_fluxes[e] = self.oxy.calc_J_M(F_M,steel_b.get(f'w([{e}])b',0.0),w_star.get(f'w([{e}])*',0.0))
                return self.oxy.calc_residual(J_fluxes,F_dict.get('F_O',1.0),a_O_b_slag,a_O_star,f_steel.get('O',1.0))-J_O2_inj

            try: a_O_star = self.solver.solve(residual_func,bracket=(1e-8,1.0))
            except: a_O_star = 1e-6

            # === Steps 5-8: Rates and mass balance (per kg of element) ===
            w_star_all = self.interface.calc_all(F_dict,steel_b,slag_b,E_M,a_O_star,G_CO,K_C)

            delta_kg = {}
            for elem in ['C','Si','Mn','P','V','Ti','Fe']:
                w_b = steel_wt_pct.get(elem,0.0)
                w_s = w_star_all.get(f'w([{elem}])*',0.0)
                F_M = F_dict.get(f'F_{elem}',0.0)
                M_i = MOLAR_MASS.get(elem,0.05)
                # Clamp oxidation direction: for Si/Mn/P ensure net oxidation when element present
                if elem in ('Si','Mn','P') and w_b > 0.01:
                    w_s = min(w_s, w_b*0.5)  # ensure interface is lower -> oxidation
                delta_kg[elem] = A*F_M*M_i*(w_b-w_s)*dt

            # Steel update
            for elem in ['C','Si','Mn','P','V','Ti']:
                steel_kg[elem] = max(0.0, steel_kg.get(elem,0.0)-delta_kg.get(elem,0.0))
            steel_kg['Fe'] = max(0.0, steel_kg.get('Fe',0.0)-delta_kg.get('Fe',0.0))
            m_steel = sum(steel_kg.values())

            # Slag update by oxide MASS (kg)
            elem_to_oxide = {
                'Si': ('SiO2',  0.06009/0.02809),  # kg_SiO2 per kg_Si
                'Mn': ('MnO',   0.07094/0.05494),
                'P':  ('P2O5',  0.14194/0.03097/2.0),
                'Fe': ('FeO',   0.07185/0.05585),
            }
            for elem, (ox, factor) in elem_to_oxide.items():
                dm_elem = delta_kg.get(elem,0.0)
                if dm_elem > 0:  # oxidation: element -> slag
                    slag_kg[ox] = slag_kg.get(ox,0.0) + dm_elem*factor
                # reduction (dm<0): skip — oxide doesn't go back to steel in practice

            # Normalize slag (cap FeO at 35%)
            m_slag = sum(slag_kg.values())
            if m_slag > 0:
                w_FeO_test = slag_kg.get('FeO',0.0)/m_slag*100.0
                if w_FeO_test > 35.0:
                    slag_kg['FeO'] = 0.35*m_slag
                    # redistribute excess to CaO and SiO2
                    excess = (w_FeO_test-35.0)/100.0*m_slag
                    slag_kg['CaO'] = slag_kg.get('CaO',0.0)+excess*0.6
                    slag_kg['SiO2'] = slag_kg.get('SiO2',0.0)+excess*0.4
            m_slag = sum(slag_kg.values())

            # q_CO for next step
            q_CO_prev = delta_kg.get('C',0.0)/MOLAR_MASS['C']/max(dt,1e-6)*0.0224 if dt>0 else 0.0

            # === Steps 9-10: Lime and scrap incremental melting ===
            lime_melt_rate = lime_total/900.0   # kg/s, linear over 15 min
            dl = lime_melt_rate*dt
            lime_melted_cum += dl
            slag_kg['CaO'] += dl*0.90
            slag_kg['MgO'] += dl*0.03
            slag_kg['SiO2'] += dl*0.02

            dolo_melt_rate = dolomite_total/900.0
            dd = dolo_melt_rate*dt
            dolomite_melted_cum += dd
            slag_kg['CaO'] += dd*0.30
            slag_kg['MgO'] += dd*0.21
            slag_kg['SiO2'] += dd*0.02

            # Sinter coolant: distributed over first 10 min (600s), Fe2O3→FeO
            sinter_rate = sinter_total/600.0 if t < 600.0 else 0.0
            dsinter = sinter_rate*dt
            sinter_melted_cum += dsinter
            slag_kg['FeO'] += dsinter*0.09   # FeO from sinter
            slag_kg['FeO'] += dsinter*0.64*0.9  # Fe2O3→FeO (partial reduction)
            slag_kg['CaO'] += dsinter*0.13
            slag_kg['SiO2'] += dsinter*0.07
            slag_kg['MgO'] += dsinter*0.02
            slag_kg['MnO'] += dsinter*0.015

            # Scrap: melt proportionally to time and temperature
            T_MP_scrap = self.scrap.calc_T_MP(0.2)
            if T > T_MP_scrap:
                scrap_rate = m_scrap_init/600.0  # melt over 10 min when T > T_MP
            else:
                scrap_rate = m_scrap_init/1200.0*0.3
            ds = min(scrap_rate*dt, m_scrap_remaining)
            m_scrap_remaining -= ds
            steel_kg['Fe'] += ds
            m_steel = sum(steel_kg.values())

            # === Steps 11-12: Heat balance (physics-based) ===

            # 热收入：元素氧化放热 + 二次燃烧 + C2S/C3P形成热
            ox_heat = {}
            for elem in ['C','Si','Mn','P','Fe']:
                dm = delta_kg.get(elem,0.0)
                if abs(dm) > 0:
                    ox_heat[elem] = abs(dm)/MOLAR_MASS[elem]*OXIDATION_THERMO[elem]['delta_H']
            Q_oxid = sum(ox_heat.values())
            Q_post = ox_heat.get('C',0.0)*0.12  # ~12% CO -> CO2 post-combustion
            # SiO2 + CaO formation heat
            dm_Si = delta_kg.get('Si',0.0)
            if dm_Si > 0:
                Q_C2S = dm_Si/MOLAR_MASS['Si']*120000.0  # J
            else:
                Q_C2S = 0.0
            Q_HG = Q_oxid + Q_post + Q_C2S

            # 热支出（物理吸热项）：
            # 1. 废钢升温+熔化: 25°C->1450°C 显热~750 + 潜热~260 = ~1010 kJ/kg
            Q_scrap = ds*1010000.0

            # 2. 石灰/白云石分解吸热
            Q_CM = dl*1780000.0 + dd*3020000.0

            # 3. 新生成炉渣升温和熔化 (25°C -> T)
            slag_new_kg = (dl*0.95 + dd*0.53 +
                           sum(max(0.0, delta_kg.get(e,0.0)) for e in ('Si','Mn','P','Fe')))
            if T > 298 and slag_new_kg > 0:
                Q_slag_heat = slag_new_kg*CP_SLAG*(T-298.0)*0.8  # 0.8: slag absorbs most sensible heat
            else:
                Q_slag_heat = 0.0

            # 4. CO炉气带走显热 (Cp_CO≈33 J/mol·K, full sensible heat)
            if q_CO_prev > 0 and T > 298:
                Q_offgas = q_CO_prev/0.0224*33.0*(T-298.0)*dt  # J
            else:
                Q_offgas = 0.0

            # 5. 辐射+炉衬+粉尘固定热损失
            Q_rad_lining = 200000.0*dt

            # 6. 烧结矿冷却吸热: 升温显热+Fe2O3还原≈1.8 MJ/kg
            Q_sinter = dsinter*1800000.0

            Q_consumed = Q_scrap + Q_CM + Q_slag_heat + Q_offgas + Q_rad_lining + Q_sinter

            delta_T = self.heat_bal.calc_delta_T(max(0,Q_HG), max(0,Q_consumed), 0.0,
                                                  CP_STEEL, m_steel, CP_SLAG, m_slag)
            T += max(-50.0, min(delta_T, 100.0))

            # === Record full state ===
            m_slag = sum(slag_kg.values())
            ts['time_s'].append(t)
            ts['T_C'].append(T-273.15)
            ts['T_K'].append(T)
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
            ts['m_steel_kg'].append(m_steel)
            ts['m_slag_kg'].append(m_slag)

            t += dt
            O2_consumed += inputs['供氧流量']*dt/60.0
            w_C_pct = steel_kg.get('C',0)/m_steel*100 if m_steel>0 else 0

            if O2_consumed >= O2_total_target*0.5 and w_C_pct <= inputs['终点碳目标']:
                break

        return {
            '终点C': round(ts['w_C'][-1],4) if ts['w_C'] else 0.0,
            '终点T': round(ts['T_C'][-1],1) if ts['T_C'] else 0.0,
            '时序数据': ts,
            '步数': len(ts['time_s']),
        }
