"""参数标定目标函数"""
import numpy as np

def calibration_objective(params, engine, heats_data):
    """目标函数: weighted MSE of C and T predictions"""
    alpha_steel, alpha_slag, G_CO0, gamma, eta_post, alpha_rad, epsilon_dust = params

    # Update engine parameters
    from config import KITAMURA_BASE
    KITAMURA_BASE['k_m0'] = alpha_steel * 2.0e-4
    KITAMURA_BASE['k_s0'] = alpha_slag * 1.0e-5
    KITAMURA_BASE['G_CO0'] = G_CO0

    errors_C = []
    errors_T = []
    for heat in heats_data:
        result = engine.run(heat)
        C_actual = heat.get('终点C_actual', heat.get('终点碳目标', 0.04))
        T_actual = heat.get('终点T_actual', heat.get('终点温度目标', 1650))
        errors_C.append((result['终点C'] - C_actual) ** 2)
        errors_T.append((result['终点T'] - T_actual) ** 2)

    w_C, w_T = 1.0, 0.01
    return np.mean(w_C * np.array(errors_C) + w_T * np.array(errors_T))
