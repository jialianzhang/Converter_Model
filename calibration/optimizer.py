"""参数标定优化器"""
import numpy as np
from scipy.optimize import minimize
from .objective import calibration_objective

class Calibrator:
    def __init__(self, engine):
        self.engine = engine

    def calibrate(self, heats_data, initial_params=None):
        if initial_params is None:
            initial_params = [1.0, 1.0, 0.01, 0.3, 0.10, 0.03, 0.01]

        bounds = [
            (0.5, 5.0),     # alpha_steel
            (0.5, 5.0),     # alpha_slag
            (0.001, 0.1),   # G_CO0
            (0.1, 0.8),     # gamma
            (0.05, 0.30),   # eta_post
            (0.01, 0.10),   # alpha_rad
            (0.005, 0.05),  # epsilon_dust
        ]

        result = minimize(
            calibration_objective,
            initial_params,
            args=(self.engine, heats_data),
            method='L-BFGS-B',
            bounds=bounds,
            options={'maxiter': 200, 'disp': True}
        )
        return result
