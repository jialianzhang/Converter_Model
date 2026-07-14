"""参数标定优化器: 70/30 训练验证划分 + 有界优化 + 命中率报告

泛化性验证内置于流程: 参数只在训练集上优化, 命中率在从未参与优化的验证集上
复核 —— 验证集命中率 ≈ 训练集 ⇒ 参数抓住的是设备特性而非炉次噪声。
支持 multiprocessing 并行 (n_jobs>1): 100炉 × 数百次评估的必要加速。
Windows 下调用方脚本必须置于 if __name__ == '__main__' 保护内。
"""
import time
import numpy as np
from scipy.optimize import minimize
from config import CALIBRATION_DEFAULTS, CALIBRATION_BOUNDS
from .objective import (PARAM_NAMES, calibration_objective, evaluate_heats,
                        vector_to_params, params_to_vector, hit_rates)


class Calibrator:
    def __init__(self, seed: int = 42):
        self.seed = seed

    def split(self, heats: list, train_frac: float = 0.7):
        """随机 70/30 划分 (固定种子保证可复现)"""
        rng = np.random.default_rng(self.seed)
        idx = rng.permutation(len(heats))
        k = max(1, int(round(len(heats) * train_frac)))
        train = [heats[i] for i in idx[:k]]
        val = [heats[i] for i in idx[k:]]
        return train, val

    def calibrate(self, heats: list, method: str = 'Nelder-Mead',
                  maxiter: int = 250, w_C: float = 1.0, w_T: float = 1.0,
                  w_FeO: float = 1.0, w_aO: float = 0.0, w_W: float = 0.5,
                  w_Mn: float = 0.0, w_P: float = 0.0,
                  train_frac: float = 0.7, n_jobs: int = 1, verbose: bool = True,
                  framing: str = 'oxygen', fixed: dict = None):
        # w_aO 默认 0: 模型 a_O* 为界面氧活度, 与 TSO 实测本体溶解氧有系统偏置(~0.21),
        #   原始差值不可直接作目标(会以 ~114 量级压倒 C/T)。防 T↔FeO 补偿的 FeO 约束
        #   改由出钢量(铁平衡)承担; TSO 氧作独立验证指标 (hit_rates 中报告 MAE_aO)。
        # framing='target_c': 目标碳终止(使用框架), C 构造性命中, 请置 w_C=0 并启用 w_Mn/w_P。
        # fixed: 固定不参与优化的参数 {名: 值}, 如 {'eta_O2_late': 0.45} (按物理定死)。
        """返回: {'参数': dict, '训练集': 命中率, '验证集': 命中率, 'loss': float}"""
        train, val = self.split(heats, train_frac)
        defaults = dict(CALIBRATION_DEFAULTS)
        if fixed:
            defaults.update(fixed)
        x0 = params_to_vector(defaults)
        bounds = [(defaults[k], defaults[k] + 1e-12) if fixed and k in fixed
                  else CALIBRATION_BOUNDS[k] for k in PARAM_NAMES]

        pool = None
        try:
            if n_jobs and n_jobs > 1:
                from multiprocessing import Pool
                pool = Pool(processes=n_jobs)
                if verbose:
                    print(f"  并行进程数: {n_jobs}", flush=True)

            counter = {'n': 0, 't0': time.time(), 'best': float('inf')}

            def obj(x):
                v = calibration_objective(x, train, w_C, w_T, w_FeO, w_aO, w_W,
                                          w_Mn, w_P, pool, framing)
                counter['n'] += 1
                counter['best'] = min(counter['best'], v)
                if verbose and counter['n'] % 20 == 0:
                    print(f"  eval {counter['n']:4d}: loss={v:.3f} (best={counter['best']:.3f}) "
                          f"用时{time.time() - counter['t0']:.0f}s", flush=True)
                return v

            result = minimize(
                obj, x0, method=method, bounds=bounds,
                options={'maxiter': maxiter, 'disp': verbose,
                         'xatol': 1e-3, 'fatol': 1e-3} if method == 'Nelder-Mead'
                        else {'maxiter': maxiter, 'disp': verbose},
            )

            best = vector_to_params(result.x)
            report = {
                '参数': {k: round(v, 5) for k, v in best.items()},
                'loss': round(float(result.fun), 4),
                '评估次数': counter['n'],
                '迭代': int(result.nit) if hasattr(result, 'nit') else None,
                '训练集': hit_rates(evaluate_heats(best, train, pool, framing)),
                '训练集炉数': len(train),
                '框架': framing,
                '固定参数': fixed or {},
            }
            if val:
                report['验证集'] = hit_rates(evaluate_heats(best, val, pool, framing))
                report['验证集炉数'] = len(val)
            # 边界触碰提示 (参数打到物理边界 → 可能存在未建模效应, 需人工审视)
            at_bounds = [k for k, v in best.items()
                         if abs(v - CALIBRATION_BOUNDS[k][0]) < 1e-6 * max(abs(CALIBRATION_BOUNDS[k][0]), 1)
                         or abs(v - CALIBRATION_BOUNDS[k][1]) < 1e-6 * max(abs(CALIBRATION_BOUNDS[k][1]), 1)]
            if at_bounds:
                report['触及边界的参数'] = at_bounds
            return report
        finally:
            if pool is not None:
                pool.close()
                pool.join()
