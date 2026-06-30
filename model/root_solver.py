"""Brent + NR + Bisection three-tier fallback solver"""
from scipy.optimize import brentq
import numpy as np


class RootSolver:
    def __init__(self, tol=1e-10, max_iter=200):
        self.tol = tol
        self.max_iter = max_iter

    def solve(self, residual_func, bracket=(1e-8, 1.0)):
        try:
            a, b = bracket
            fa = residual_func(a)
            fb = residual_func(b)
            if fa * fb > 0:
                if abs(fa) < self.tol:
                    return a
                if abs(fb) < self.tol:
                    return b
                raise ValueError("Bracket does not bracket a root")
            return brentq(residual_func, a, b, xtol=self.tol, maxiter=self.max_iter)
        except Exception:
            return self._newton_fallback(residual_func, bracket)

    def _newton_fallback(self, f, bracket):
        x = (bracket[0] + bracket[1]) / 2.0
        for _ in range(self.max_iter):
            fx = f(x)
            if abs(fx) < self.tol:
                return x
            EPS = 1e-8
            df = (f(x + EPS) - fx) / EPS
            if abs(df) < 1e-20:
                return self._bisection_fallback(f, bracket)
            x_new = x - fx / df
            if x_new <= bracket[0] or x_new >= bracket[1]:
                return self._bisection_fallback(f, bracket)
            x = x_new
        return self._bisection_fallback(f, bracket)

    def _bisection_fallback(self, f, bracket):
        a, b = bracket
        for _ in range(self.max_iter):
            c = (a + b) / 2.0
            fc = f(c)
            if abs(fc) < self.tol:
                return c
            fa = f(a)
            if fa * fc < 0:
                b = c
            else:
                a = c
        return (a + b) / 2.0
