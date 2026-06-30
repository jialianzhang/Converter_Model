import pytest
import numpy as np
from model.root_solver import RootSolver


def f_linear(x):
    return x - 0.5


def test_brent_solves_linear():
    rs = RootSolver()
    root = rs.solve(f_linear, bracket=(0.0, 1.0))
    assert abs(root - 0.5) < 1e-8


def f_quadratic(x):
    return x**2 - 0.25


def test_brent_solves_quadratic():
    rs = RootSolver()
    root = rs.solve(f_quadratic, bracket=(0.0, 1.0))
    assert abs(root - 0.5) < 1e-8


def test_newton_fallback():
    rs = RootSolver(tol=1e-10)
    root = rs.solve(lambda x: x - 0.3, bracket=(0.0, 1.0))
    assert abs(root - 0.3) < 1e-8


def test_bisection_fallback():
    rs = RootSolver(tol=1e-10)
    root = rs.solve(lambda x: x - 0.7, bracket=(0.0, 1.0))
    assert abs(root - 0.7) < 1e-8
