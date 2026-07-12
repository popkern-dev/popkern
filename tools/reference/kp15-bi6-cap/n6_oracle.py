"""Faithful sympy port of the KP15 n=6 Maple script (docs/papers/n6.txt).

Reference oracle for the Rust bi-6 polar-cap construction
(crates/popkern-core/src/geometry/gregory/kp15_bi6_*.rs). It reproduces
the paper "Minimal bi-6 G^2 completion of bicubic spline surfaces"
(Karciauskas & Peters, 2015) end to end and re-runs the script's own
smoothness verification: for every sector the residuals

    G1[s] = Fv + Gv - 2c(1-u)^2 Fu
    G2[s] = -Gvv + Fvv - 4c(1-u)^2 Fuv + 4c^2(1-u)^4 Fuu
            + c(1-u)^3 (A u - 8c - A) Fu + A(1-u)^2 Fv,   A = 24c(tau-1)/(1-u+tau u)

(all evaluated at v=0) must vanish, together with the deferred `sal`
residual lyc[0].

CRITICAL INDEXING FLIP (script header): script `p[i,j]` = paper
`p[6-i,6-j]`; script `[6,6]` is the extraordinary point. All CP indices
in this file are in SCRIPT convention.

Run: see README.md.  This module exposes `build_cap()` and residual
helpers; CLI modes live in __main__ at the bottom.
"""

from __future__ import annotations

import argparse
import random
from fractions import Fraction

import sympy as sp

U, V = sp.symbols("u v")


# --------------------------------------------------------------------------
# BtoBB : B-spline 4x4 net -> bi-3 BB-form (script lines 78-113)
# --------------------------------------------------------------------------
def BtoBB(P):
    """P is a 4x4 list, P[row][col].  Returns b(i,j) accessor as a dict
    keyed (i,j), i,j=0..3.  Script: p[i,j] := P[j+1][i+1]."""

    def p(i, j):
        return P[j][i]

    b = {}
    b[(0, 0)] = (
        16 * p(1, 1) + 4 * (p(1, 0) + p(2, 1) + p(1, 2) + p(0, 1)) + p(0, 0) + p(2, 0) + p(2, 2) + p(0, 2)
    ) / 36
    b[(3, 0)] = (
        16 * p(2, 1) + 4 * (p(2, 0) + p(3, 1) + p(2, 2) + p(1, 1)) + p(1, 0) + p(3, 0) + p(3, 2) + p(1, 2)
    ) / 36
    b[(0, 3)] = (
        16 * p(1, 2) + 4 * (p(1, 1) + p(2, 2) + p(1, 3) + p(0, 2)) + p(0, 1) + p(2, 1) + p(2, 3) + p(0, 3)
    ) / 36
    b[(3, 3)] = (
        16 * p(2, 2) + 4 * (p(2, 1) + p(3, 2) + p(2, 3) + p(1, 2)) + p(1, 1) + p(3, 1) + p(3, 3) + p(1, 3)
    ) / 36
    b[(1, 0)] = sp.Rational(4, 9) * p(1, 1) + sp.Rational(2, 9) * p(2, 1) + sp.Rational(1, 9) * (p(1, 0) + p(1, 2)) + sp.Rational(1, 18) * (p(2, 0) + p(2, 2))
    b[(2, 0)] = sp.Rational(2, 9) * p(1, 1) + sp.Rational(4, 9) * p(2, 1) + sp.Rational(1, 18) * (p(1, 0) + p(1, 2)) + sp.Rational(1, 9) * (p(2, 0) + p(2, 2))
    b[(1, 3)] = sp.Rational(4, 9) * p(1, 2) + sp.Rational(2, 9) * p(2, 2) + sp.Rational(1, 9) * (p(1, 1) + p(1, 3)) + sp.Rational(1, 18) * (p(2, 1) + p(2, 3))
    b[(2, 3)] = sp.Rational(2, 9) * p(1, 2) + sp.Rational(4, 9) * p(2, 2) + sp.Rational(1, 18) * (p(1, 1) + p(1, 3)) + sp.Rational(1, 9) * (p(2, 1) + p(2, 3))
    b[(0, 1)] = sp.Rational(4, 9) * p(1, 1) + sp.Rational(2, 9) * p(1, 2) + sp.Rational(1, 9) * (p(0, 1) + p(2, 1)) + sp.Rational(1, 18) * (p(0, 2) + p(2, 2))
    b[(0, 2)] = sp.Rational(2, 9) * p(1, 1) + sp.Rational(4, 9) * p(1, 2) + sp.Rational(1, 18) * (p(0, 1) + p(2, 1)) + sp.Rational(1, 9) * (p(0, 2) + p(2, 2))
    b[(3, 1)] = sp.Rational(4, 9) * p(2, 1) + sp.Rational(2, 9) * p(2, 2) + sp.Rational(1, 9) * (p(1, 1) + p(3, 1)) + sp.Rational(1, 18) * (p(1, 2) + p(3, 2))
    b[(3, 2)] = sp.Rational(2, 9) * p(2, 1) + sp.Rational(4, 9) * p(2, 2) + sp.Rational(1, 18) * (p(1, 1) + p(3, 1)) + sp.Rational(1, 9) * (p(1, 2) + p(3, 2))
    b[(1, 1)] = sp.Rational(4, 9) * p(1, 1) + sp.Rational(2, 9) * p(2, 1) + sp.Rational(2, 9) * p(1, 2) + sp.Rational(1, 9) * p(2, 2)
    b[(2, 1)] = sp.Rational(4, 9) * p(2, 1) + sp.Rational(2, 9) * p(1, 1) + sp.Rational(2, 9) * p(2, 2) + sp.Rational(1, 9) * p(1, 2)
    b[(1, 2)] = sp.Rational(4, 9) * p(1, 2) + sp.Rational(2, 9) * p(1, 1) + sp.Rational(2, 9) * p(2, 2) + sp.Rational(1, 9) * p(2, 1)
    b[(2, 2)] = sp.Rational(4, 9) * p(2, 2) + sp.Rational(2, 9) * p(2, 1) + sp.Rational(2, 9) * p(1, 2) + sp.Rational(1, 9) * p(1, 1)
    return b


# --------------------------------------------------------------------------
# TBdata : reparameterised tensor border (script lines 41-74)
#   pin : dict (i,j)->value for i=0..3, j=0..2   ; returns rows[j][i], j=0..2, i=0..6
# --------------------------------------------------------------------------
def TBdata(pin, c, d2, tau):
    def p(i, j):
        return pin[(i, j)]

    R = sp.Rational
    row0 = [
        p(0, 0),
        R(1, 2) * p(0, 0) + R(1, 2) * p(1, 0),
        R(1, 5) * p(0, 0) + R(1, 5) * p(2, 0) + R(3, 5) * p(1, 0),
        R(1, 20) * p(0, 0) + R(9, 20) * p(1, 0) + R(9, 20) * p(2, 0) + R(1, 20) * p(3, 0),
        R(1, 5) * p(3, 0) + R(1, 5) * p(1, 0) + R(3, 5) * p(2, 0),
        R(1, 2) * p(2, 0) + R(1, 2) * p(3, 0),
        p(3, 0),
    ]
    row1 = [
        R(1, 2) * p(0, 0) + R(1, 2) * p(0, 1),
        R(1, 4) * p(0, 0) + R(1, 4) * p(1, 0) + R(1, 4) * p(0, 1) + R(1, 4) * p(1, 1),
        R(1, 10) * p(0, 0) + R(1, 10) * p(0, 1) + R(3, 10) * p(1, 0) + R(1, 10) * p(2, 0) + R(3, 10) * p(1, 1) + R(1, 10) * p(2, 1),
        R(1, 40) * p(0, 0) + R(9, 40) * p(1, 0) + R(9, 40) * p(2, 0) + R(1, 40) * p(3, 0) + R(1, 40) * p(0, 1) + R(9, 40) * p(1, 1) + R(9, 40) * p(2, 1) + R(1, 40) * p(3, 1),
        R(1, 10) * p(1, 0) + R(3, 10) * p(2, 0) + R(1, 10) * p(3, 0) + R(1, 10) * p(1, 1) + R(3, 10) * p(2, 1) + R(1, 10) * p(3, 1),
        R(1, 4) * p(2, 0) + R(1, 4) * p(3, 0) + R(1, 4) * p(2, 1) + R(1, 4) * p(3, 1),
        R(1, 2) * p(3, 0) + R(1, 2) * p(3, 1),
    ]
    row2 = [
        R(3, 5) * p(0, 1) + R(1, 5) * p(0, 0) + R(1, 5) * p(0, 2),
        R(1, 10) * p(1, 2) + R(1, 10) * p(0, 0) + R(1, 10) * p(1, 0) + R(3, 10) * p(0, 1) + R(3, 10) * p(1, 1) + R(1, 10) * p(0, 2),
        R(3, 25) * p(1, 2) + R(1, 25) * p(2, 2) + R(1, 25) * p(0, 0) + R(3, 25) * p(1, 0) + R(1, 25) * p(2, 0)
        + R(3, 25) * p(0, 1) + R(9, 25) * p(1, 1) + R(3, 25) * p(2, 1) + R(1, 25) * p(0, 2)
        + R(1, 50) * p(0, 1) * d2 - R(1, 25) * d2 * p(0, 0) + R(1, 50) * d2 * p(1, 0),
        R(1, 600) * (
            18 * d2 * p(2, 0) * tau + 3 * p(0, 1) * tau * d2 + 2 * p(0, 1) * tau * c - 12 * p(1, 0) * c * tau
            + 10 * p(0, 0) * c * tau - 6 * p(0, 0) * d2 * tau + 27 * p(1, 1) * d2 * tau + 12 * p(1, 0) * c
            + 54 * p(1, 2) * tau + 162 * p(2, 1) * tau + 18 * p(0, 1) * tau + 162 * p(1, 1) * tau
            + 54 * p(2, 0) * tau + 6 * p(3, 0) * tau + 6 * p(0, 0) * tau + 54 * p(1, 0) * tau
            + 54 * p(2, 2) * tau + 6 * p(3, 2) * tau + 18 * p(3, 1) * tau - 42 * d2 * p(1, 0) * tau
            + 6 * p(0, 2) * tau - 12 * p(0, 0) * c
        ) / tau,
        R(1, 150) * (
            6 * p(1, 2) * tau + 18 * p(2, 2) * tau + 6 * p(3, 2) * tau + 6 * p(1, 0) * tau + 18 * p(2, 0) * tau
            + 6 * p(3, 0) * tau + 18 * p(1, 1) * tau + 54 * p(2, 1) * tau + 18 * p(3, 1) * tau + 9 * p(2, 1) * d2 * tau
            + 2 * p(1, 1) * tau * c + 3 * p(1, 1) * d2 * tau - 5 * d2 * p(1, 0) * tau + 6 * p(1, 0) * c * tau
            - 8 * p(1, 0) * c - 10 * d2 * p(2, 0) * tau - 8 * p(2, 0) * c * tau + 8 * p(2, 0) * c + 3 * d2 * p(3, 0) * tau
        ) / tau,
        R(1, 60) * (
            6 * p(2, 2) * tau + 6 * p(3, 2) * tau + 6 * p(2, 0) * tau + 6 * p(3, 0) * tau + 18 * p(2, 1) * tau
            + 18 * p(3, 1) * tau + 3 * p(3, 1) * d2 * tau + 2 * p(2, 1) * tau * c + 3 * p(2, 1) * d2 * tau
            - 4 * d2 * p(2, 0) * tau + 2 * p(2, 0) * c * tau - 4 * p(2, 0) * c - 2 * d2 * p(3, 0) * tau
            - 4 * p(3, 0) * c * tau + 4 * p(3, 0) * c
        ) / tau,
        R(1, 5) * p(3, 2) + R(1, 5) * p(3, 0) + R(3, 5) * p(3, 1) + R(1, 15) * p(3, 1) * c + R(1, 10) * p(3, 1) * d2
        - R(1, 15) * p(3, 0) * c - R(1, 10) * d2 * p(3, 0),
    ]
    return [row0, row1, row2]


# --------------------------------------------------------------------------
# mu coefficient tables (script lines 157-176, 191-202, 270-307)
# --------------------------------------------------------------------------
def mu1_tab(c, tau):
    m = {1: sp.Rational(2, 5) * c, 2: -c + 2, 3: c, 4: -sp.Rational(2, 5) * c, 5: sp.Rational(1, 15) * c}
    m[0] = 2 - (m[1] + m[2] + m[3] + m[4] + m[5])
    return m


def mu2_tab(c, tau):
    m = {
        1: -sp.Rational(4, 5) * c * tau**2,
        2: sp.Rational(24, 5) * c**2 * tau + sp.Rational(4, 5) * c**2 - 4 * c + sp.Rational(4, 5) * c**2 * tau**2 + 4 - sp.Rational(24, 5) * c * tau,
        3: sp.Rational(24, 5) * c * tau - 4,
        4: -sp.Rational(2, 3) * c * (-11 + 6 * c * tau + 3 * c),
        5: -4 * c,
        6: 2 * c**2,
        7: -sp.Rational(4, 5) * c**2,
        8: sp.Rational(2, 15) * c**2,
    }
    m[0] = -sum(m[k] for k in range(1, 9))
    return m


def mu3_tab(c, tau):
    m = {
        1: sp.Rational(4, 25) * tau**3 * c,
        2: -sp.Rational(4, 25) * c * (6 * c * tau**2 - 6 * tau**2 - 5 + 3 * c * tau + tau**3 * c),
        3: -sp.Rational(24, 25) * c * tau**2,
        4: -2 * c + sp.Rational(4, 5) * c**2 * tau**2 + sp.Rational(6, 5) * c**2 * tau + 4 - sp.Rational(12, 5) * c * tau,
        5: sp.Rational(12, 5) * c * tau - 4,
        6: -sp.Rational(2, 5) * c * (3 * c * tau - 5),
        7: sp.Rational(4, 25) * c * (3 * c * tau - 5),
        8: -sp.Rational(2, 75) * c * (3 * c * tau - 5),
        9: -sp.Rational(3, 5) * c / tau,
        10: -sp.Rational(6, 125) * c / tau**2 * (-5 + c * tau),
        11: sp.Rational(1, 125) * c * (-5 + 6 * c * tau) / tau**3,
    }
    m[0] = -sum(m[k] for k in range(1, 9))
    return m


def mut_tab(c, tau):
    m = {
        1: -sp.Rational(11, 20),
        2: sp.Rational(1, 8) * (-10 + 11 * c) / c,
        3: sp.Rational(5, 4) / c,
        4: -sp.Rational(1, 4) * c + sp.Rational(1, 8),
        5: sp.Rational(3, 10) * c - sp.Rational(1, 20),
        6: -sp.Rational(1, 20) * c + sp.Rational(1, 120),
        7: sp.Rational(3, 8) / tau,
        8: sp.Rational(3, 100) / tau**2 * (c * tau - 5),
        9: -sp.Rational(1, 200) * (-5 + 6 * c * tau) / tau**3,
        10: -sp.Rational(5, 16) / c,
    }
    m[0] = 1 - sum(m[k] for k in range(1, 7))
    return m


def mub_tab(c, tau):
    m = {
        1: 1,
        2: -sp.Rational(5, 2) * (-2 + c) / c,
        3: -5 / c,
        4: sp.Rational(5, 2),
        5: -1,
        6: sp.Rational(1, 6),
        7: sp.Rational(3, 10) / tau * (c * tau - 5),
        8: sp.Rational(3, 125) / tau**2 * (-20 * c * tau + 25 + c**2 * tau**2),
        9: -sp.Rational(1, 250) * (-45 * c * tau + 25 + 6 * c**2 * tau**2) / tau**3,
        10: sp.Rational(5, 4) / c,
    }
    m[0] = -sum(m[k] for k in range(1, 7))
    return m


# --------------------------------------------------------------------------
# Full cap construction.
# net : list of n 4x4 arrays (sympy values), one per sector.
# free : dict with keys P00, TAN0, TAN1, TWO0, TWO1, TWO2, P220, P210,
#        P23[r] (r=0..n-1), P33[r] (r=0..n-1).
# Returns (P, sal, meta) where P[s] is dict (i,j)->expr.
# --------------------------------------------------------------------------
def build_cap(n, c, tau, net, free, simplify=False):
    sf = lambda k: k % n
    d2 = sp.together((sp.Rational(1, 2) * (3 * tau**2 + tau**2 * c - 3 * tau - 2 * tau * c + 1) / tau**2))

    P = [dict() for _ in range(n)]

    def simp(x):
        return sp.factor(x) if simplify else sp.expand(x)

    # ---- Border data first (independent of interior). script 225-255 ----
    for s in range(n):
        b = BtoBB(net[s])  # b[(i,j)] i,j=0..3
        # tad: TBdata over pp[i,j]=b[i,j]  (i=0..3, j=0..2)
        pin_d = {(i, j): b[(i, j)] for i in range(4) for j in range(3)}
        # tak: transposed, pp[j,i]=b[j,i] -> TBdata over p(i,j)=b[j,i]
        pin_k = {(i, j): b[(j, i)] for i in range(4) for j in range(3)}
        tad = TBdata(pin_d, c, d2, tau)  # rows[j][i], j=0..2, i=0..6
        tak = TBdata(pin_k, c, d2, tau)
        for i in range(7):
            for j in range(3):
                P[s][(i, j)] = tad[j][i]  # p[s][i,j]
                P[s][(j, i)] = tak[j][i]  # p[s][j,i]

    # ---- (S1) central quadratic  script 128-153 ----
    for s in range(n):
        P[s][(6, 6)] = free["P00"]
    P[0][(5, 6)] = free["TAN0"]
    P[0][(6, 5)] = free["TAN1"]
    P[0][(4, 6)] = free["TWO0"]
    P[0][(5, 5)] = free["TWO1"]
    P[0][(6, 4)] = free["TWO2"]
    # pre-assign the free interior CPs (renamed at script end).
    # p[0][4,5] = P210 is handled inside the eq (7) block because it
    # is free only at n in {3, 6}.
    P[0][(4, 4)] = free["P220"]
    for r in range(n):
        P[r][(4, 3)] = free["P23"][r]
        P[r][(3, 3)] = free["P33"][r]

    for s in range(n - 1):
        t = sf(s + 1)
        P[t][(5, 6)] = P[s][(6, 5)]
        P[t][(4, 6)] = P[s][(6, 4)]
        P[t][(6, 5)] = -P[s][(5, 6)] + 2 * c * P[s][(6, 5)] + 2 * (1 - c) * P[s][(6, 6)]
        P[t][(5, 5)] = (
            sp.Rational(5, 3) * P[s][(6, 4)] * c + sp.Rational(1, 3) * c * P[s][(6, 6)]
            - 2 * P[s][(6, 5)] * c - P[s][(5, 5)] + 2 * P[s][(6, 5)]
        )
        P[t][(6, 4)] = (
            sp.Rational(24, 5) * tau * c * P[s][(5, 6)] - sp.Rational(24, 5) * tau * c**2 * P[s][(6, 5)]
            - sp.Rational(24, 5) * tau * c * P[s][(6, 6)] + sp.Rational(24, 5) * tau * c**2 * P[s][(6, 6)]
            - 4 * P[s][(5, 6)] + sp.Rational(44, 5) * P[s][(6, 5)] * c + 4 * P[s][(6, 6)] - 4 * c * P[s][(6, 6)]
            + P[s][(4, 6)] + 4 * P[s][(6, 4)] * c**2 + sp.Rational(4, 5) * c**2 * P[s][(6, 6)]
            - sp.Rational(24, 5) * P[s][(6, 5)] * c**2 - sp.Rational(24, 5) * P[s][(5, 5)] * c
        )

    # ---- (S2) eq (7): circulant for p[r][5,4] (r=0..n-1), p[r][4,5] ----
    #
    # Valence-aware free-set: the full 2n x 2n eq (7) circulant is
    # rank-deficient (rank 2n-1) ONLY at n in {3, 6}; there p[0][4,5]
    # = P210 is a genuine free parameter and only 2n-1 unknowns are
    # solved.  At every other n the circulant is full rank, so
    # p[0][4,5] is DETERMINED (P210 is unused).  This matches the
    # Rust `needs_pin = n == 3 || n == 6` in kp15_bi6_circulant.rs.
    mu1 = mu1_tab(c, tau)
    mu2 = mu2_tab(c, tau)
    p0_45_free = n in (3, 6)
    x54 = [sp.Symbol(f"x54_{r}") for r in range(n)]
    x45 = [sp.Symbol(f"x45_{r}") for r in range(n)]
    for r in range(n):
        P[r][(5, 4)] = x54[r]
    r_start_45 = 1 if p0_45_free else 0
    if p0_45_free:
        P[0][(4, 5)] = free["P210"]  # keep pinned free value
    for r in range(r_start_45, n):
        P[r][(4, 5)] = x45[r]
    lyg = [None] * (2 * n)
    for r in range(n):
        t = sf(r + 1)
        lyg[2 * r + 1] = P[t][(4, 5)] - (
            -P[r][(5, 4)] + mu1[0] * P[r][(6, 6)] + mu1[1] * P[r][(6, 5)] + mu1[2] * P[r][(6, 4)]
            + mu1[3] * P[r][(6, 2)] + mu1[4] * P[r][(6, 1)] + mu1[5] * P[r][(6, 0)]
        )
        lyg[2 * r] = P[t][(5, 4)] - (
            P[r][(4, 5)] + mu2[0] * P[r][(6, 6)] + mu2[1] * P[r][(5, 6)] + mu2[2] * P[r][(6, 5)]
            + mu2[3] * P[r][(5, 5)] + mu2[4] * P[r][(6, 4)] + mu2[5] * P[r][(5, 4)]
            + mu2[6] * P[r][(6, 2)] + mu2[7] * P[r][(6, 1)] + mu2[8] * P[r][(6, 0)]
        )
    unknowns7 = x54 + [x45[r] for r in range(r_start_45, n)]
    sol7 = _linsolve([sp.expand(e) for e in lyg], unknowns7)
    for r in range(n):
        P[r][(5, 4)] = sol7[x54[r]]
    for r in range(r_start_45, n):
        P[r][(4, 5)] = sol7[x45[r]]

    # ---- (S2) eq (8): circulant for p[r][4,4] (r=1..n-1); lyc[0]=sal ----
    mu3 = mu3_tab(c, tau)
    y44 = [None] + [sp.Symbol(f"y44_{r}") for r in range(1, n)]
    for r in range(1, n):
        P[r][(4, 4)] = y44[r]
    lyc = [None] * n
    for r in range(n):
        t = sf(r + 1)
        lyc[r] = P[t][(4, 4)] - (
            P[r][(4, 4)] + mu3[0] * P[r][(6, 6)] + mu3[1] * P[r][(5, 6)] + mu3[2] * P[r][(6, 5)]
            + mu3[3] * P[r][(5, 5)] + mu3[4] * P[r][(6, 4)] + mu3[5] * P[r][(5, 4)]
            + mu3[6] * P[r][(6, 2)] + mu3[7] * P[r][(6, 1)] + mu3[8] * P[r][(6, 0)]
            + mu3[9] * (P[r][(4, 2)] - P[t][(2, 4)]) + mu3[10] * (P[r][(4, 1)] - P[t][(1, 4)])
            + mu3[11] * (P[r][(4, 0)] - P[t][(0, 4)])
        )
    sol8 = _linsolve([sp.expand(lyc[r]) for r in range(1, n)], [y44[r] for r in range(1, n)])
    for r in range(1, n):
        P[r][(4, 4)] = sol8[y44[r]]
    # sal = lyc[0] re-evaluated with all data assigned
    sal = sp.expand(
        P[sf(1)][(4, 4)] - (
            P[0][(4, 4)] + mu3[0] * P[0][(6, 6)] + mu3[1] * P[0][(5, 6)] + mu3[2] * P[0][(6, 5)]
            + mu3[3] * P[0][(5, 5)] + mu3[4] * P[0][(6, 4)] + mu3[5] * P[0][(5, 4)]
            + mu3[6] * P[0][(6, 2)] + mu3[7] * P[0][(6, 1)] + mu3[8] * P[0][(6, 0)]
            + mu3[9] * (P[0][(4, 2)] - P[sf(1)][(2, 4)]) + mu3[10] * (P[0][(4, 1)] - P[sf(1)][(1, 4)])
            + mu3[11] * (P[0][(4, 0)] - P[sf(1)][(0, 4)])
        )
    )

    # ---- reduced degree sector curve eq (5)  script 259-263 ----
    for r in range(n):
        P[r][(6, 3)] = (
            sp.Rational(3, 4) * P[r][(6, 4)] + sp.Rational(3, 4) * P[r][(6, 2)] + sp.Rational(1, 20) * P[r][(6, 6)]
            - sp.Rational(3, 10) * P[r][(6, 5)] - sp.Rational(3, 10) * P[r][(6, 1)] + sp.Rational(1, 20) * P[r][(6, 0)]
        )
    # glue sectors: p[sf(r+1)][i,6] = p[r][6,i]  script 264-267
    for r in range(n):
        t = sf(r + 1)
        for i in range(7):
            P[t][(i, 6)] = P[r][(6, i)]

    # ---- (S3) eq (9)  script 269-315 ----
    mut = mut_tab(c, tau)
    mub = mub_tab(c, tau)
    p53 = {}
    p35 = {}
    p34 = {}
    for r in range(n):
        t = sf(r + 1)
        p53[r] = (
            mut[0] * P[r][(6, 6)] + mut[1] * P[r][(6, 5)] + mut[2] * P[r][(6, 4)] + mut[3] * P[r][(5, 4)]
            + mut[4] * P[r][(6, 2)] + mut[5] * P[r][(6, 1)] + mut[6] * P[r][(6, 0)]
            + mut[7] * (P[r][(4, 2)] - P[t][(2, 4)]) + mut[8] * (P[r][(4, 1)] - P[t][(1, 4)])
            + mut[9] * (P[r][(4, 0)] - P[t][(0, 4)]) + mut[10] * (P[r][(4, 4)] - P[t][(4, 4)])
        )
        p35[t] = (
            mut[0] * P[r][(6, 6)] + mut[1] * P[r][(6, 5)] + mut[2] * P[r][(6, 4)] + mut[3] * P[t][(4, 5)]
            + mut[4] * P[r][(6, 2)] + mut[5] * P[r][(6, 1)] + mut[6] * P[r][(6, 0)]
            - mut[7] * (P[r][(4, 2)] - P[t][(2, 4)]) - mut[8] * (P[r][(4, 1)] - P[t][(1, 4)])
            - mut[9] * (P[r][(4, 0)] - P[t][(0, 4)]) - mut[10] * (P[r][(4, 4)] - P[t][(4, 4)])
        )
    for r in range(n):
        t = sf(r + 1)
        p34[t] = (
            P[r][(4, 3)] + mub[0] * P[r][(6, 6)] + mub[1] * P[r][(6, 5)] + mub[2] * P[r][(6, 4)]
            + mub[3] * P[r][(5, 4)] + mub[4] * P[r][(6, 2)] + mub[5] * P[r][(6, 1)] + mub[6] * P[r][(6, 0)]
            + mub[7] * (P[r][(4, 2)] - P[t][(2, 4)]) + mub[8] * (P[r][(4, 1)] - P[t][(1, 4)])
            + mub[9] * (P[r][(4, 0)] - P[t][(0, 4)]) + mub[10] * (P[r][(4, 4)] - P[t][(4, 4)])
        )
    for r in range(n):
        P[r][(5, 3)] = p53[r]
    for t in range(n):
        P[t][(3, 5)] = p35[t]
        P[t][(3, 4)] = p34[t]

    if simplify:
        for s in range(n):
            for key in list(P[s].keys()):
                P[s][key] = sp.factor(P[s][key])

    return P, sal, {"d2": d2, "c": c, "tau": tau}


def _linsolve(eqs, unknowns):
    """Solve an (over/exactly)determined consistent linear system, return
    a dict unknown->expr.

    Exact (rational/symbolic) systems go through sympy.linsolve. Systems
    carrying floats are ill-suited to linsolve's exact pivoting, so they
    are solved by numpy least-squares on the (well-conditioned)
    coefficient matrix — appropriate for the numeric valence probe."""
    A, b = sp.linear_eq_to_matrix(eqs, unknowns)
    has_float = any(getattr(a, "is_Float", False) for a in A) or any(
        getattr(x, "is_Float", False) for x in b
    )
    if has_float:
        import numpy as np

        An = np.array(A.tolist(), dtype=float)
        bn = np.array(b.tolist(), dtype=float).reshape(-1)
        sol, *_ = np.linalg.lstsq(An, bn, rcond=None)
        return {u: sp.Float(val) for u, val in zip(unknowns, sol)}
    solset = sp.linsolve((A, b), unknowns)
    if not solset:
        raise RuntimeError("linear system inconsistent / no solution")
    tup = next(iter(solset))
    return {u: val for u, val in zip(unknowns, tup)}


# --------------------------------------------------------------------------
# Smoothness residuals G1[s], G2[s]  (script 318-353)
# --------------------------------------------------------------------------
def sector_residuals(P, s, n, c, tau):
    sf = lambda k: k % n
    t = sf(s + 1)
    C = sp.binomial
    F = 0
    G = 0
    for i in range(7):
        for j in range(3):
            F += P[s][(6 - j, 6 - i)] * C(6, i) * C(6, j) * (1 - U) ** (6 - i) * U**i * (1 - V) ** (6 - j) * V**j
            G += P[t][(6 - i, 6 - j)] * C(6, i) * C(6, j) * (1 - U) ** (6 - i) * U**i * (1 - V) ** (6 - j) * V**j
    at0 = lambda expr: expr.subs(V, 0)
    Fu = at0(sp.diff(F, U))
    Fv = at0(sp.diff(F, V))
    Gv = at0(sp.diff(G, V))
    Fuu = at0(sp.diff(F, U, U))
    Fuv = at0(sp.diff(F, U, V))
    Fvv = at0(sp.diff(F, V, V))
    Gvv = at0(sp.diff(G, V, V))
    g1 = Fv + Gv - 2 * c * (1 - U) ** 2 * Fu
    A = 24 * c * (tau - 1) / (1 - U + tau * U)
    g2 = (
        -Gvv + Fvv - 4 * c * (1 - U) ** 2 * Fuv + 4 * c**2 * (1 - U) ** 4 * Fuu
        + c * (1 - U) ** 3 * (A * U - 8 * c - A) * Fu + A * (1 - U) ** 2 * Fv
    )
    return g1, g2


# --------------------------------------------------------------------------
# Input builders
# --------------------------------------------------------------------------
def zero_net(n):
    return [[[0 for _ in range(4)] for _ in range(4)] for _ in range(n)]


def shared_net(n, rng):
    """C2-consistent net from per-spoke a1..a6 + global E (script comment
    lines 227-230).  net[s] laid out as:
        [a1[s], a2[s], a3[s], a4[s+1]]
        [a4[s], a5[s], a6[s], a5[s+1]]
        [a3[s-1], a6[s-1], E, a6[s+1]]
        [a2[s-1], a5[s-1], a6[s-2], 0]
    """

    def rv():
        return sp.Rational(rng.randint(-9, 9), rng.randint(1, 7))

    a = {name: [rv() for _ in range(n)] for name in ["a1", "a2", "a3", "a4", "a5", "a6"]}
    E = rv()
    sf = lambda k: k % n
    nets = []
    for s in range(n):
        net = [
            [a["a1"][s], a["a2"][s], a["a3"][s], a["a4"][sf(s + 1)]],
            [a["a4"][s], a["a5"][s], a["a6"][s], a["a5"][sf(s + 1)]],
            [a["a3"][sf(s - 1)], a["a6"][sf(s - 1)], E, a["a6"][sf(s + 1)]],
            [a["a2"][sf(s - 1)], a["a5"][sf(s - 1)], a["a6"][sf(s - 2)], 0],
        ]
        nets.append(net)
    return nets


def rational_free(n, rng):
    def rv():
        return sp.Rational(rng.randint(-9, 9), rng.randint(1, 7))

    return {
        "P00": rv(),
        "TAN0": rv(),
        "TAN1": rv(),
        "TWO0": rv(),
        "TWO1": rv(),
        "TWO2": rv(),
        "P220": rv(),
        "P210": rv(),
        "P23": [rv() for _ in range(n)],
        "P33": [rv() for _ in range(n)],
    }


def symbolic_free(n):
    S = sp.Symbol
    return {
        "P00": S("P00"),
        "TAN0": S("TAN0"),
        "TAN1": S("TAN1"),
        "TWO0": S("TWO0"),
        "TWO1": S("TWO1"),
        "TWO2": S("TWO2"),
        "P220": S("P220"),
        "P210": S("P210"),
        "P23": [S(f"P23_{r}") for r in range(n)],
        "P33": [S(f"P33_{r}") for r in range(n)],
    }


def c_of(n, exact):
    if exact:
        if n == 6:
            return sp.Rational(1, 2)
        if n == 3:
            return sp.Rational(-1, 2)
        if n == 4:
            return sp.Integer(0)
        return sp.cos(2 * sp.pi / n)
    return sp.Float(sp.cos(2 * sp.pi / sp.Integer(n)).evalf(30))


TAU_TABLE = {3: sp.Rational(4, 3), 5: sp.Rational(87, 100), 6: sp.Rational(17, 20),
             7: sp.Rational(21, 25), 8: sp.Rational(83, 100), 9: sp.Rational(41, 50), 10: sp.Rational(81, 100)}


# --------------------------------------------------------------------------
# CLI modes
# --------------------------------------------------------------------------
def run_symbolic(net_mode):
    n = 6
    c = sp.Rational(1, 2)
    tau = sp.Rational(17, 20)
    free = symbolic_free(n)
    if net_mode == "zero":
        net = zero_net(n)
    else:
        # symbolic shared net
        sf = lambda k: k % n
        a = {name: [sp.Symbol(f"{name}_{r}") for r in range(n)] for name in ["a1", "a2", "a3", "a4", "a5", "a6"]}
        E = sp.Symbol("E")
        net = []
        for s in range(n):
            net.append([
                [a["a1"][s], a["a2"][s], a["a3"][s], a["a4"][sf(s + 1)]],
                [a["a4"][s], a["a5"][s], a["a6"][s], a["a5"][sf(s + 1)]],
                [a["a3"][sf(s - 1)], a["a6"][sf(s - 1)], E, a["a6"][sf(s + 1)]],
                [a["a2"][sf(s - 1)], a["a5"][sf(s - 1)], a["a6"][sf(s - 2)], 0],
            ])
    print(f"[symbolic] n={n} c=1/2 tau=17/20 net={net_mode}  (building...)")
    P, sal, meta = build_cap(n, c, tau, net, free, simplify=False)
    print(f"  d2 = {sp.nsimplify(meta['d2'])} = {sp.factor(meta['d2'])}")
    print(f"  sal (deferred lyc[0]) simplifies to: {sp.simplify(sal)}")
    ok = True
    for s in range(n):
        g1, g2 = sector_residuals(P, s, n, c, tau)
        g1s = sp.simplify(g1)
        g2s = sp.simplify(sp.together(g2))
        z1 = g1s == 0
        z2 = g2s == 0
        ok = ok and z1 and z2
        print(f"  sector {s}: G1==0 {z1}   G2==0 {z2}")
        if not z1:
            print("    G1 =", g1s)
        if not z2:
            print("    G2 =", g2s)
    print(f"[symbolic] ACCEPTED: {ok and sp.simplify(sal) == 0}")
    return ok and sp.simplify(sal) == 0


def run_numeric(n, tau, seed, exact, net_mode, samples=None, verbose=True):
    rng = random.Random(seed)
    c = c_of(n, exact)
    tau = sp.nsimplify(tau) if exact else tau
    free = rational_free(n, rng)
    net = shared_net(n, rng) if net_mode == "shared" else zero_net(n)
    if not exact:
        # float everything
        c = sp.Float(c)
        tau = sp.Float(tau)
        free = {k: (sp.Float(v) if not isinstance(v, list) else [sp.Float(x) for x in v]) for k, v in free.items()}
        net = [[[sp.Float(x) for x in row] for row in blk] for blk in net]
    P, sal, meta = build_cap(n, c, tau, net, free, simplify=False)
    max_g1 = 0.0
    max_g2 = 0.0
    worst = None
    us = samples or [sp.Rational(i, 10) for i in range(1, 10)]
    sal_val = float(abs(sp.N(sal)))
    for s in range(n):
        g1, g2 = sector_residuals(P, s, n, c, tau)
        for uu in us:
            v1 = abs(float(sp.N(g1.subs(U, uu))))
            v2 = abs(float(sp.N(g2.subs(U, uu))))
            if v1 > max_g1:
                max_g1 = v1
            if v2 > max_g2:
                max_g2, worst = v2, (s, uu)
    if verbose:
        print(f"[numeric] n={n} tau={float(tau):.4f} c={float(c):.6f} net={net_mode} exact={exact} seed={seed}")
        print(f"  max|G1|={max_g1:.3e}  max|G2|={max_g2:.3e}  |sal|={sal_val:.3e}  worst_G2@{worst}")
    return max_g1, max_g2, sal_val


def run_valence():
    print("=== valence-generalization probe (float, shared net) ===")
    cases = [(5, TAU_TABLE[5]), (6, TAU_TABLE[6]), (7, TAU_TABLE[7]), (8, TAU_TABLE[8]), (9, TAU_TABLE[9]), (10, TAU_TABLE[10])]
    verdict = {}
    for n, tau in cases:
        g1, g2, sal = run_numeric(n, float(tau), seed=100 + n, exact=False, net_mode="shared")
        good = g1 < 1e-9 and g2 < 1e-8 and sal < 1e-9
        verdict[n] = (good, g1, g2, sal)
        print(f"  --> n={n}: {'CLOSES' if good else 'FAILS'}")
    print("VALENCE VERDICT:")
    for n, (good, g1, g2, sal) in verdict.items():
        print(f"  n={n:2d}: {'PASS' if good else 'FAIL'}  G1={g1:.2e} G2={g2:.2e} sal={sal:.2e}")
    return verdict


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="KP15 n6 Maple-script oracle")
    ap.add_argument("--mode", choices=["symbolic", "numeric", "valence", "all"], default="all")
    ap.add_argument("--net", choices=["zero", "shared"], default="shared")
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    if args.mode in ("symbolic", "all"):
        run_symbolic(args.net if args.mode == "symbolic" else "zero")
    if args.mode in ("numeric", "all"):
        tau = float(TAU_TABLE.get(args.n, sp.Rational(4, 5)))
        run_numeric(args.n, tau, args.seed, exact=(args.n in (3, 4, 6)), net_mode=args.net)
        # exact rational proof at n=6
        if args.n == 6:
            print("[numeric-exact] n=6 exact-rational coefficient check:")
            g1max, g2max, sal = run_numeric(6, sp.Rational(17, 20), 7, exact=True, net_mode="shared", verbose=True)
    if args.mode in ("valence", "all"):
        run_valence()
