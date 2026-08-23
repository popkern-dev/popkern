#!/usr/bin/env python3
"""Clean-room generator for the PnS polar mask tables (PK-317).

Implements the bi-degree (2,3) polar cap of Karciauskas & Peters,
"Smooth polar caps for locally quad-dominant meshes", CAGD 81 (2020)
(docs/papers/20polar.pdf), Section 4: one patch of circular degree 2
and radial degree 3 per sector, with the radial degree-3 layer
collapsed to the pole.

The construction is fully closed-form (no tuned parameters):

  pole      C        = 3/4 P + 1/(4m) sum_i p_i
  outer     b^i_00   = 1/4(p_{i-1} + p_i) + 1/2 P
            b^i_10   = 1/2 p_i + 1/2 P
            b^i_20   = 1/4(p_i + p_{i+1}) + 1/2 P
  next      b^i_j1   = 1/3 P + 2/3 b^i_j0
  tangent   b^i_02   = C + 1/(3m) sum_s c_s 1/2 (p_{i-1+s} + p_{i+s})
            b^i_22   = b^{i+1}_02
            b^i_12   = (b^i_02 + b^i_22)/(c+1) + (c-1)/(c+1) C
  (c_s = cos(2 pi s / m), c = c_1)

Each mask row is a coefficient vector over the stencil
[P, p_0, ..., p_{m-1}]; every row is affine (sums to 1). Row layout
(measured against the package tables, which serve as verification
oracles only): per sector, radial layers from the POLE outward
(k = 3, 2, 1, 0), each layer j = 0, 1, 2 in circular order.
"""

import math

import numpy as np


def polar_mask(m: int) -> np.ndarray:
    """(12*m) x (m+1) mask matrix for a valence-m polar cap."""
    cols = m + 1  # [P, p_0..p_{m-1}]

    def vec(p_coef=None, big_p=0.0):
        v = np.zeros(cols)
        v[0] = big_p
        if p_coef:
            for idx, w in p_coef.items():
                v[1 + (idx % m)] += w
        return v

    c = math.cos(2.0 * math.pi / m)

    pole = vec({i: 1.0 / (4.0 * m) for i in range(m)}, big_p=0.75)

    def b0(i, j):
        if j == 0:
            return vec({i - 1: 0.25, i: 0.25}, big_p=0.5)
        if j == 1:
            return vec({i: 0.5}, big_p=0.5)
        return vec({i: 0.25, i + 1: 0.25}, big_p=0.5)

    def b1(i, j):
        return vec(big_p=1.0 / 3.0) + (2.0 / 3.0) * b0(i, j)

    def b02(i):
        v = pole.copy()
        for s in range(m):
            cs = math.cos(2.0 * math.pi * s / m)
            v += (cs / (3.0 * m)) * 0.5 * (vec({i - 1 + s: 1.0}) + vec({i + s: 1.0}))
        return v

    def b2(i, j):
        if j == 0:
            return b02(i)
        if j == 2:
            return b02(i + 1)
        return (b02(i) + b02(i + 1)) / (c + 1.0) + ((c - 1.0) / (c + 1.0)) * pole

    rows = []
    for i in range(m):
        for k in (3, 2, 1, 0):
            for j in range(3):
                if k == 3:
                    rows.append(pole)
                elif k == 2:
                    rows.append(b2(i, j))
                elif k == 1:
                    rows.append(b1(i, j))
                else:
                    rows.append(b0(i, j))
    return np.array(rows)


def load_package_csv(path: str) -> np.ndarray:
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            rows.append([float(x) for x in line.split(",")])
    return np.array(rows)


def main() -> None:
    import sys

    pkg_dir = "docs/code/ACM-ToMS-PnS-1032/1032/Src/Source/src/Patch/Table"
    worst = 0.0
    for m in range(3, 9):
        gen = polar_mask(m)
        ref = load_package_csv(f"{pkg_dir}/polarSct{m}.csv")
        if gen.shape != ref.shape:
            print(f"m={m}: SHAPE {gen.shape} vs {ref.shape}")
            sys.exit(1)
        dev = float(np.abs(gen - ref).max())
        worst = max(worst, dev)
        print(f"m={m}: shape {gen.shape} max|dev| {dev:.3e}")
    print(f"worst over all m: {worst:.3e}")
    # The package tables carry ~5e-10 of their OWN numeric noise (they
    # were produced by numeric solves and printed at 10 digits; e.g.
    # polarSct3 prints 0.1666666661 where the closed form is exactly
    # 1/6) — the generator's exact closed forms deviate only within
    # that noise.
    sys.exit(0 if worst < 2e-9 else 1)


if __name__ == "__main__":
    main()
