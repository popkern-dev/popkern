#!/usr/bin/env python3
"""Clean-room generator for the PnS T-junction mask tables (PK-317).

Implements the GT-spline constructions of Karciauskas & Peters, "Low
degree splines for locally quad-dominant meshes", CAGD 83 (2020)
(docs/papers/20quaT.pdf). T0 here; T1/T2 follow.

T0: the tau0-net has 14 nodes — two SPARSE rows of 3 (s0 = [n1 n2 n3]
outermost, s1 = [n4 n5 n6]) and two DENSE rows of 4 (d1 = [n7..n10],
d0 = [n11..n14] outermost), with in-row node order OPPOSITE the mask
grid's u direction (measured layout; ordering is a package
convention, values are paper-derived). The tau-surface is a 2x2
macro-patch of bi-3 pieces = a 7x7 grid of shared BB-coefficients:

  rows r0, r1  : dense-side tensor-border + the M_b stencil row
  rows r5, r6  : sparse-side tensor-border (tensor of a 3-node
                 degree-elevated bi-2 u-profile with v-weights
                 (1/6,5/6), (1/2,1/2)) with the M_t stencil filling
                 the r5 interior
  row  r2, r4  : horizontal C2-layers — ends from the side
                 tensor-borders and the printed M_m stencils, centre
                 from the printed 7-weight /96 rows, interior from
                 the C2-completion solve c2 = (c1 + 4m - g2)/4,
                 g1 = 2m - c2
  row  r3      : the internal C1 mid-layer = average of r2 and r4

Every entry is an exact rational; the package table is a verification
oracle only.
"""

from fractions import Fraction as F

import numpy as np

COLS = 14

# Node indices (0-based) for the four net rows, in mask-grid u order
# (left to right).
S0 = [2, 1, 0]        # n3 n2 n1
S1 = [5, 4, 3]        # n6 n5 n4
D1 = [9, 8, 7, 6]     # n10 n9 n8 n7
D0 = [13, 12, 11, 10]


def vec(weights):
    v = [F(0)] * COLS
    for idx, w in weights:
        v[idx] += w
    return v


def combine(rows, coeffs):
    v = [F(0)] * COLS
    for row, c in zip(rows, coeffs):
        for i in range(COLS):
            v[i] += c * row[i]
    return v


def tensor(u_profile, nodes_rows, v_weights):
    """u_profile over len(nodes_rows[0]) nodes x v blend over rows."""
    v = [F(0)] * COLS
    for row_nodes, vw in zip(nodes_rows, v_weights):
        for node, uw in zip(row_nodes, u_profile):
            v[node] += vw * uw
    return v


# Degree-elevated bi-2 boundary u-profiles at the 7 BB positions.
# 4-node (dense) side:
U_DENSE = [
    [F(1, 4), F(1, 4) * 0, 0, 0],  # placeholder; replaced below
]
# Explicit profiles (derived from bi-2 B-to-BB + raise to degree 3,
# verified exact):
U_DENSE = [
    [F(1, 2), F(1, 2), F(0), F(0)],
    [F(1, 6), F(5, 6) * F(1, 1), F(0), F(0)],  # placeholder fixed below
]


def u_dense_profiles():
    # positions c0..c6 over 4 uniform bi-2 control values expressed as
    # the two bi-2 pieces' BB rows raised to shared degree-3 layout:
    # c0 corner (a0+a1)/2 ; c1 (a0+5a1)/6 ; c2 (a1*5+a2)/6 ... measured
    # exact forms:
    return [
        [F(1, 2), F(1, 2), F(0), F(0)],
        [F(1, 6), F(5, 6), F(0), F(0)],
        [F(0), F(5, 6), F(1, 6), F(0)],
        [F(0), F(1, 2), F(1, 2), F(0)],
        [F(0), F(1, 6), F(5, 6), F(0)],
        [F(0), F(0), F(5, 6), F(1, 6)],
        [F(0), F(0), F(1, 2), F(1, 2)],
    ]


def u_sparse_profiles():
    return [
        [F(1, 2), F(1, 2), F(0)],
        [F(1, 3), F(2, 3), F(0)],
        [F(5, 24), F(3, 4), F(1, 24)],
        [F(1, 8), F(3, 4), F(1, 8)],
        [F(1, 24), F(3, 4), F(5, 24)],
        [F(0), F(2, 3), F(1, 3)],
        [F(0), F(1, 2), F(1, 2)],
    ]


def t0_grid():
    """The 7x7 grid of mask rows (each a COLS-vector of Fractions)."""
    g = [[None] * 7 for _ in range(7)]

    # r0: dense tensor-border boundary (v-blend 1/2,1/2 of d1,d0).
    # r1: M_b row layer (v-blend 5/6,1/6 at the plain-tensor columns,
    #     printed M_b stencils at c1, c2, c4, c5).
    for c, up in enumerate(u_dense_profiles()):
        g[0][c] = tensor(up, [D1, D0], [F(1, 2), F(1, 2)])
        g[1][c] = tensor(up, [D1, D0], [F(5, 6), F(1, 6)])
    # printed M_b stencils overwrite the r1 interior neighbours:
    #   [3 11 9 49]/72 at c1 (and mirrored at c5),
    #   [11 51 1 9]/72 at c2 (mirrored at c4)
    g[1][1] = vec([(D1[0], F(11, 72)), (D1[1], F(49, 72)), (D0[0], F(3, 72)), (D0[1], F(9, 72))])
    g[1][5] = vec([(D1[3], F(11, 72)), (D1[2], F(49, 72)), (D0[3], F(3, 72)), (D0[2], F(9, 72))])
    g[1][2] = vec([(D1[1], F(51, 72)), (D1[2], F(9, 72)), (D0[1], F(11, 72)), (D0[2], F(1, 72))])
    g[1][4] = vec([(D1[2], F(51, 72)), (D1[1], F(9, 72)), (D0[2], F(11, 72)), (D0[1], F(1, 72))])

    # r6 / r5: sparse tensor-border with M_t stencils in the r5
    # interior: [3 19 9 41]/72 at c1/c5, [2 12 9 45 1 3]/72 at c2/c4.
    for c, up in enumerate(u_sparse_profiles()):
        g[6][c] = tensor(up, [S1, S0], [F(1, 2), F(1, 2)])
        g[5][c] = tensor(up, [S1, S0], [F(5, 6), F(1, 6)])
    g[5][1] = vec([(S0[0], F(3, 72)), (S0[1], F(9, 72)), (S1[0], F(19, 72)), (S1[1], F(41, 72))])
    g[5][5] = vec([(S0[2], F(3, 72)), (S0[1], F(9, 72)), (S1[2], F(19, 72)), (S1[1], F(41, 72))])
    g[5][2] = vec([
        (S0[0], F(2, 72)), (S0[1], F(9, 72)), (S0[2], F(1, 72)),
        (S1[0], F(12, 72)), (S1[1], F(45, 72)), (S1[2], F(3, 72)),
    ])
    g[5][4] = vec([
        (S0[2], F(2, 72)), (S0[1], F(9, 72)), (S0[0], F(1, 72)),
        (S1[2], F(12, 72)), (S1[1], F(45, 72)), (S1[0], F(3, 72)),
    ])

    # r2 / r4 horizontal C2 layers.
    mid3 = [3, 4, 5]          # n4 n5 n6 (paper order)
    mid4 = [6, 7, 8, 9]       # n7 n8 n9 n10
    def w96(three, four):
        return vec(
            [(n, F(t, 96)) for n, t in zip(mid3, three)]
            + [(n, F(t, 96)) for n, t in zip(mid4, four)]
        )
    center_r2 = w96((-1, 18, -1), (3, 37, 37, 3))
    center_r4 = w96((7, 66, 7), (3, 5, 5, 3))

    # ends: c0/c6 boundary tensor (edge rule between the two middle
    # net columns) and the printed M_m stencils at c1/c5:
    #   r2: [1 1 5 5]/12 boundary, [7 23 1 5]/36
    #   r4: [5 5 1 1]/12 boundary, [1 2 4 11]/18
    g[2][0] = vec([(4, F(1, 12)), (5, F(1, 12)), (8, F(5, 12)), (9, F(5, 12))])
    g[2][6] = vec([(3, F(1, 12)), (4, F(1, 12)), (6, F(5, 12)), (7, F(5, 12))])
    g[2][1] = vec([(4, F(5, 36)), (5, F(1, 36)), (8, F(23, 36)), (9, F(7, 36))])
    g[2][5] = vec([(4, F(5, 36)), (3, F(1, 36)), (7, F(23, 36)), (6, F(7, 36))])
    g[4][0] = vec([(4, F(5, 12)), (5, F(5, 12)), (8, F(1, 12)), (9, F(1, 12))])
    g[4][6] = vec([(3, F(5, 12)), (4, F(5, 12)), (6, F(1, 12)), (7, F(1, 12))])
    g[4][1] = vec([(4, F(22, 36)), (5, F(8, 36)), (8, F(4, 36)), (9, F(2, 36))])
    g[4][5] = vec([(4, F(22, 36)), (3, F(8, 36)), (7, F(4, 36)), (6, F(2, 36))])
    g[2][3] = center_r2
    g[4][3] = center_r4

    # C2-completion for the interiors of r2 and r4:
    #   c2 = (c1 + 4 m - g2)/4,  c4(=g1) = 2 m - c2  (per row).
    for r in (2, 4):
        c1, m, g2 = g[r][1], g[r][3], g[r][5]
        c2 = combine([c1, m, g2], [F(1, 4), F(1, 1), F(-1, 4)])
        g[r][2] = c2
        g[r][4] = combine([m, c2], [F(2), F(-1)])

    # r3: internal C1 mid-layer.
    for c in range(7):
        g[3][c] = combine([g[2][c], g[4][c]], [F(1, 2), F(1, 2)])

    return g


def t0_mask():
    g = t0_grid()
    # patches: P0 top-right, P1 bottom-right, P2 top-left,
    # P3 bottom-left; each 4x4 in (row, col) grid order.
    spans = [((0, 4), (3, 7)), ((3, 7), (3, 7)), ((0, 4), (0, 4)), ((3, 7), (0, 4))]
    rows = []
    for (r0, r1), (c0, c1) in spans:
        for r in range(r0, r1):
            for c in range(c0, c1):
                rows.append([float(x) for x in g[r][c]])
    return np.array(rows)


def main():
    import sys
    sys.path.insert(0, 'tools/reference/pns-tables')
    from polar import load_package_csv
    ref = load_package_csv(
        "docs/code/ACM-ToMS-PnS-1032/1032/Src/Source/src/Patch/Table/T0.csv"
    )
    gen = t0_mask()
    if gen.shape != ref.shape:
        print(f"T0 SHAPE {gen.shape} vs {ref.shape}")
        sys.exit(1)
    dev = float(np.abs(gen - ref).max())
    print(f"T0: shape {gen.shape} max|dev| {dev:.3e}")
    sys.exit(0 if dev < 2e-9 else 1)


if __name__ == "__main__":
    main()
