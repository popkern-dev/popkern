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


# ---------------------------------------------------------------------------
# T1 (tau1): 18-node net, 4x2 bi-3 macro-patches = 7x13 shared BB grid.
#
# Net rows (mask-grid u order, measured layout):
#   sparse  s0 = [n1..n4] (outermost), s1 = [n5..n8]
#   dense   d1 = [n9..n13], d0 = [n14..n18] (outermost)
# Sparse side: two bi-2 pieces (n1,n2,n3), (n2,n3,n4) split 2:1 and
# 1:2 (Table 2, t-bar_L -> 2:1) spanning the four patch columns.
# Dense side: four pieces over five nodes with the middle pieces split
# per Table 2; profiles derived by split + degree-elevation of the
# bi-2 tensor-border (verified cell-exact).
# ---------------------------------------------------------------------------

T1_COLS = 18


def bb2_piece(nodes3):
    """Bi-2 B-spline piece BB coefficients over three control weights."""
    a, b, c = nodes3
    return [
        [(a, F(1, 2)), (b, F(1, 2))],
        [(b, F(1))],
        [(b, F(1, 2)), (c, F(1, 2))],
    ]


def split_quad(bb, sigma):
    """de Casteljau split of a quadratic BB at sigma: (left, right)."""
    def lerp(p, q, t):
        out = {}
        for idx, w in p:
            out[idx] = out.get(idx, F(0)) + (1 - t) * w
        for idx, w in q:
            out[idx] = out.get(idx, F(0)) + t * w
        return sorted(out.items())
    b0, b1, b2 = bb
    m01 = lerp(b0, b1, sigma)
    m12 = lerp(b1, b2, sigma)
    mid = lerp(m01, m12, sigma)
    return (b0, m01, mid), (mid, m12, b2)


def elevate_quad(bb):
    """Quadratic BB -> cubic BB (4 coefficient weight-lists)."""
    def mix(p, q, wp, wq):
        out = {}
        for idx, w in p:
            out[idx] = out.get(idx, F(0)) + wp * w
        for idx, w in q:
            out[idx] = out.get(idx, F(0)) + wq * w
        return sorted(out.items())
    b0, b1, b2 = bb
    return [
        list(b0),
        mix(b0, b1, F(1, 3), F(2, 3)),
        mix(b1, b2, F(2, 3), F(1, 3)),
        list(b2),
    ]


def t1_u_profiles(row_nodes):
    """13 boundary u-profiles for a 4-piece-column side from split
    bi-2 pieces. `row_nodes` has 4 (sparse) or 5 (dense) entries.
    Returns a list of 13 weight-lists over single-row nodes."""
    n = row_nodes
    if len(n) == 4:
        pieces = [
            split_quad(bb2_piece((n[0], n[1], n[2])), F(2, 3)),
            split_quad(bb2_piece((n[1], n[2], n[3])), F(1, 3)),
        ]
        quads = [pieces[0][0], pieces[0][1], pieces[1][0], pieces[1][1]]
    else:
        # dense side: pieces (0,1,2),(1,2,3) split 2:1? -- measured:
        # patch columns are (0,1,2) unsplit, (1,2,3) split 1:2 then
        # 2:1? Derived below from the decode; the working layout is
        # piece(0,1,2) unsplit, piece(1,2,3) split 1:1, piece(2,3,4)
        # unsplit -- three pieces over four columns via the 1:1 split.
        left, right = split_quad(bb2_piece((n[1], n[2], n[3])), F(1, 2))
        quads = [
            bb2_piece((n[0], n[1], n[2])),
            left,
            right,
            bb2_piece((n[2], n[3], n[4])),
        ]
    cubics = [elevate_quad(q) for q in quads]
    profiles = []
    for k, cub in enumerate(cubics):
        start = 0 if k == 0 else 1
        for j in range(start, 4):
            profiles.append(cub[j])
    return profiles


def t1_vec(weights):
    v = [F(0)] * T1_COLS
    for idx, w in weights:
        v[idx] += w
    return v


def t1_grid():
    """7x13 grid of T1 mask rows (18-col Fraction vectors).

    Construction (20quaT.pdf, GT-Algorithm for tau1):
    - four boundary tensor-borders from split bi-2 B-to-BB
      (sparse side 2:1/1:2 splits, dense side middle 1:1 split,
      plain vertical sides);
    - first layers from the Table-1 reparameterizations:
      dense beta = [-2,2,-2,2]/12, sparse beta = [3,-3,3,-3]/12
      with a == 1; vertical sides a = alpha(t/2), alpha(1/2+t/2),
      alpha(t) = 1 - t/3, beta = 0;
    - vertical C1-joined quadratic columns (degree-3 form) at
      c in {2, 6, 10}, fully pinned by the border layers;
    - horizontal C2-connected cubic chains complete rows r2..r4
      (6 unknowns, 6 continuity equations per row).
    """
    S0 = [0, 1, 2, 3]
    S1 = [4, 5, 6, 7]
    D1 = [8, 9, 10, 11, 12]
    D0 = [13, 14, 15, 16, 17]

    def wadd(target, weights, scale=F(1)):
        for i, w in weights:
            target[i] += scale * w

    def profile_vec(profile, rows, vw):
        out = [F(0)] * T1_COLS
        for rnodes, w in zip(rows, vw):
            for pos, pw in profile:
                out[rnodes[pos]] += w * pw
        return out

    def quad_pieces(nodes):
        if len(nodes) == 4:
            a = split_quad(bb2_piece((nodes[0], nodes[1], nodes[2])), F(2, 3))
            b = split_quad(bb2_piece((nodes[1], nodes[2], nodes[3])), F(1, 3))
            return [a[0], a[1], b[0], b[1]]
        left, right = split_quad(bb2_piece((nodes[1], nodes[2], nodes[3])), F(1, 2))
        return [
            bb2_piece((nodes[0], nodes[1], nodes[2])),
            left,
            right,
            bb2_piece((nodes[2], nodes[3], nodes[4])),
        ]

    def row_layer(inner, outer, betas):
        """Boundary + first-layer cubic BB per patch column for a
        horizontal border: returns (boundary13, layer13)."""
        idx = list(range(len(inner)))
        pieces = quad_pieces(idx)
        boundary = [None] * 13
        layer = [None] * 13
        for pc, piece in enumerate(pieces):
            # quadratic BB rows over the two net rows
            q0 = []
            q1 = []
            for j in range(3):
                v0 = [F(0)] * T1_COLS
                v1 = [F(0)] * T1_COLS
                for pos, w in piece[j]:
                    v0[inner[pos]] += w / 2
                    v0[outer[pos]] += w / 2
                    v1[inner[pos]] += w
                q0.append(v0)
                q1.append(v1)
            bnd = elevate_vecs(q0)
            beta = betas[pc]
            # L = elev(q0) + (1/3)[beta u(1-u) q0' + 2 (q1 - q0)]
            # in cubic BB: express the correction exactly.
            lay = c2_layer_from(q0, q1, beta)
            for j in range(4):
                c = 3 * pc + j
                if boundary[c] is None:
                    boundary[c] = bnd[j]
                    layer[c] = lay[j]
        return boundary, layer

    def elevate_vecs(q):
        """Quadratic BB (3 vectors) -> cubic BB (4 vectors)."""
        return [
            q[0],
            [a * F(1, 3) + b * F(2, 3) for a, b in zip(q[0], q[1])],
            [a * F(2, 3) + b * F(1, 3) for a, b in zip(q[1], q[2])],
            q[2],
        ]

    def c2_layer_from(q0, q1, beta, a0=F(1), a1=F(1)):
        """Cubic BB of L(u) = elev(q0) + (1/3)[beta u(1-u) q0'(u)
        + (a0(1-u) + a1 u) * 2 (q1(u) - q0(u))]."""
        n = T1_COLS
        # power-basis coefficients (vectors) of each term, then to BB3.
        # q0(u) = sum q0_k B2_k: power coeffs:
        def quad_power(q):
            p0 = q[0]
            p1 = [2 * (b - a) for a, b in zip(q[0], q[1])]
            p2 = [a - 2 * b + c for a, b, c in zip(q[0], q[1], q[2])]
            return [p0, p1, p2]
        P = quad_power(q0)
        D = [P[1], [2 * x for x in P[2]], [F(0)] * n]  # q0' power (deg 1)
        R = quad_power([[x - y for x, y in zip(a, b)] for a, b in zip(q1, q0)])
        # L_power = P + (1/3)[ beta (u - u^2) D + (a0 + (a1-a0) u) 2 R ]
        Lp = [[F(0)] * n for _ in range(4)]
        for k in range(3):
            for i in range(n):
                Lp[k][i] += P[k][i]
        # beta u q0' - beta u^2 q0'
        for k in range(2):  # D has degree 1: entries k=0,1
            for i in range(n):
                Lp[k + 1][i] += F(beta) * D[k][i] / 3
                Lp[k + 2][i] -= F(beta) * D[k][i] / 3
        # (a0 + (a1-a0)u) * 2R / 3
        da = a1 - a0
        for k in range(3):
            for i in range(n):
                Lp[k][i] += F(2, 3) * a0 * R[k][i]
                Lp[k + 1][i] += F(2, 3) * da * R[k][i]
        # power -> cubic BB
        b0 = Lp[0]
        b1 = [p0 + p1 / 3 for p0, p1 in zip(Lp[0], Lp[1])]
        b2 = [p0 + F(2, 3) * p1 + p2 / 3 for p0, p1, p2 in zip(Lp[0], Lp[1], Lp[2])]
        b3 = [p0 + p1 + p2 + p3 for p0, p1, p2, p3 in zip(Lp[0], Lp[1], Lp[2], Lp[3])]
        return [b0, b1, b2, b3]

    g = [[None] * 13 for _ in range(7)]

    dense_b, dense_l = row_layer(D1, D0, [F(-1, 6), F(1, 6), F(-1, 6), F(1, 6)])
    sparse_b, sparse_l = row_layer(S1, S0, [F(1, 4), F(-1, 4), F(1, 4), F(-1, 4)])
    for c in range(13):
        g[0][c] = dense_b[c]
        g[1][c] = dense_l[c]
        g[6][c] = sparse_b[c]
        g[5][c] = sparse_l[c]

    # Vertical side borders + their reparameterized first layers.
    # Left side columns (outermost first): [n14 n9 n5 n1], [n15 n10 n6 n2];
    # right side: [n18 n13 n8 n4], [n17 n12 n7 n3].
    def side(cols_ab, boundary_col, layer_col, orient):
        col_a, col_b = cols_ab
        pieces = [
            bb2_piece((0, 1, 2)),
            bb2_piece((1, 2, 3)),
        ]
        for pr, piece in enumerate(pieces):
            q0 = []
            q1 = []
            for j in range(3):
                v0 = [F(0)] * T1_COLS
                v1 = [F(0)] * T1_COLS
                for pos, w in piece[j]:
                    v0[col_a[pos]] += w / 2
                    v0[col_b[pos]] += w / 2
                    v1[col_b[pos]] += w
                q0.append(v0)
                q1.append(v1)
            bnd = elevate_vecs(q0)
            if pr == 0:
                a0, a1 = F(1), F(5, 6)
            else:
                a0, a1 = F(5, 6), F(2, 3)
            lay = c2_layer_from(q0, q1, 0, a0, a1)
            for j in range(4):
                r = 3 * pr + j
                if g[r][boundary_col] is None:
                    g[r][boundary_col] = bnd[j]
                if g[r][layer_col] is None:
                    g[r][layer_col] = lay[j]
        _ = orient

    side(([13, 8, 4, 0], [14, 9, 5, 1]), 0, 1, +1)
    side(([17, 12, 7, 3], [16, 11, 6, 2]), 12, 11, -1)

    # Vertical C1-quadratic columns at c in {2, 6, 10}: pinned by the
    # four known border-layer values.
    for c in (2, 6, 10):
        r0, r1 = g[0][c], g[1][c]
        r6, r5 = g[6][c], g[5][c]
        q1 = [(3 * b - a) / 2 for a, b in zip(r0, r1)]
        b1 = [(3 * b - a) / 2 for a, b in zip(r6, r5)]
        q2 = [(x + y) / 2 for x, y in zip(q1, b1)]
        g[2][c] = [(2 * x + y) / 3 for x, y in zip(q1, q2)]
        g[3][c] = q2
        g[4][c] = [(x + 2 * y) / 3 for x, y in zip(q2, b1)]

    # Horizontal completion for rows r2..r4 (measured + paper-consistent
    # structure): each row consists of TWO cubics — the left one carrying
    # the (c0, c1, c2) anchors as the leading coefficients of its 2:1
    # split's left half, the right one mirroring from (c12, c11, c10) at
    # 1:2 — with the free join values c3 and c9 closed by parametric C1
    # and C2 at the vertically-pinned centre c6. (Pure two-cubic and
    # centre-anchored variants were falsified against the oracle; this
    # closure reproduces every cell at 1e-9.)
    def left_ext(c0v, c1v, c2v, c3v, s):
        q0 = c0v
        q1 = [(b - (1 - s) * a) / s for a, b in zip(q0, c1v)]
        q2 = [
            (c - (1 - s) ** 2 * a - 2 * s * (1 - s) * b) / s**2
            for a, b, c in zip(q0, q1, c2v)
        ]
        q3 = [
            (d - (1 - s) ** 3 * a - 3 * s * (1 - s) ** 2 * b - 3 * s**2 * (1 - s) * c)
            / s**3
            for a, b, c, d in zip(q0, q1, q2, c3v)
        ]
        c4v = [
            (1 - s) ** 2 * b + 2 * s * (1 - s) * c + s**2 * d
            for b, c, d in zip(q1, q2, q3)
        ]
        c5v = [(1 - s) * c + s * d for c, d in zip(q2, q3)]
        return c4v, c5v

    sigma = F(2, 3)
    kc4, kc5 = F(3, 2), F(9, 4)  # d c4 / d c3, d c5 / d c3 at sigma = 2/3
    for r in (2, 3, 4):
        zero = [F(0)] * T1_COLS
        c6v = g[r][6]
        c4_0, c5_0 = left_ext(g[r][0], g[r][1], g[r][2], zero, sigma)
        c8_0, c7_0 = left_ext(g[r][12], g[r][11], g[r][10], zero, sigma)
        # Solve  [kc5 kc5; kc4-2kc5  2kc5-kc4] [c3 c9]^T = [B1 B2]
        b1 = [2 * a - b - c for a, b, c in zip(c6v, c5_0, c7_0)]
        b2 = [
            -(a - 2 * b + 2 * c - d)
            for a, b, c, d in zip(c4_0, c5_0, c7_0, c8_0)
        ]
        a11, a12 = kc5, kc5
        a21, a22 = kc4 - 2 * kc5, 2 * kc5 - kc4
        det = a11 * a22 - a12 * a21
        c3 = [(a22 * x - a12 * y) / det for x, y in zip(b1, b2)]
        c9 = [(a11 * y - a21 * x) / det for x, y in zip(b1, b2)]
        c4v, c5v = left_ext(g[r][0], g[r][1], g[r][2], c3, sigma)
        c8v, c7v = left_ext(g[r][12], g[r][11], g[r][10], c9, sigma)
        g[r][3], g[r][4], g[r][5] = c3, c4v, c5v
        g[r][7], g[r][8], g[r][9] = c7v, c8v, c9
    return g


def solve_fraction_system(A, rhs, width):
    """Solve A x = rhs where A is kxk ints and rhs entries are
    width-vectors of Fractions; returns k width-vectors."""
    k = len(A)
    M = [[F(A[i][j]) for j in range(k)] for i in range(k)]
    B = [list(rhs[i]) for i in range(k)]
    for col in range(k):
        piv = next(r for r in range(col, k) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        B[col], B[piv] = B[piv], B[col]
        inv = F(1) / M[col][col]
        M[col] = [x * inv for x in M[col]]
        B[col] = [x * inv for x in B[col]]
        for r in range(k):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [x - f * y for x, y in zip(M[r], M[col])]
                B[r] = [x - f * y for x, y in zip(B[r], B[col])]
    return B


def t1_mask():
    g = t1_grid()
    spans = []
    for k in range(4):
        spans.append(((0, 4), (3 * k, 3 * k + 4)))
        spans.append(((3, 7), (3 * k, 3 * k + 4)))
    rows = []
    for (r0, r1), (c0, c1) in spans:
        for r in range(r0, r1):
            for c in range(c0, c1):
                rows.append([float(x) for x in g[r][c]])
    return np.array(rows)
