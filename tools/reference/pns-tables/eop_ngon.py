"""Clean-room regeneration of the PnS ngon/eop mask tables (star and n-gon caps).

Source: J. Fan, J. Peters, "Smooth multi-sided blending of biquadratic
splines" (docs/papers/14bi234.pdf) — Section 4.2 (bi-3 caps, n in {3,5}),
Appendix C (2x2 split bi-3 macro caps, n >= 6), Appendix A/Fig. 5 (CC-net
virtual refinement), eq. (5) central point, eq. (1)/(3^3) G1 join and border
reparameterization.

Input configurations (paper Fig. 4/6):
  ngon: DS-net — 4n points, sector-major groups of 4 (this package layout,
        decoded from the tables): within sector k the points are
        d0 (diagonal), d1, d2 (edge-adjacent), d3 (central n-gon vertex);
        sector face is the quad {d0,d1,d2,d3} (cyclic d2-d3-d1-d0), the
        central face is {d^k_3}_k, and the face across the original edge
        between sectors k-1,k is {d^k_2, d^k_3, d^{k-1}_1, d^{k-1}_3}.
  eop:  CC-net — 2n+1 points (E + one layer); equals ngon composed with the
        bi-2 virtual refinement stencil (9,3,3,1)/16 (paper Fig. 5/6a,
        Appendix A h-bar combination).

Cap layout (decoded package row order): row = 16*s + 4*a + b for n in {3,5}
(one bi-3 patch per sector), with (a,b)=(0,0) the outer corner over the
sector face centre, b-edge running to the corner shared with sector k-1,
a-edge to the corner shared with sector k+1, and (3,3) the cap centre.

The package tables were computed with ~1e-9 of solver/print noise, so the
verification gate is 2e-9 (same rationale as the polar/T tables).
"""

import numpy as np

BASE = 'docs/code/ACM-ToMS-PnS-1032/1032/Src/Source/src/Patch/Table/'


def load_package_csv(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            rows.append([float(x) for x in line.split(',')])
    return np.array(rows)


# ---------------------------------------------------------------------------
# small BB helpers (quadratic/cubic Bezier arithmetic on stencil vectors)
# ---------------------------------------------------------------------------

def elevate_quad(q):
    """Quadratic BB (3 stencil rows) -> cubic BB (4 rows)."""
    q0, q1, q2 = q
    return [q0, (q0 + 2 * q1) / 3, (2 * q1 + q2) / 3, q2]


def lin_times_quad(l, q):
    """(linear BB l) * (quadratic BB q) -> cubic BB."""
    l0, l1 = l
    q0, q1, q2 = q
    return [l0 * q0,
            (2 * l0 * q1 + l1 * q0) / 3,
            (2 * l1 * q1 + l0 * q2) / 3,
            l1 * q2]


# ---------------------------------------------------------------------------
# single bi-3 patch per sector (Section 4.2, n in {3,5})
# ---------------------------------------------------------------------------

def ngon_small(n):
    """Generating-function table (16n x 4n) for the bi-3 n-gon cap, n odd."""
    c = np.cos(2 * np.pi / n)
    m = 4 * n

    def d(k, j):
        v = np.zeros(m)
        v[4 * (k % n) + j] = 1.0
        return v

    # tensor-border data per sector: border curves and ring transversal
    # derivative, from the bi-2 patches at nodes d1/d2 (dual: patch per node).
    cells = np.zeros((n, 4, 4, m))
    for k in range(n):
        d0, d1, d2, d3 = (d(k, j) for j in range(4))
        d1m, d3m = d(k - 1, 1), d(k - 1, 3)   # sector k-1
        d2p, d3p = d(k + 1, 2), d(k + 1, 3)   # sector k+1
        X00 = (d0 + d1 + d2 + d3) / 4
        X03 = (d2 + d3 + d1m + d3m) / 4
        X30 = (d1 + d3 + d2p + d3p) / 4

        # b-edge: border with the ring patch at node d2 (3x3 node grid rows:
        # hole side (d1,d3,d3m) / centre (d0,d2,d1m)); v runs X00 -> X03.
        pos_b = [X00, (d2 + d3) / 2, X03]
        D_b = [2 * ((d0 + d2) / 2 - X00), d2 - d3, 2 * ((d2 + d1m) / 2 - X03)]
        # a-edge: border with the ring patch at node d1; v runs X00 -> X30.
        pos_a = [X00, (d1 + d3) / 2, X30]
        D_a = [2 * ((d0 + d1) / 2 - X00), d1 - d3, 2 * ((d1 + d2p) / 2 - X30)]

        # reparameterization (3^3): a(v) = (1 - c(1-v))/(1-c) (linear, a(0)=1),
        # b(v) = c v(1-v)/(1-c) (quadratic bump).
        a_lin = (1.0, 1.0 / (1 - c))

        def layer(pos, D):
            """Second cap layer along a border: elev(pos) + (1/3)(b*pos' - a*D)."""
            dpos = [2 * (pos[1] - pos[0]), 2 * (pos[2] - pos[1])]
            # b(v) has quadratic BB [0, c/(2(1-c)), 0]; product with the
            # linear-BB derivative dpos gives a cubic with zero end coeffs.
            b_bump = c / (2 * (1 - c))
            term_b = [0 * pos[0],
                      (2 * b_bump * dpos[0]) / 3,
                      (2 * b_bump * dpos[1]) / 3,
                      0 * pos[0]]
            term_a = lin_times_quad(a_lin, D)
            e = elevate_quad(pos)
            return [e[i] + (term_b[i] - term_a[i]) / 3 for i in range(4)]

        lay_b = layer(pos_b, D_b)
        lay_a = layer(pos_a, D_a)
        eb = elevate_quad(pos_b)
        ea = elevate_quad(pos_a)
        for j in range(4):
            cells[k, 0, j] = eb[j]
            cells[k, j, 0] = ea[j]
            cells[k, 1, j] = lay_b[j]
            cells[k, j, 1] = lay_a[j]

    # unknowns per coordinate: centre C, tangent ring t_k (cell (2,3) of
    # patch k == cell (3,2) of patch k-1), interior w_k (cell (2,2)).
    # x = [C, t_0..t_{n-1}, w_0..w_{n-1}]  (1+2n unknowns), each a 4n-stencil.
    nx = 1 + 2 * n
    cell_x = np.zeros((n, 4, 4, nx))
    cell_d = cells.copy()
    for k in range(n):
        xC = np.zeros(nx); xC[0] = 1
        xt = np.zeros(nx); xt[1 + k] = 1               # t_k  at (2,3)
        xt1 = np.zeros(nx); xt1[1 + (k + 1) % n] = 1   # t_{k+1} at (3,2)
        xw = np.zeros(nx); xw[1 + n + k] = 1
        cell_d[k, 3, 3] = 0; cell_x[k, 3, 3] = xC
        cell_d[k, 2, 3] = 0; cell_x[k, 2, 3] = xt
        cell_d[k, 3, 2] = 0; cell_x[k, 3, 2] = xt1
        cell_d[k, 2, 2] = 0; cell_x[k, 2, 2] = xw

    # G1 rows (eq. (1)): edge k between patch k and patch k-1, cells along
    # the shared radial edge from the centre: e0=C, e1=t_k, e2=(1,3), e3=(0,3).
    A_rows, B_rows = [], []
    for k in range(n):
        e = [(cell_x[k, 3, 3], cell_d[k, 3, 3]),
             (cell_x[k, 2, 3], cell_d[k, 2, 3]),
             (cell_x[k, 1, 3], cell_d[k, 1, 3]),
             (cell_x[k, 0, 3], cell_d[k, 0, 3])]
        mh = [(cell_x[k, 3, 2], cell_d[k, 3, 2]),      # patch k side, i=0
              (cell_x[k, 2, 2], cell_d[k, 2, 2]),
              (cell_x[k, 1, 2], cell_d[k, 1, 2]),
              (cell_x[k, 0, 2], cell_d[k, 0, 2])]
        km = (k - 1) % n
        mg = [(cell_x[km, 2, 3], cell_d[km, 2, 3]),    # patch k-1 side
              (cell_x[km, 2, 2], cell_d[km, 2, 2]),
              (cell_x[km, 2, 1], cell_d[km, 2, 1]),
              (cell_x[km, 2, 0], cell_d[km, 2, 0])]
        # rhs cubic coeffs of 2c (1-u) p_u:  [3(e1-e0), 2(e2-e1), (e3-e2), 0]*2c
        for i, wgt in ((0, 3.0), (1, 2.0)):
            lhs_x = 3 * (mh[i][0] + mg[i][0] - 2 * e[i][0])
            lhs_d = 3 * (mh[i][1] + mg[i][1] - 2 * e[i][1])
            rhs_x = 2 * c * wgt * (e[i + 1][0] - e[i][0])
            rhs_d = 2 * c * wgt * (e[i + 1][1] - e[i][1])
            A_rows.append(lhs_x - rhs_x)
            B_rows.append(rhs_d - lhs_d)
    A = np.array([r for r in A_rows])
    B = np.array([r for r in B_rows])

    # central point: n=5 pins C by eq. (5) with e^def = 0.22 (d3 = central
    # n-gon vertex = paper d1; d1,d2 = paper d2,d4 edge-adjacent).
    if n == 5:
        e_def = 0.22
        row_x = np.zeros(nx); row_x[0] = 1
        row_d = np.zeros(m)
        for k in range(n):
            row_d[4 * k + 3] += e_def
            row_d[4 * k + 1] += (1 - n * e_def) / (2 * n)
            row_d[4 * k + 2] += (1 - n * e_def) / (2 * n)
        A = np.vstack([A, row_x[None, :]])
        B = np.vstack([B, row_d[None, :]])

    # F5 objective: sum over patches of the two surviving quintic terms;
    # for bi-3, d^2_a d^3_b p has cubic-diff structure ->
    #   q(c0,c1) = c0^2 + c0*c1 + c1^2 over the two second-differences.
    G_list, H_list = [], []
    for k in range(n):
        for direction in range(2):
            # third difference along one axis, second difference along other:
            # values v[a] = sum_b D3b at row a, then D2 over a -> 2 numbers
            def cx(a, b):
                return cell_x[k, a, b] if direction == 0 else cell_x[k, b, a]

            def cd(a, b):
                return cell_d[k, a, b] if direction == 0 else cell_d[k, b, a]

            for a0 in range(2):
                # c_{a0} = Delta^2_a Delta^3_b at offset a0
                gx = np.zeros(nx); gd = np.zeros(m)
                for da, wa in ((0, 1.0), (1, -2.0), (2, 1.0)):
                    for db, wb in ((0, -1.0), (1, 3.0), (2, -3.0), (3, 1.0)):
                        gx = gx + wa * wb * cx(a0 + da, db)
                        gd = gd + wa * wb * cd(a0 + da, db)
                G_list.append(gx)
                H_list.append(gd)
    G = np.array(G_list)   # (#terms, nx): pairs (c0,c1) consecutive
    H = np.array(H_list)

    # E = sum_pairs [c0,c1]^T M [c0,c1],  M = [[1,.5],[.5,1]]
    npairs = G.shape[0] // 2
    Q = np.zeros((nx, nx)); L = np.zeros((nx, m))
    for p in range(npairs):
        g0, g1 = G[2 * p], G[2 * p + 1]
        h0, h1 = H[2 * p], H[2 * p + 1]
        Q += np.outer(g0, g0) + np.outer(g1, g1) + 0.5 * (np.outer(g0, g1) + np.outer(g1, g0))
        L += np.outer(g0, h0) + np.outer(g1, h1) + 0.5 * (np.outer(g0, h1) + np.outer(g1, h0))

    # KKT: minimize x^T Q x + 2 x^T L d  s.t.  A x = B d   (per input column)
    nc = A.shape[0]
    K = np.zeros((nx + nc, nx + nc))
    K[:nx, :nx] = 2 * Q
    K[:nx, nx:] = A.T
    K[nx:, :nx] = A
    R = np.zeros((nx + nc, m))
    R[:nx] = -2 * L
    R[nx:] = B
    sol, *_ = np.linalg.lstsq(K, R, rcond=None)
    X = sol[:nx]

    # substitute solved unknowns into cells and emit in package row order
    out = np.zeros((16 * n, m))
    for k in range(n):
        for a in range(4):
            for b in range(4):
                out[16 * k + 4 * a + b] = cell_d[k, a, b] + cell_x[k, a, b] @ X
    return out


def verify(name, mine, path):
    pkg = load_package_csv(path)
    err = np.abs(mine - pkg).max()
    print(f"{name}: shape {mine.shape} vs {pkg.shape}  max err {err:.3e}")
    return err




# ---------------------------------------------------------------------------
# 2x2 split bi-3 macro cap (Appendix C, n >= 6)
# ---------------------------------------------------------------------------

def split_cubic(cub):
    """de Casteljau split of a cubic BB at 1/2 -> (left, right)."""
    a, b, c, d = cub
    ab, bc, cd = (a + b) / 2, (b + c) / 2, (c + d) / 2
    abc, bcd = (ab + bc) / 2, (bc + cd) / 2
    mid = (abc + bcd) / 2
    return [a, ab, abc, mid], [mid, bcd, cd, d]


def _gram(deg):
    """Gram matrix of Bernstein basis of given degree on [0,1]."""
    from math import comb
    g = np.zeros((deg + 1, deg + 1))
    for i in range(deg + 1):
        for j in range(deg + 1):
            g[i, j] = comb(deg, i) * comb(deg, j) / comb(2 * deg, i + j) / (2 * deg + 1)
    return g


def ngon_macro(n, sigma_x, sigma_z):
    """INCOMPLETE (PK-317 frontier): 2x2 bi-3 macro cap for n >= 6.

    Verified against the package tables: split borders and transversal
    layers under beta = (3^{2x2}) with gamma*c (4.7e-11); G1 radial weights
    (2c, gamma/2) inner / (gamma/2, 0) outer in per-sub-patch derivatives
    (2.5e-9); centre by eq. (5) with e = (0.02n+1.5)/(n+1.75) (exact); the
    central cells follow a quadratic-jet-through-sigma model whose sigma_ss
    is pinned by assignment (2-bar) to gamma/2 - 1 + sigma_z cos(a/2)/c.
    NOT yet reproduced: the interior-cell selection (the F3 minimization of
    14bi234 Fig. 23c is falsified — the package interiors track a smooth
    guide, see the PK-317 ledger notes); this function's F3 KKT fills the
    interiors to ~2e-2 of the package values.

    sigma_x, sigma_z: the central 2-jet of the characteristic-like map sigma
    (Appendix B/C) after gauge normalization |sigma_s| = 1: sigma_ss =
    (sigma_x, 0), sigma_st = sigma_z * (cos a/2, sin a/2), sigma_tt the
    mirror of sigma_ss across the sector bisector.
    """
    alpha = 2 * np.pi / n
    c = np.cos(alpha)
    gam = (n - 2.3) / (2 * n - 7.72)
    gc = gam * c
    e_def = (0.02 * n + 1.5) / (n + 1.75)
    m = 4 * n

    def d(k, j):
        v = np.zeros(m)
        v[4 * (k % n) + j] = 1.0
        return v

    # ---- pinned border data (values in the 4n input stencil space) ----
    # per sector: split position rows and transversal layers for both edges.
    border = {}   # (k, a, b) -> stencil, for a<=1 or b<=1 cells
    for k in range(n):
        d0, d1, d2, d3 = (d(k, j) for j in range(4))
        d1m, d3m = d(k - 1, 1), d(k - 1, 3)
        d2p, d3p = d(k + 1, 2), d(k + 1, 3)
        X00 = (d0 + d1 + d2 + d3) / 4
        X03 = (d2 + d3 + d1m + d3m) / 4
        X30 = (d1 + d3 + d2p + d3p) / 4
        for edge, pos, D in (
            ('b', [X00, (d2 + d3) / 2, X03],
                  [2 * ((d0 + d2) / 2 - X00), d2 - d3, 2 * ((d2 + d1m) / 2 - X03)]),
            ('a', [X00, (d1 + d3) / 2, X30],
                  [2 * ((d0 + d1) / 2 - X00), d1 - d3, 2 * ((d1 + d2p) / 2 - X30)]),
        ):
            dpos = [2 * (pos[1] - pos[0]), 2 * (pos[2] - pos[1])]
            bump = gc / (2 * (1 - gc))
            term_b = [0 * X00, (2 * bump * dpos[0]) / 3, (2 * bump * dpos[1]) / 3, 0 * X00]
            term_a = lin_times_quad((1.0, 1.0 / (1 - gc)), D)
            Tf = [term_b[i] - term_a[i] for i in range(4)]
            Ppieces = split_cubic(elevate_quad(pos))
            Tpieces = split_cubic(Tf)
            for half in range(2):
                for j in range(4):
                    p0 = Ppieces[half][j]
                    l1 = Ppieces[half][j] + Tpieces[half][j] / 6
                    if edge == 'b':
                        border[(k, 0, 4 * half + j)] = p0
                        border[(k, 1, 4 * half + j)] = l1
                    else:
                        border[(k, 4 * half + j, 0)] = p0
                        border[(k, 4 * half + j, 1)] = l1

    # center by eq. (5)
    C_row = np.zeros(m)
    for k in range(n):
        C_row[4 * k + 3] += e_def
        C_row[4 * k + 1] += (1 - n * e_def) / (2 * n)
        C_row[4 * k + 2] += (1 - n * e_def) / (2 * n)

    # ---- unknown vector: per sector 9 cells + theta (5 quadratic-jet DOF) ----
    # per-sector unknown cells: (2,2),(2,5),(2,6),(5,2),(5,5),(5,6),(6,2),(6,5)
    # (interior; (6,6) is jet-determined) + the radial edge cell (2,7).
    FREE = [(2, 2), (2, 5), (2, 6), (5, 2), (5, 5), (5, 6), (6, 2), (6, 5), (2, 7)]
    nfree = len(FREE)
    N = nfree * n + 5
    TH = nfree * n   # theta offset

    # sigma jet, per-sector composition:
    #   f_s^k = grad . rho_k s_s     f_ss^k = s_s^T rho_k^T H rho_k s_s + grad . rho_k s_ss
    s_s = np.array([1.0, 0.0])
    s_t = np.array([np.cos(alpha), np.sin(alpha)])
    s_ss = np.array([sigma_x, 0.0])
    s_st = sigma_z * np.array([np.cos(alpha / 2), np.sin(alpha / 2)])
    s_tt = np.array([sigma_x * np.cos(alpha), sigma_x * np.sin(alpha)])

    def jet_rows(k):
        """Return (fs, ft, fss, fst, ftt) as rows over theta=(g1,g2,h11,h12,h22)."""
        ca, sa = np.cos(k * alpha), np.sin(k * alpha)
        rho = np.array([[ca, -sa], [sa, ca]])
        vs, vt = rho @ s_s, rho @ s_t
        ws, wt, wst = rho @ s_ss, rho @ s_tt, rho @ s_st

        def hq(u, v):
            return np.array([u[0] * v[0], u[0] * v[1] + u[1] * v[0], u[1] * v[1]])

        fs = np.concatenate([vs, np.zeros(3)])
        ft = np.concatenate([vt, np.zeros(3)])
        fss = np.concatenate([ws, hq(vs, vs)])
        fst = np.concatenate([wst, hq(vs, vt)])
        ftt = np.concatenate([wt, hq(vt, vt)])
        return fs, ft, fss, fst, ftt

    jets = [jet_rows(k) for k in range(n)]

    # ---- cell map: (k,a,b) -> (x_row (N,), d_row (m,)) with a,b in 0..7 ----
    from functools import lru_cache

    def cell(k, a, b):
        k = k % n
        # deduplicate junction indices
        if a == 4:
            a = 3
        if b == 4:
            b = 3
        if a == 3 and b == 3:
            xa, da = cell(k, 2, 3)
            xb, db = cell(k, 5, 3)
            return (xa + xb) / 2, (da + db) / 2
        if a == 3:
            xa, da = cell(k, 2, b)
            xb, db = cell(k, 5, b)
            return (xa + xb) / 2, (da + db) / 2
        if b == 3:
            xa, da = cell(k, a, 2)
            xb, db = cell(k, a, 5)
            return (xa + xb) / 2, (da + db) / 2
        # radial sharing: a=7 row belongs to patch k+1's b=7 column
        if a == 7 and b != 7:
            return cell(k + 1, b, 7)
        x = np.zeros(N)
        zero = np.zeros(m)
        if a <= 1 or b <= 1:
            if (k, a, b) in border:
                return x, border[(k, a, b)]
        if (a, b) == (7, 7):
            return x, C_row.copy()
        fs, ft, fss, fst, ftt = jets[k]
        if (a, b) == (6, 7):
            x[TH:] = -fs / 3
            return x, C_row.copy()
        if (a, b) == (7, 6):
            x[TH:] = -ft / 3
            return x, C_row.copy()
        if (a, b) == (5, 7):
            x[TH:] = -2 * fs / 3 + fss / 6
            return x, C_row.copy()
        if (a, b) == (7, 5):
            x[TH:] = -2 * ft / 3 + ftt / 6
            return x, C_row.copy()
        if (a, b) == (6, 6):
            x[TH:] = -fs / 3 - ft / 3 + fst / 9
            return x, C_row.copy()
        if (a, b) in FREE:
            x[nfree * k + FREE.index((a, b))] = 1.0
            return x, zero
        raise AssertionError(f"unmapped cell {(k, a, b)}")

    # ---- G1 constraints along each radial edge ----
    A_rows, B_rows = [], []
    for k in range(n):
        # edge column b=7 of patch k, from centre outward
        E = [cell(k, 7 - i, 7) for i in range(8)]  # index by distance from centre
        Nk = [cell(k, 7 - i, 6) for i in range(8)]
        Mk = [cell(k - 1, 6, 7 - i) for i in range(8)]
        # pieces: inner (dist 0..3, junction at dist 3.5 -> BB [E0,E1,E2,J]),
        # outer (BB [J,E5?,..]) -- use dedup: distances 0,1,2 -> cells
        # (7,7),(6,7),(5,7); junction J = ((5,7)+(2,7))/2; outer: J,(2,7),(1,7),(0,7)
        def bb_col(cells_by_dist):
            J = tuple((cells_by_dist[2][t] + cells_by_dist[5][t]) / 2 for t in range(2))
            inner = [cells_by_dist[0], cells_by_dist[1], cells_by_dist[2], J]
            outer = [J, cells_by_dist[5], cells_by_dist[6], cells_by_dist[7]]
            return inner, outer
        Ei, Eo = bb_col(E)
        Ni, No = bb_col(Nk)
        Mi, Mo = bb_col(Mk)
        # (1-bar)/(1) weights in per-sub-patch local derivatives: the printed
        # gamma*u / gamma*(1-u) terms use the global edge parameter, so the
        # junction value is gamma/2 (decoded and pinned against the package).
        for piece, (e, nn, mm, w) in enumerate((
            (Ei, Ni, Mi, (2 * c, gam / 2)),
            (Eo, No, Mo, (gam / 2, 0.0)),
        )):
            eu = [3 * (e[1][t] - e[0][t]) for t in range(2)], \
                 [3 * (e[2][t] - e[1][t]) for t in range(2)], \
                 [3 * (e[3][t] - e[2][t]) for t in range(2)]
            rhs = [None] * 4
            # cubic BB of w(u)*e_u: lin x quad product per component
            for t in range(2):
                q = [eu[0][t], eu[1][t], eu[2][t]]
                p = lin_times_quad(w, q)
                for i in range(4):
                    if rhs[i] is None:
                        rhs[i] = [None, None]
                    rhs[i][t] = p[i]
            for i in range(4):
                lhs_x = 3 * (nn[i][0] + mm[i][0] - 2 * e[i][0])
                lhs_d = 3 * (nn[i][1] + mm[i][1] - 2 * e[i][1])
                A_rows.append(lhs_x - rhs[i][0])
                B_rows.append(rhs[i][1] - lhs_d)
    A = np.array(A_rows)
    B = np.array(B_rows)

    # ---- F3 over all 4n sub-patches ----
    G1d, G2d, G3d = _gram(1), _gram(2), _gram(3)
    Q = np.zeros((N, N))
    L = np.zeros((N, m))

    def accumulate(coeff_x, coeff_d, gram_pair, weight):
        # coeff arrays: list-of-lists BB grid of (x_row, d_row)
        Gi, Gj = gram_pair
        ni, nj = len(coeff_x), len(coeff_x[0])
        for i1 in range(ni):
            for j1 in range(nj):
                for i2 in range(ni):
                    for j2 in range(nj):
                        wgt = weight * Gi[i1, i2] * Gj[j1, j2]
                        if wgt == 0.0:
                            continue
                        Q_ = np.outer(coeff_x[i1][j1], coeff_x[i2][j2])
                        Q[:, :] += wgt * (Q_ + Q_.T) / 2 * 2
                        L[:, :] += wgt * (np.outer(coeff_x[i1][j1], coeff_d[i2][j2])
                                          + np.outer(coeff_x[i2][j2], coeff_d[i1][j1]))

    for k in range(n):
        for I in range(2):
            for J in range(2):
                cx = [[None] * 4 for _ in range(4)]
                cd = [[None] * 4 for _ in range(4)]
                for i in range(4):
                    for j in range(4):
                        cx[i][j], cd[i][j] = cell(k, 4 * I + i, 4 * J + j)
                # third differences
                def diff3(grid, axis):
                    out = []
                    for i in range(4):
                        row = []
                        for j in range(1):
                            pass
                        out.append(row)
                    return out
                # (3,0): 6 * D3 along i -> row vector over j (cubic in j)
                d3i_x = [sum(w * cx[i0 + di][j] for di, w in ((0, -1), (1, 3), (2, -3), (3, 1))) for i0 in range(1) for j in range(4)]
                d3i_d = [sum(w * cd[i0 + di][j] for di, w in ((0, -1), (1, 3), (2, -3), (3, 1))) for i0 in range(1) for j in range(4)]
                accumulate([d3i_x], [d3i_d], (np.array([[1.0]]), G3d), 36.0)
                d3j_x = [[sum(w * cx[i][j0 + dj] for dj, w in ((0, -1), (1, 3), (2, -3), (3, 1)))] for i in range(4) for j0 in range(1)]
                d3j_d = [[sum(w * cd[i][j0 + dj] for dj, w in ((0, -1), (1, 3), (2, -3), (3, 1)))] for i in range(4) for j0 in range(1)]
                accumulate(d3j_x, d3j_d, (G3d, np.array([[1.0]])), 36.0)
                # (2,1): 6 D2_i x 3 D1_j -> deg (1,2) grid
                c21_x = [[18 * sum(wi * wj * cx[i0 + di][j0 + dj]
                                   for di, wi in ((0, 1), (1, -2), (2, 1))
                                   for dj, wj in ((0, -1), (1, 1)))
                          for j0 in range(3)] for i0 in range(2)]
                c21_d = [[18 * sum(wi * wj * cd[i0 + di][j0 + dj]
                                   for di, wi in ((0, 1), (1, -2), (2, 1))
                                   for dj, wj in ((0, -1), (1, 1)))
                          for j0 in range(3)] for i0 in range(2)]
                accumulate(c21_x, c21_d, (G1d, G2d), 3.0)
                c12_x = [[18 * sum(wi * wj * cx[i0 + di][j0 + dj]
                                   for di, wi in ((0, -1), (1, 1))
                                   for dj, wj in ((0, 1), (1, -2), (2, 1)))
                          for j0 in range(2)] for i0 in range(3)]
                c12_d = [[18 * sum(wi * wj * cd[i0 + di][j0 + dj]
                                   for di, wi in ((0, -1), (1, 1))
                                   for dj, wj in ((0, 1), (1, -2), (2, 1)))
                          for j0 in range(2)] for i0 in range(3)]
                accumulate(c12_x, c12_d, (G2d, G1d), 3.0)

    # ---- KKT solve ----
    nc = A.shape[0]
    K = np.zeros((N + nc, N + nc))
    K[:N, :N] = 2 * Q
    K[:N, N:] = A.T
    K[N:, :N] = A
    R = np.zeros((N + nc, m))
    R[:N] = -2 * L
    R[N:] = B
    sol, *_ = np.linalg.lstsq(K, R, rcond=None)
    X = sol[:N]

    out = np.zeros((64 * n, m))
    for k in range(n):
        for I in range(2):
            for J in range(2):
                for i in range(4):
                    for j in range(4):
                        xr, dr = cell(k, 4 * I + i, 4 * J + j)
                        out[64 * k + 16 * (I + 2 * J) + 4 * i + j] = dr + xr @ X
    return out


def _macro_err(n, x, z, pkg):
    return np.abs(ngon_macro(n, x, z) - pkg).max()


def fit_sigma_jet(n):
    """Fit the two sigma-jet scalars against the package table (diagnostic)."""
    pkg = load_package_csv(BASE + f'ngonSct{n}.csv')
    # Appendix B candidate as start
    d1 = (0.24 * n - 0.75) / (n - 3.16)
    d2 = (0.44 * n - 0.75) / (n - 4 + 2.02 * np.sqrt(2))
    al = 2 * np.pi / n
    c = np.cos(al)
    s20 = (1 - 4 / (3 * c)) * d1 + (4 * np.cos(al / 2) / (3 * c)) * d2
    x0 = 3 * (s20 - 2 * d1) / d1
    z0 = 4 * (d2 - 2 * d1 * np.cos(al / 2)) / d1
    print(f"n={n}: appendix-B candidate x={x0:.6f} z={z0:.6f} -> err {_macro_err(n, x0, z0, pkg):.3e}")
    print(f"      halved x/2,z/2 -> err {_macro_err(n, x0/2, z0/2, pkg):.3e}")
    # simple coordinate-descent fit (no scipy in this venv)
    x, z = x0, z0
    best = _macro_err(n, x, z, pkg)
    step = 0.1
    while step > 1e-11:
        improved = False
        for dx, dz in ((step, 0), (-step, 0), (0, step), (0, -step),
                       (step, step), (step, -step), (-step, step), (-step, -step)):
            e = _macro_err(n, x + dx, z + dz, pkg)
            if e < best:
                best, x, z = e, x + dx, z + dz
                improved = True
        if not improved:
            step /= 2
    print(f"      fitted x={x:.9f} z={z:.9f} -> err {best:.3e}")
    return x, z, best


# ---------------------------------------------------------------------------
# eop (star/CC-net) tables via virtual refinement composition
# ---------------------------------------------------------------------------

def refinement_matrix(n):
    """Bi-2 virtual refinement (paper Fig. 5/6a): DS points from CC points.

    CC columns: (q1^0, q2^0, q1^1, q2^1, ..., E) — 2n+1. DS rows in the ngon
    table's own sector order: ngon sector s pairs with the CC quad
    (E, q1^{-s}, q2^{-s}, q1^{-s+1}) (decoded; sector orientations of the two
    package layouts run opposite).
    """
    R = np.zeros((4 * n, 2 * n + 1))
    E = 2 * n
    for s in range(n):
        qa = 2 * ((-s) % n)        # q1^{-s}
        qd = 2 * ((-s) % n) + 1    # q2^{-s}
        qb = 2 * ((-s + 1) % n)    # q1^{-s+1}
        for m_, (wa, wd, wb, we) in enumerate(((3, 9, 3, 1), (9, 3, 1, 3),
                                               (1, 3, 9, 3), (3, 1, 3, 9))):
            R[4 * s + m_, qa] = wa / 16
            R[4 * s + m_, qd] = wd / 16
            R[4 * s + m_, qb] = wb / 16
            R[4 * s + m_, E] = we / 16
    return R


def eop_small(n):
    """eop table for n in {3,5}: ngon table composed with virtual refinement.

    Package eop row (s,a,b) holds the composed ngon row (s, 3-b, 3-a) — the
    two layouts transpose the patch grid and count the radial layers from the
    opposite end (decoded against the package tables; same convention for
    n = 3 and n = 5).
    """
    ngon = ngon_small(n)
    comp = ngon @ refinement_matrix(n)
    out = np.zeros_like(comp)
    for s in range(n):
        for a in range(4):
            for b in range(4):
                out[16 * s + 4 * a + b] = comp[16 * s + 4 * (3 - b) + (3 - a)]
    return out


if __name__ == '__main__':
    worst = 0.0
    for n in (3, 5):
        worst = max(worst, verify(f'ngonSct{n}', ngon_small(n), BASE + f'ngonSct{n}.csv'))
        worst = max(worst, verify(f'eopSct{n}', eop_small(n), BASE + f'eopSct{n}.csv'))
    print('worst:', worst)
