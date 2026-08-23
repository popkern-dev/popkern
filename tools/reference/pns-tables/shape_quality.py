"""Shape-quality battery for the PnS macro cap tables (PK-318).

Measures the surface fairness of the 2x2 bi-3 macro star/ngon caps
(n in {6,7,8}) produced by a mask-table set, on a canonical battery of
curved input nets (dome / saddle / wavy), and compares table sets:

  - mean-curvature jump across the internal C1 junctions of each sector
    (the quantity guided interiors improve; sector G1 seams are reported
    separately),
  - total variation of mean curvature over the cap (fairness proxy),
  - worst normal-angle difference between two table sets' surfaces,
  - extension-agreement fingerprint: how far neighbouring sub-patches are
    from being polynomial extensions of each other across the internal
    junctions (coefficient-space, per unit net scale).

Run from the repo root (compares committed tables against the reference
package oracle when it is present, else prints committed-only metrics):

    ~/venvs/ocp/bin/python tools/reference/pns-tables/shape_quality.py
    ~/venvs/ocp/bin/python tools/reference/pns-tables/shape_quality.py --render OUTDIR

The committed acceptance gate (asserted by emit_all.py after regeneration,
and by this script's exit code when the oracle is present): on every
(family, n, net) battery cell the committed tables' internal-junction
curvature jump stays within GATE_JUMP_FACTOR and curvature total-variation
within GATE_TV_FACTOR of the reference package's, and the worst
normal-angle difference stays below NORMAL_GATE_DEG.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
REF = os.path.join(REPO, 'docs/code/ACM-ToMS-PnS-1032/1032/Src/Source/src/Patch/Table')
OURS = os.path.join(REPO, 'crates/popkern-core/data/pns-masks')

# Gate levels: worst measured for the committed guided interiors are
# jump 1.66x / TV 1.09x of the reference and 0.64 deg normal deviation
# (PK-318 battery); an F3-class interior regression sits at 5-10x jump and
# ~1.5 deg, so these bounds have both honest headroom and detection power.
GATE_JUMP_FACTOR = 1.8
GATE_TV_FACTOR = 1.2
NORMAL_GATE_DEG = 0.8


def load_csv(path):
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            rows.append([float(x) for x in line.split(',')])
    return np.array(rows)


# ---------------------------------------------------------------------------
# canonical input nets (see eop_ngon.py for the stencil layouts)
# ---------------------------------------------------------------------------

def shape_z(x, y, kind):
    if kind == 'dome':
        return 0.35 * (x * x + y * y)
    if kind == 'saddle':
        return 0.35 * (x * x - y * y)
    if kind == 'wavy':
        r2 = x * x + y * y
        th = np.arctan2(y, x)
        return 0.3 * r2 + 0.08 * np.cos(2 * th) * r2
    raise ValueError(kind)


def cc_net(n, kind):
    """CC-net (eop input): columns q1^0,q2^0,...,q1^{n-1},q2^{n-1},E."""
    al = 2 * np.pi / n
    pts = []
    for s in range(n):
        a1, a2 = s * al, s * al + al / 2
        q1 = np.array([np.cos(a1), np.sin(a1)])
        q2 = 1.4 * np.array([np.cos(a2), np.sin(a2)])
        pts.append([q1[0], q1[1], shape_z(q1[0], q1[1], kind)])
        pts.append([q2[0], q2[1], shape_z(q2[0], q2[1], kind)])
    pts.append([0.0, 0.0, shape_z(0.0, 0.0, kind)])
    return np.array(pts)


def ds_net(n, kind):
    """DS-net (ngon input): sector groups d0 (diagonal), d1, d2, d3."""
    al = 2 * np.pi / n
    de = al / 3
    pts = []
    for k in range(n):
        th = k * al
        for p in (2.4 * np.array([np.cos(th), np.sin(th)]),
                  1.9 * np.array([np.cos(th + de), np.sin(th + de)]),
                  1.9 * np.array([np.cos(th - de), np.sin(th - de)]),
                  1.0 * np.array([np.cos(th), np.sin(th)])):
            pts.append([p[0], p[1], shape_z(p[0], p[1], kind)])
    return np.array(pts)


# ---------------------------------------------------------------------------
# bi-3 evaluation
# ---------------------------------------------------------------------------

def bern3(t):
    t = np.asarray(t)
    mt = 1 - t
    return np.stack([mt ** 3, 3 * t * mt ** 2, 3 * t ** 2 * mt, t ** 3], -1)


def dbern3(t):
    t = np.asarray(t)
    mt = 1 - t
    return np.stack([-3 * mt ** 2, 3 * mt ** 2 - 6 * t * mt,
                     6 * t * mt - 3 * t ** 2, 3 * t ** 2], -1)


def ddbern3(t):
    t = np.asarray(t)
    return np.stack([6 * (1 - t), -12 + 18 * t, 6 - 18 * t, 6 * t], -1)


def patch_grids(coeffs, n):
    """(64n,3) macro rows -> dict (k,I,J) -> 4x4x3 grid (row 64k+16(I+2J)+4i+j)."""
    P = {}
    for k in range(n):
        for I in range(2):
            for J in range(2):
                g = np.zeros((4, 4, 3))
                for i in range(4):
                    for j in range(4):
                        g[i, j] = coeffs[64 * k + 16 * (I + 2 * J) + 4 * i + j]
                P[(k, I, J)] = g
    return P


def normals_curvature(g, us, vs):
    Bu, Bv, Du, Dv = bern3(us), bern3(vs), dbern3(us), dbern3(vs)
    DDu, DDv = ddbern3(us), ddbern3(vs)
    pos = np.einsum('ua,vb,abc->uvc', Bu, Bv, g)
    Su = np.einsum('ua,vb,abc->uvc', Du, Bv, g)
    Sv = np.einsum('ua,vb,abc->uvc', Bu, Dv, g)
    Suu = np.einsum('ua,vb,abc->uvc', DDu, Bv, g)
    Suv = np.einsum('ua,vb,abc->uvc', Du, Dv, g)
    Svv = np.einsum('ua,vb,abc->uvc', Bu, DDv, g)
    Nr = np.cross(Su, Sv)
    nn = np.linalg.norm(Nr, axis=-1, keepdims=True)
    nn = np.where(nn == 0, 1, nn)
    N = Nr / nn
    E = np.einsum('uvc,uvc->uv', Su, Su)
    F = np.einsum('uvc,uvc->uv', Su, Sv)
    G = np.einsum('uvc,uvc->uv', Sv, Sv)
    L = np.einsum('uvc,uvc->uv', Suu, N)
    M = np.einsum('uvc,uvc->uv', Suv, N)
    Nn = np.einsum('uvc,uvc->uv', Svv, N)
    den = E * G - F * F
    den = np.where(np.abs(den) < 1e-300, 1e-300, den)
    H = (E * Nn - 2 * F * M + G * L) / (2 * den)
    return pos, N, H


def cap_field(P, res=33):
    us = np.linspace(0, 1, res)
    return {key: normals_curvature(g, us, us) for key, g in P.items()}


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------

def normal_angle_diff(fa, fb):
    worst = 0.0
    for key in fa:
        dot = np.clip(np.einsum('uvc,uvc->uv', fa[key][1], fb[key][1]), -1, 1)
        worst = max(worst, np.arccos(dot).max())
    return worst


def junction_H_jumps(P, n, res=65):
    """(max internal C1-junction jump, max sector G1-seam jump) of H."""
    us = np.linspace(0, 1, res)
    max_int, max_sec = 0.0, 0.0

    def H_line(g, axis, val):
        if axis == 'u':
            return normals_curvature(g, np.array([val]), us)[2][0]
        return normals_curvature(g, us, np.array([val]))[2][:, 0]

    for k in range(n):
        for J in range(2):
            max_int = max(max_int, np.abs(
                H_line(P[(k, 0, J)], 'u', 1.0) - H_line(P[(k, 1, J)], 'u', 0.0)).max())
        for I in range(2):
            max_int = max(max_int, np.abs(
                H_line(P[(k, I, 0)], 'v', 1.0) - H_line(P[(k, I, 1)], 'v', 0.0)).max())
        for half in range(2):
            # radial sharing: patch k cells (a,7) == patch k+1 cells (7,a)
            max_sec = max(max_sec, np.abs(
                H_line(P[(k, half, 1)], 'v', 1.0)
                - H_line(P[((k + 1) % n, 1, half)], 'u', 1.0)).max())
    return max_int, max_sec


def curvature_tv(field):
    tv = 0.0
    for key in field:
        H = field[key][2]
        tv += np.abs(np.diff(H, axis=0)).sum() + np.abs(np.diff(H, axis=1)).sum()
    return tv


def _extend_cubic(p):
    """BB coeffs of the cubic's extension from [0,1] onto [1,2]."""
    return np.array([p[3],
                     2 * p[3] - p[2],
                     4 * p[3] - 4 * p[2] + p[1],
                     8 * p[3] - 12 * p[2] + 6 * p[1] - p[0]])


def extension_agreement(P, n, scale):
    """Worst coefficient distance between a sub-patch and its neighbour's
    polynomial extension across the internal junctions, per unit net scale."""
    worst = 0.0
    for k in range(n):
        for J in range(2):
            gL, gR = P[(k, 0, J)], P[(k, 1, J)]
            for j in range(4):
                worst = max(worst, np.abs(_extend_cubic(gL[:, j]) - gR[:, j]).max())
        for I in range(2):
            gB, gT = P[(k, I, 0)], P[(k, I, 1)]
            for i in range(4):
                worst = max(worst, np.abs(_extend_cubic(gB[i, :]) - gT[i, :]).max())
    return worst / scale


# ---------------------------------------------------------------------------
# battery
# ---------------------------------------------------------------------------

BATTERY = [(fam, n, kind)
           for fam in ('eop', 'ngon')
           for n in (6, 7, 8)
           for kind in ('dome', 'saddle', 'wavy')]


def measure(table, fam, n, kind):
    net = cc_net(n, kind) if fam == 'eop' else ds_net(n, kind)
    P = patch_grids(table @ net, n)
    f = cap_field(P)
    ji, js = junction_H_jumps(P, n)
    scale = np.abs(net).max()
    return {'P': P, 'field': f, 'jump_int': ji, 'jump_sec': js,
            'tv': curvature_tv(f), 'ext': extension_agreement(P, n, scale)}


def render_pair(Pa, Pb, title, path, res=97, bands=14):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    us = np.linspace(0, 1, res)
    Ld = np.array([0.3, 0.25, 0.92])
    Ld = Ld / np.linalg.norm(Ld)
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.5))
    for ax, P, name in ((axes[0], Pa, 'committed'), (axes[1], Pb, 'reference')):
        for key in P:
            pos, N, _ = normals_curvature(P[key], us, us)
            s = np.einsum('uvc,c->uv', N, Ld)
            ax.pcolormesh(pos[:, :, 0], pos[:, :, 1],
                          (np.sin(np.pi * bands * s) > 0).astype(float),
                          cmap='gray', vmin=0, vmax=1, shading='gouraud')
        ax.set_aspect('equal')
        ax.set_title(name)
        ax.axis('off')
    ax = axes[2]
    vmax, polys = 0.0, []
    for key in Pa:
        pos, Na, _ = normals_curvature(Pa[key], us, us)
        _, Nb, _ = normals_curvature(Pb[key], us, us)
        ang = np.degrees(np.arccos(np.clip(
            np.einsum('uvc,uvc->uv', Na, Nb), -1, 1)))
        vmax = max(vmax, ang.max())
        polys.append((pos, ang))
    for pos, ang in polys:
        pc = ax.pcolormesh(pos[:, :, 0], pos[:, :, 1], ang, cmap='inferno',
                           vmin=0, vmax=max(vmax, 1e-9), shading='gouraud')
    fig.colorbar(pc, ax=ax, shrink=0.8, label='normal diff (deg)')
    ax.set_aspect('equal')
    ax.set_title(f'normal deviation (max {vmax:.3f} deg)')
    ax.axis('off')
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def run_battery(render_dir=None):
    """Compare committed tables against the reference oracle (if present).
    Returns (rows, failures); failures non-empty means a gate broke."""
    have_ref = os.path.isdir(REF)
    rows, failures = [], []
    for fam, n, kind in BATTERY:
        ours = load_csv(os.path.join(OURS, f'{fam}Sct{n}.csv'))
        a = measure(ours, fam, n, kind)
        row = {'fam': fam, 'n': n, 'kind': kind, 'ours': a}
        if have_ref:
            ref = load_csv(os.path.join(REF, f'{fam}Sct{n}.csv'))
            b = measure(ref, fam, n, kind)
            row['ref'] = b
            row['nang'] = np.degrees(normal_angle_diff(a['field'], b['field']))
            if a['jump_int'] > GATE_JUMP_FACTOR * max(b['jump_int'], 1e-6):
                failures.append(f"{fam}{n} {kind}: jump_int {a['jump_int']:.4f} "
                                f"> {GATE_JUMP_FACTOR} x ref {b['jump_int']:.4f}")
            if a['tv'] > GATE_TV_FACTOR * b['tv']:
                failures.append(f"{fam}{n} {kind}: H-TV {a['tv']:.3f} "
                                f"> {GATE_TV_FACTOR} x ref {b['tv']:.3f}")
            if row['nang'] > NORMAL_GATE_DEG:
                failures.append(f"{fam}{n} {kind}: normal dev {row['nang']:.3f} deg "
                                f"> {NORMAL_GATE_DEG}")
            if render_dir:
                os.makedirs(render_dir, exist_ok=True)
                render_pair(a['P'], b['P'], f'{fam} n={n} {kind}',
                            os.path.join(render_dir, f'{fam}{n}_{kind}.png'))
        rows.append(row)
    return rows, failures


def main(argv):
    render_dir = None
    if '--render' in argv:
        render_dir = argv[argv.index('--render') + 1]
    rows, failures = run_battery(render_dir)
    have_ref = 'ref' in rows[0]
    for r in rows:
        a = r['ours']
        line = (f"{r['fam']}{r['n']} {r['kind']:7s} "
                f"Hjump-int {a['jump_int']:7.4f} sec {a['jump_sec']:7.4f} "
                f"H-TV {a['tv']:8.3f} ext {a['ext']:.4f}")
        if have_ref:
            b = r['ref']
            line += (f" | ref Hjump-int {b['jump_int']:7.4f} H-TV {b['tv']:8.3f} "
                     f"ext {b['ext']:.4f} | Ndev {r['nang']:6.3f} deg")
        print(line)
    if not have_ref:
        print('reference package not present: committed-only metrics, no gate')
        return 0
    if failures:
        print('\nGATE FAILURES:')
        for f in failures:
            print(' ', f)
        return 1
    print(f'\nAll shape gates pass (jump within {GATE_JUMP_FACTOR}x, TV within '
          f'{GATE_TV_FACTOR}x of reference, normal dev < {NORMAL_GATE_DEG} deg).')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
