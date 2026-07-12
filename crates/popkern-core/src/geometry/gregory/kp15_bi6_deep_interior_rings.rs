//! KP15 bi-6 deep-interior ring solvers and (S3) rules.
//!
//! The interior-ring stages of the deep-interior population (see
//! [`super::kp15_bi6_deep_interior`] for the orchestration and the
//! full construction narrative): the eq (7) `2n×2n` ring, the
//! eq (8) `n×n` ring, and the (S3) `µ̃`/`µ̄` rules. Split out of the
//! orchestrator module to keep each file within the 500-line
//! budget. All indexing is in **script** convention.

use super::super::direct_linalg::solve_dense;
use super::super::kp15_bi6_mu::{mu_bar, mu_ddot, mu_dot, mu_tdot, mu_tilde};
use super::super::kp15_bi6_reparam::Bi6Reparams;
use super::Bi6DeepInteriorFree;
use crate::geometry::Point3;

pub(super) type Grid = [[Point3; 7]; 7];
const ZERO: Point3 = Point3 {
    x: 0.0,
    y: 0.0,
    z: 0.0,
};

/// `µ̇` with the reconstructed sum-entry `[0] = 2 − Σµ̇_{1..5}`
/// prepended (oracle `mu1_tab`). Result indexed `[0..=5]`.
fn mu1_full(p: &Bi6Reparams) -> [f64; 6] {
    let m = mu_dot(p);
    let s: f64 = m.iter().sum();
    [2.0 - s, m[0], m[1], m[2], m[3], m[4]]
}

/// `µ̈` with `[0] = −Σµ̈_{1..8}` prepended (oracle `mu2_tab`).
fn mu2_full(p: &Bi6Reparams) -> [f64; 9] {
    let m = mu_ddot(p);
    let s: f64 = m.iter().sum();
    [-s, m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7]]
}

/// `⃛µ` with `[0] = −Σ⃛µ_{1..8}` prepended (oracle `mu3_tab`; the
/// sum runs over indices 1..8 only, **not** 1..11).
fn mu3_full(p: &Bi6Reparams) -> [f64; 12] {
    let m = mu_tdot(p);
    let s: f64 = m.iter().take(8).sum();
    [
        -s, m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9], m[10],
    ]
}

/// `µ̃` with `[0] = 1 − Σµ̃_{1..6}` prepended (oracle `mut_tab`).
fn mut_full(p: &Bi6Reparams) -> [f64; 11] {
    let m = mu_tilde(p);
    let s: f64 = m.iter().take(6).sum();
    [
        1.0 - s,
        m[0],
        m[1],
        m[2],
        m[3],
        m[4],
        m[5],
        m[6],
        m[7],
        m[8],
        m[9],
    ]
}

/// `µ̄` with `[0] = −Σµ̄_{1..6}` prepended (oracle `mub_tab`).
fn mub_full(p: &Bi6Reparams) -> [f64; 11] {
    let m = mu_bar(p);
    let s: f64 = m.iter().take(6).sum();
    [
        -s, m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7], m[8], m[9],
    ]
}

/// Least-squares solution of `A x = b` (consistent, full column
/// rank) via the normal equations `AᵀA x = Aᵀb`. `a` is
/// `rows × cols` with `rows ≥ cols`. Returns `None` if `AᵀA` is
/// singular. Used for the eq (7) ring, which is `2n × (2n − 1)`
/// (over-determined but consistent) at `n ∈ {3, 6}` and square
/// otherwise.
fn solve_normal_equations(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let rows = a.len();
    let cols = a.first().map_or(0, Vec::len);
    if rows == 0 || cols == 0 || rows < cols {
        return None;
    }
    let mut ata = vec![vec![0.0; cols]; cols];
    let mut atb = vec![0.0; cols];
    for i in 0..cols {
        for j in 0..cols {
            let mut s = 0.0;
            for row in a.iter() {
                s += row[i] * row[j];
            }
            ata[i][j] = s;
        }
        let mut s = 0.0;
        for (row, &bv) in a.iter().zip(b.iter()) {
            s += row[i] * bv;
        }
        atb[i] = s;
    }
    solve_dense(ata, atb)
}

/// Read the `axis`-th component (0 = x, 1 = y, 2 = z) of `p`.
#[inline]
fn comp(p: Point3, axis: usize) -> f64 {
    match axis {
        0 => p.x,
        1 => p.y,
        _ => p.z,
    }
}

/// Build a `Point3` from a per-axis closure.
#[inline]
fn point_from_axes(f: impl Fn(usize) -> f64) -> Point3 {
    Point3::new(f(0), f(1), f(2))
}

/// Solve the (S2) eq (7) ring for `p[r][5,4]` and `p[r][4,5]`.
pub(super) fn solve_eq7_ring(
    g: &mut [Grid],
    params: &Bi6Reparams,
    free: &Bi6DeepInteriorFree,
    n: usize,
) -> Result<(), String> {
    let mu1 = mu1_full(params);
    let mu2 = mu2_full(params);
    let sf = |k: usize| k % n;
    let p0_45_free = n == 3 || n == 6;
    let r_start = if p0_45_free { 1 } else { 0 };
    // Unknown layout: x54[r] -> col r (r = 0..n); x45[r] -> col
    // n + (r - r_start) for r = r_start..n.
    let col_x54 = |r: usize| r;
    let col_x45 = |r: usize| n + (r - r_start);
    let ncols = n + (n - r_start);
    let nrows = 2 * n;

    let mut a = vec![vec![0.0; ncols]; nrows];
    let coord_of = [|p: Point3| p.x, |p: Point3| p.y, |p: Point3| p.z];
    let mut rhs = [vec![0.0; nrows], vec![0.0; nrows], vec![0.0; nrows]];

    for r in 0..n {
        let t = sf(r + 1);
        // Row 2r+1: x45[t] + x54[r] = K1r.
        let row_hi = 2 * r + 1;
        a[row_hi][col_x54(r)] += 1.0;
        let t_pinned = p0_45_free && t == 0;
        if !t_pinned {
            a[row_hi][col_x45(t)] += 1.0;
        }
        // Row 2r: x54[t] − x45[r] − µ̈_5·x54[r] = K2r.
        let row_lo = 2 * r;
        a[row_lo][col_x54(t)] += 1.0;
        a[row_lo][col_x54(r)] += -mu2[5];
        let r_pinned = p0_45_free && r == 0;
        if !r_pinned {
            a[row_lo][col_x45(r)] += -1.0;
        }

        for (axis, coord) in coord_of.iter().enumerate() {
            let k1 = mu1[0] * coord(g[r][6][6])
                + mu1[1] * coord(g[r][6][5])
                + mu1[2] * coord(g[r][6][4])
                + mu1[3] * coord(g[r][6][2])
                + mu1[4] * coord(g[r][6][1])
                + mu1[5] * coord(g[r][6][0]);
            let mut b_hi = k1;
            if t_pinned {
                b_hi -= coord(free.p210);
            }
            rhs[axis][row_hi] = b_hi;

            let k2 = mu2[0] * coord(g[r][6][6])
                + mu2[1] * coord(g[r][5][6])
                + mu2[2] * coord(g[r][6][5])
                + mu2[3] * coord(g[r][5][5])
                + mu2[4] * coord(g[r][6][4])
                + mu2[6] * coord(g[r][6][2])
                + mu2[7] * coord(g[r][6][1])
                + mu2[8] * coord(g[r][6][0]);
            let mut b_lo = k2;
            if r_pinned {
                b_lo += coord(free.p210);
            }
            rhs[axis][row_lo] = b_lo;
        }
    }

    let sol: Vec<Vec<f64>> = rhs
        .iter()
        .map(|b| {
            solve_normal_equations(&a, b)
                .ok_or_else(|| format!("KP15 eq (7) ring singular for n = {n}"))
        })
        .collect::<Result<_, _>>()?;

    for r in 0..n {
        g[r][5][4] = Point3::new(sol[0][col_x54(r)], sol[1][col_x54(r)], sol[2][col_x54(r)]);
    }
    if p0_45_free {
        g[0][4][5] = free.p210;
    }
    for r in r_start..n {
        g[r][4][5] = Point3::new(sol[0][col_x45(r)], sol[1][col_x45(r)], sol[2][col_x45(r)]);
    }
    Ok(())
}

/// Solve the (S2) eq (8) ring for `p[r][4,4]` (`r = 1..n`); return
/// the deferred `sal = lyc[0]` residual.
pub(super) fn solve_eq8_ring(
    g: &mut [Grid],
    params: &Bi6Reparams,
    n: usize,
) -> Result<Point3, String> {
    let mu3 = mu3_full(params);
    let sf = |k: usize| k % n;
    // m_r(coord) — known part of lyc[r]; equation is
    // p[t][4,4] − p[r][4,4] = m_r.
    let m_r = |g: &[Grid], r: usize, coord: fn(Point3) -> f64| -> f64 {
        let t = sf(r + 1);
        mu3[0] * coord(g[r][6][6])
            + mu3[1] * coord(g[r][5][6])
            + mu3[2] * coord(g[r][6][5])
            + mu3[3] * coord(g[r][5][5])
            + mu3[4] * coord(g[r][6][4])
            + mu3[5] * coord(g[r][5][4])
            + mu3[6] * coord(g[r][6][2])
            + mu3[7] * coord(g[r][6][1])
            + mu3[8] * coord(g[r][6][0])
            + mu3[9] * (coord(g[r][4][2]) - coord(g[t][2][4]))
            + mu3[10] * (coord(g[r][4][1]) - coord(g[t][1][4]))
            + mu3[11] * (coord(g[r][4][0]) - coord(g[t][0][4]))
    };

    // Unknowns y44[r], r = 1..n, mapped to column r − 1. p[0][4,4]
    // = P220 (already set). Equation index r − 1 for r = 1..n.
    let ncols = n - 1;
    let col = |r: usize| r - 1;
    let mut a = vec![vec![0.0; ncols]; ncols];
    let coord_of = [|p: Point3| p.x, |p: Point3| p.y, |p: Point3| p.z];
    let p220 = g[0][4][4];
    let mut rhs = [vec![0.0; ncols], vec![0.0; ncols], vec![0.0; ncols]];

    for r in 1..n {
        let t = sf(r + 1);
        let eq = col(r);
        a[eq][col(r)] += -1.0;
        if t == 0 {
            // p[t][4,4] = P220 known -> move to RHS.
            for (axis, coord) in coord_of.iter().enumerate() {
                rhs[axis][eq] = m_r(g, r, *coord) - coord(p220);
            }
        } else {
            a[eq][col(t)] += 1.0;
            for (axis, coord) in coord_of.iter().enumerate() {
                rhs[axis][eq] = m_r(g, r, *coord);
            }
        }
    }

    let sol: Vec<Vec<f64>> = rhs
        .iter()
        .map(|b| {
            solve_dense(a.clone(), b.clone())
                .ok_or_else(|| format!("KP15 eq (8) ring singular for n = {n}"))
        })
        .collect::<Result<_, _>>()?;

    for r in 1..n {
        g[r][4][4] = Point3::new(sol[0][col(r)], sol[1][col(r)], sol[2][col(r)]);
    }

    // sal = lyc[0] re-evaluated with all p[·][4,4] now assigned:
    //   p[1][4,4] − (p[0][4,4] + m_0).
    let t1 = sf(1);
    Ok(Point3::new(
        g[t1][4][4].x - (g[0][4][4].x + m_r(g, 0, |p| p.x)),
        g[t1][4][4].y - (g[0][4][4].y + m_r(g, 0, |p| p.y)),
        g[t1][4][4].z - (g[0][4][4].z + m_r(g, 0, |p| p.z)),
    ))
}

/// Apply the (S3) `µ̃`/`µ̄` rules for `p[r][5,3]`, `p[t][3,5]`,
/// `p[t][3,4]`.
pub(super) fn apply_s3_rules(g: &mut [Grid], params: &Bi6Reparams, n: usize) {
    let mut_ = mut_full(params);
    let mub = mub_full(params);
    let sf = |k: usize| k % n;
    let mut p53 = vec![ZERO; n];
    let mut p35 = vec![ZERO; n];
    let mut p34 = vec![ZERO; n];

    // Common block: the µ-weighted combination shared by (S3)'s
    // three outputs. `mirror = false` is the `p53`/`p34` side
    // (+diff block, midpoint = p[r][5,4]); `mirror = true` is the
    // `p35` side (−diff block, midpoint = p[t][4,5]).
    let block = |g: &[Grid], r: usize, t: usize, a: usize, tab: &[f64; 11], mirror: bool| -> f64 {
        let mid = if mirror {
            comp(g[t][4][5], a)
        } else {
            comp(g[r][5][4], a)
        };
        let sign = if mirror { -1.0 } else { 1.0 };
        let d42 = comp(g[r][4][2], a) - comp(g[t][2][4], a);
        let d41 = comp(g[r][4][1], a) - comp(g[t][1][4], a);
        let d40 = comp(g[r][4][0], a) - comp(g[t][0][4], a);
        let d44 = comp(g[r][4][4], a) - comp(g[t][4][4], a);
        tab[0] * comp(g[r][6][6], a)
            + tab[1] * comp(g[r][6][5], a)
            + tab[2] * comp(g[r][6][4], a)
            + tab[3] * mid
            + tab[4] * comp(g[r][6][2], a)
            + tab[5] * comp(g[r][6][1], a)
            + tab[6] * comp(g[r][6][0], a)
            + sign * (tab[7] * d42 + tab[8] * d41 + tab[9] * d40 + tab[10] * d44)
    };

    for r in 0..n {
        let t = sf(r + 1);
        p53[r] = point_from_axes(|a| block(g, r, t, a, &mut_, false));
        p35[t] = point_from_axes(|a| block(g, r, t, a, &mut_, true));
        p34[t] = point_from_axes(|a| comp(g[r][4][3], a) + block(g, r, t, a, &mub, false));
    }
    for r in 0..n {
        g[r][5][3] = p53[r];
    }
    for t in 0..n {
        g[t][3][5] = p35[t];
        g[t][3][4] = p34[t];
    }
}
