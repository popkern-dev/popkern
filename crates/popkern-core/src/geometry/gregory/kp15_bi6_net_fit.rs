//! PK-087 stage 5 — fit the shared cap-net control ring to the real
//! fillet tensor border, inverting the B-spline → Bézier shrinkage.
//!
//! # Why
//!
//! The coarse [`convert_polar_to_bi6_cap_net`](super::convert_polar_to_bi6_cap_net)
//! ring places control symbols directly at geometric samples
//! (spoke samples, rim corners). The `BtoBB`
//! ([`bspline_net_to_bi3_bb`](super::kp15_bi6_deep_interior::...)) B-spline →
//! Bézier averaging then pulls the assembled cap's outer boundary
//! *inward* from those samples (the corner `b[0][0]` is a weighted
//! average dominated by the interior control point `a5`), so the cap
//! does not interpolate the surrounding surface: `cap_outer_g0 ≫ 5 mm`.
//!
//! # What
//!
//! The cap's outer-boundary row (script row 0) is a pure degree raise
//! of the `BtoBB` boundary coefficients `b[i][0]`. Those are a **linear**
//! function of the shared ring symbols `a1..a6, E`
//! ([`shared_net_layout`](super::shared_net_layout) ∘ `BtoBB`). So we
//! solve, per coordinate, for the ring that makes
//! `b[i][0] = tb.sectors[s][0][i]` — the real bi-3 tensor-border row 0
//! from [`extract_bi3_tensor_border_from_polar_input`](super::super::kp15_bi6_tensor_border_extract::extract_bi3_tensor_border_from_polar_input).
//! Because the sharing is preserved, the assembled net is still one
//! C²-consistent B-spline surface, so every eq (2)/(3) cross-sector
//! constraint still closes by construction (oracle-proven, independent
//! of the control values). The boundary interpolation and the
//! continuity closure hold **simultaneously** — the trade-off the
//! coarse ring could not escape.
//!
//! # Regularisation (fold containment)
//!
//! The ring → boundary map is well-conditioned (cond ≈ 60 at `n ∈
//! {5, 6}`) but the boundary corner is dominated by the interior
//! symbol `a5`, so an *exact* interpolation moves interior control
//! points far enough to fold the deep interior. We therefore solve a
//! Tikhonov least-squares centred on the coarse ring: minimise
//! `‖A·Δ − r‖² + λ²‖Δ‖²`, `r = target − A·coarse`, and take
//! `ring = coarse + Δ`. `λ` trades residual boundary error against
//! ring deviation; the caller sweeps down until the cap is fold-free
//! *and* within the cap-outer gate, and declines with a structured
//! diagnostic if no `λ` clears both (mechanism honesty).

use super::super::kp_polar_boundary_input::{PolarBoundaryInput, TensorBorder};
use super::super::kp15_bi6_deep_interior::bspline_net_to_bi3_bb;
use super::super::kp15_bi6_tensor_border_extract::extract_bi3_tensor_border_from_polar_input;
use super::{Bi6CapNetInputs, SharedControlRing, build_coarse_ring_and_free, shared_net_layout};
use crate::geometry::Point3;

/// Per-sector explicit fit target: each sector's two star-opposite
/// outer edges as exact bi-3 (cubic) BB-CP rows. `target[s][0]` is
/// sector `s`'s `tad` edge (script `g[s][i][0]`), `target[s][1]` its
/// `tak` edge (`g[s][0][i]`), each 4 CPs oriented from the mid-sector
/// junction toward the spoke foot. Supplied by the fillet wire from
/// the real `boundary.segments` (`FilletEnd`/`FaceEdge`) curves, which
/// are cubic-exact, so the assembled cap's outer boundary interpolates
/// the surrounding surface rather than a crude spoke-derived proxy.
///
/// A transparent alias for `[[Point3; 4]; 2]`; callers outside this
/// crate build the concrete array directly (the alias is intentionally
/// not re-exported at the crate root, to keep the frozen `gregory`
/// module registration untouched).
pub type SectorOuterEdgeTarget = [[Point3; 4]; 2];

/// Convert a [`PolarBoundaryInput`] into the shared cap-net frame with
/// the control ring **fitted** to the fillet's outer boundary, so the
/// assembled cap interpolates the surrounding surface *and* the
/// eq (2)/(3) cross-sector constraints still close by construction.
///
/// `lambda` is the Tikhonov regularisation weight (see
/// [`fit_shared_ring_to_border`]): `0` interpolates the boundary
/// exactly (larger interior deviation → fold risk), larger `λ` keeps
/// the ring near the fold-free coarse seed at the cost of residual
/// boundary error.
///
/// `target`: `None` falls back to the spoke-derived tensor-border
/// extractor (retained for the core trade-off tests, whose crude
/// on-curve degree reduction is ~2.5e-2 off the real edge);
/// `Some(target)` fits the two outer edges of each sector to the
/// supplied real-boundary bi-3 rows (the fillet wire's exact
/// `boundary.segments` extraction — PK-087 stage 6).
///
/// # Errors
/// - `build_coarse_ring_and_free` fails (`polar.validate()`, `n < 5`).
/// - The tensor-border extraction fails (`None` path only).
/// - `target.len() != n` (`Some` path).
/// - The regularised solve is singular (`λ = 0` on a rank-deficient map).
pub fn convert_polar_to_bi6_cap_fitted(
    polar: &PolarBoundaryInput,
    lambda: f64,
    target: Option<&[SectorOuterEdgeTarget]>,
) -> Result<Bi6CapNetInputs, String> {
    let (params, coarse, free) = build_coarse_ring_and_free(polar)?;
    let n = params.n;
    let target_edges: Vec<SectorOuterEdgeTarget> = if let Some(t) = target {
        if t.len() != n {
            return Err(format!("explicit target length {} ≠ n = {n}", t.len()));
        }
        t.to_vec()
    } else {
        let TensorBorder::Bi3(tb) = extract_bi3_tensor_border_from_polar_input(polar)? else {
            return Err("bi-6 shared-net fit requires a bi-3 tensor border".to_string());
        };
        if tb.sectors.len() != n {
            return Err(format!(
                "tensor-border sector count {} ≠ n = {n}",
                tb.sectors.len()
            ));
        }
        // Legacy fallback: the crude extractor shares one row-0 between
        // adjacent sectors (tad of s, tak of s−1), reproduced here.
        (0..n)
            .map(|s| [tb.sectors[s][0], tb.sectors[(s + n - 1) % n][0]])
            .collect()
    };
    let ring = fit_shared_ring_to_border(&coarse, &target_edges, n, lambda)
        .ok_or_else(|| "shared-net boundary fit solve was singular".to_string())?;
    let net: Vec<[[Point3; 4]; 4]> = (0..n).map(|s| shared_net_layout(&ring, s, n)).collect();
    Ok((params, net, free))
}

/// Number of scalar unknowns in the shared ring for valence `n`:
/// six per-spoke symbol arrays (`a1..a6`) plus the single global `E`.
fn n_unknowns(n: usize) -> usize {
    6 * n + 1
}

/// Read scalar ring unknown `u` from `ring` (coordinate `coord`).
fn ring_get(ring: &SharedControlRing, u: usize, n: usize, coord: usize) -> f64 {
    let comp = |p: Point3| match coord {
        0 => p.x,
        1 => p.y,
        _ => p.z,
    };
    if u == 6 * n {
        return comp(ring.e);
    }
    let grp = u / n;
    let k = u % n;
    comp(match grp {
        0 => ring.a1[k],
        1 => ring.a2[k],
        2 => ring.a3[k],
        3 => ring.a4[k],
        4 => ring.a5[k],
        _ => ring.a6[k],
    })
}

/// Build a unit ring with scalar unknown `u` set to `1` on the `x`
/// channel and everything else zero — used to read the linear map's
/// coefficient column.
fn unit_ring(u: usize, n: usize) -> SharedControlRing {
    let zero = Point3::new(0.0, 0.0, 0.0);
    let mut ring = SharedControlRing {
        a1: vec![zero; n],
        a2: vec![zero; n],
        a3: vec![zero; n],
        a4: vec![zero; n],
        a5: vec![zero; n],
        a6: vec![zero; n],
        e: zero,
    };
    let one = Point3::new(1.0, 0.0, 0.0);
    if u == 6 * n {
        ring.e = one;
    } else {
        let grp = u / n;
        let k = u % n;
        match grp {
            0 => ring.a1[k] = one,
            1 => ring.a2[k] = one,
            2 => ring.a3[k] = one,
            3 => ring.a4[k] = one,
            4 => ring.a5[k] = one,
            _ => ring.a6[k] = one,
        }
    }
    ring
}

/// Number of scalar boundary constraints for valence `n`: each sector
/// contributes both of its star-opposite outer edges — the `tad` edge
/// `b[i][0]` (spoke `s`) and the `tak` edge `b[0][i]` (spoke `s − 1`),
/// four BB-CPs each.
fn n_constraints(n: usize) -> usize {
    8 * n
}

/// The design matrix `A` of the boundary map. Row `s * 8 + i` (`i ∈
/// 0..4`) is sector `s`'s `tad`-edge `BtoBB` coefficient `b[i][0]`; row
/// `s * 8 + 4 + i` is its `tak`-edge coefficient `b[0][i]`. Column `u`
/// is ring unknown `u`. The map is coordinate-independent, so a single
/// scalar matrix serves all three channels.
///
/// Fitting **both** outer edges (rather than only `tad`) is what makes
/// the cap interpolate the surrounding surface on *both* of each
/// sector's boundary edges: a `tad`-only fit nails one edge exactly but
/// lets the other drift (and needs a far larger, fold-prone ring move),
/// whereas the balanced two-edge fit shares the residual evenly and
/// stays near the fold-free coarse ring.
fn boundary_design_matrix(n: usize) -> Vec<Vec<f64>> {
    let nu = n_unknowns(n);
    let mut a = vec![vec![0.0; nu]; n_constraints(n)];
    // `u` selects the column (ring unknown) being probed; the loop fills
    // that column across every constraint row, so index-by-`u` is the
    // natural form here.
    #[allow(clippy::needless_range_loop)]
    for u in 0..nu {
        let ring = unit_ring(u, n);
        for s in 0..n {
            let net = shared_net_layout(&ring, s, n);
            let b = bspline_net_to_bi3_bb(&net);
            for i in 0..4 {
                a[s * 8 + i][u] = b[i][0].x; // tad edge: b[i][0]
                a[s * 8 + 4 + i][u] = b[0][i].x; // tak edge: b[0][i]
            }
        }
    }
    a
}

/// The fit target for constraint row `row`, coordinate `coord`.
/// `tad` rows (`local ∈ 0..4`) interpolate sector `s`'s `tad` outer
/// edge (`target[s][0]`); `tak` rows (`local ∈ 4..8`) its `tak` outer
/// edge (`target[s][1]`). Both edges belong to sector `s` — matching
/// the real 2n-segment fillet boundary, where each cap sub-patch owns
/// its own `FilletEnd` and `FaceEdge` meeting at the mid-sector
/// junction (unlike the legacy shared-rim proxy, which fit `tak` from
/// the previous sector).
fn constraint_target(
    target_edges: &[SectorOuterEdgeTarget],
    row: usize,
    _n: usize,
    coord: usize,
) -> f64 {
    let comp = |p: Point3| match coord {
        0 => p.x,
        1 => p.y,
        _ => p.z,
    };
    let s = row / 8;
    let local = row % 8;
    let (edge, i) = if local < 4 {
        (0, local)
    } else {
        (1, local - 4)
    };
    comp(target_edges[s][edge][i])
}

/// Solve the Tikhonov-regularised normal equations
/// `(AᵀA + λ²I)·Δ = Aᵀ·r` for one coordinate channel.
fn solve_tikhonov(a: &[Vec<f64>], r: &[f64], lambda: f64, nu: usize) -> Option<Vec<f64>> {
    // Normal matrix AᵀA + λ²I.
    let mut ata = vec![vec![0.0; nu]; nu];
    for (col_i, row) in ata.iter_mut().enumerate() {
        for (col_j, cell) in row.iter_mut().enumerate() {
            let mut acc = 0.0;
            for a_row in a {
                acc += a_row[col_i] * a_row[col_j];
            }
            *cell = acc;
            if col_i == col_j {
                *cell += lambda * lambda;
            }
        }
    }
    // Right-hand side Aᵀ·r.
    let mut atr = vec![0.0; nu];
    for (u, atr_u) in atr.iter_mut().enumerate() {
        let mut acc = 0.0;
        for (row_idx, a_row) in a.iter().enumerate() {
            acc += a_row[u] * r[row_idx];
        }
        *atr_u = acc;
    }
    super::super::direct_linalg::solve_dense(ata, atr)
}

/// Fit the shared ring to the two real outer-edge bi-3 rows of each
/// sector.
///
/// `target_edges[s][e][i]` is the surrounding surface's boundary BB-CP
/// the cap must interpolate along sector `s`'s `tad` (`e = 0`) / `tak`
/// (`e = 1`) outer edge. `lambda` is the Tikhonov weight (≥ 0; `0` =
/// exact min-norm-from-coarse least squares). Returns the fitted ring
/// (coarse + regularised correction), or `None` if the normal-equation
/// solve is singular (should not happen for `λ > 0`).
pub fn fit_shared_ring_to_border(
    coarse: &SharedControlRing,
    target_edges: &[SectorOuterEdgeTarget],
    n: usize,
    lambda: f64,
) -> Option<SharedControlRing> {
    let nu = n_unknowns(n);
    let n_rows = n_constraints(n);
    let a = boundary_design_matrix(n);

    // Per-coordinate residual r = target − A·coarse, then Δ, then ring.
    let mut fitted = coarse.clone();
    for coord in 0..3 {
        // A·coarse for this channel.
        let coarse_vec: Vec<f64> = (0..nu).map(|u| ring_get(coarse, u, n, coord)).collect();
        let mut a_coarse = vec![0.0; n_rows];
        for (row_idx, a_row) in a.iter().enumerate() {
            let mut acc = 0.0;
            for (u, &c) in coarse_vec.iter().enumerate() {
                acc += a_row[u] * c;
            }
            a_coarse[row_idx] = acc;
        }
        // Residual against the target.
        let mut r = vec![0.0; n_rows];
        for (row_idx, r_cell) in r.iter_mut().enumerate() {
            *r_cell = constraint_target(target_edges, row_idx, n, coord) - a_coarse[row_idx];
        }
        let delta = solve_tikhonov(&a, &r, lambda, nu)?;
        let solved: Vec<f64> = (0..nu)
            .map(|u| ring_get(coarse, u, n, coord) + delta[u])
            .collect();
        write_channel(&mut fitted, &solved, n, coord);
    }
    Some(fitted)
}

/// Write one coordinate channel of a solved unknown vector back into
/// the ring.
fn write_channel(ring: &mut SharedControlRing, x: &[f64], n: usize, coord: usize) {
    let set = |p: &mut Point3, v: f64| match coord {
        0 => p.x = v,
        1 => p.y = v,
        _ => p.z = v,
    };
    for (u, &v) in x.iter().enumerate() {
        if u == 6 * n {
            set(&mut ring.e, v);
            continue;
        }
        let grp = u / n;
        let k = u % n;
        match grp {
            0 => set(&mut ring.a1[k], v),
            1 => set(&mut ring.a2[k], v),
            2 => set(&mut ring.a3[k], v),
            3 => set(&mut ring.a4[k], v),
            4 => set(&mut ring.a5[k], v),
            _ => set(&mut ring.a6[k], v),
        }
    }
}
