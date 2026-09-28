//! Spoke re-reconcile for the setback cap fan (PK-582).
//!
//! The rim passes land each rim on its strip. At a SPOKE junction (a
//! rim end meeting the next sector's rim across a spoke seam), the
//! landing needs a translation of the junction neighbourhood, and that
//! translation moves the spoke's own boundary CPs near the rim end
//! (the rim's depth-1/2 cells there ARE the spoke's boundary cells 1
//! and 2). The spoke's cross rows do not follow, so the spoke's G1
//! relation `c_A + λ c_B = μ t` broke by `μ Δt` near the junction:
//! pentagon5_uniform's five spokes certified g1 1.2e-5 against 3.4e-7
//! before the rim passes. Measured and rejected: dropping the shared
//! translation (the rim cannot land: flagship rim g2 4e-5 -> 2.6e-2),
//! and continuing it along the spoke as a per-station translation
//! (translating cross cells keeps `c` but not `t`: g1 7.6e-5).
//!
//! This pass restores the relation directly. It reads the spoke
//! pair's `λ(s)`, `μ(s)` pointwise from the fan as built before the
//! rim passes (where the relation held), then solves for the smallest
//! correction of both pieces' spoke depth-1 cells (their own claims
//! at the dual corner, the grid further in) that zeroes
//! `c'_A + λ c'_B − μ t'` along the rebuilt spoke. The spoke's
//! boundary row and the junction's rim cells are never written: the
//! corner cell needs no change, because the rim is C1 through the
//! junction, so `μ(0) = 0`.

use popkern_core::geometry::gregory::{GregoryPatch, dual_jet_cell_of};
use popkern_core::geometry::{Point3, SurfaceGeometry, Vec3};

use super::dense_spd_solve::VecLsq;
use super::setback_corner_dual_jet::{add_claim_delta, row_for_side};

/// Spoke depth-1 cells re-fitted, counted from the junction (the far
/// end belongs to the fan centre's own relations).
const REACH: usize = 6;
/// Samples along the spoke for the relation fit.
const SAMPLES: usize = 41;
/// Damping of the correction (min-norm among the relation's null
/// directions, which splits it between the two pieces).
const DAMPING: f64 = 1.0e-3;
/// Two spoke rows are one seam when their end points agree to this.
const SAME_POINT: f64 = 1.0e-6;

/// Measured outcome of [`reconcile_spokes`].
#[derive(Debug, Default, Clone, Copy)]
pub(super) struct SpokeReconcile {
    /// Spoke seams re-fitted.
    pub(super) spokes: usize,
    /// Worst relation residual `|c_A + λ c_B − μ t|` relative to `|t|`
    /// over the re-fitted spokes, before and after.
    pub(super) worst_before: f64,
    pub(super) worst_after: f64,
}

/// One piece's view of a spoke: patch, side, and whether its row index
/// runs away from the junction (`k = 0` at the junction).
#[derive(Clone, Copy)]
struct SpokeSide {
    pi: usize,
    side: usize,
    from_junction: bool,
}

fn bern(n: usize, k: usize, t: f64) -> f64 {
    let mut c = 1.0;
    for i in 0..k {
        c = c * (n - i) as f64 / (i + 1) as f64;
    }
    c * t.powi(k as i32) * (1.0 - t).powi((n - k) as i32)
}

/// A side's rows re-indexed from the junction.
fn rows(patches: &[GregoryPatch], s: SpokeSide) -> Option<([Point3; 10], [Point3; 10])> {
    let mut r0 = row_for_side(patches, s.pi, s.side, 0)?;
    let mut r1 = row_for_side(patches, s.pi, s.side, 1)?;
    if !s.from_junction {
        r0.reverse();
        r1.reverse();
    }
    Some((r0, r1))
}

/// Cross derivative `9 Σ (r1_k − r0_k) B_k(s)` and boundary tangent.
fn jets(r0: &[Point3; 10], r1: &[Point3; 10], s: f64) -> (Vec3, Vec3) {
    let mut c = Vec3::zero();
    let mut t = Vec3::zero();
    for k in 0..10 {
        c = c + (r1[k] - r0[k]) * (9.0 * bern(9, k, s));
    }
    for k in 0..9 {
        t = t + (r0[k + 1] - r0[k]) * (9.0 * bern(8, k, s));
    }
    (c, t)
}

/// `(λ, μ)` minimizing `|c_A + λ c_B − μ t|`.
fn relation(ca: Vec3, cb: Vec3, t: Vec3) -> Option<(f64, f64)> {
    // Normal equations of [c_B, −t]·(λ, μ) = −c_A.
    let (a11, a12, a22) = (cb.dot(&cb), -cb.dot(&t), t.dot(&t));
    let (b1, b2) = (-ca.dot(&cb), ca.dot(&t));
    let det = a11 * a22 - a12 * a12;
    if det.abs() <= 1e-24 * (a11 * a22).max(1e-300) {
        return None;
    }
    Some(((b1 * a22 - a12 * b2) / det, (a11 * b2 - a12 * b1) / det))
}

/// The spoke pairs at activated spoke-junction corners: each rim side
/// ending at a dual corner whose other side is NOT a rim side.
fn spoke_pairs(
    patches: &[GregoryPatch],
    rim_sides: &[(usize, usize)],
) -> Vec<(SpokeSide, SpokeSide)> {
    use super::setback_cluster_field::fan_side;
    use super::setback_corner_dual_jet::active_corner_at_end;
    let mut ends: Vec<(Point3, SpokeSide)> = Vec::new();
    for &(pi, side) in rim_sides {
        for at_k9 in [false, true] {
            let spoke = fan_side(side, at_k9);
            if active_corner_at_end(patches, pi, side, at_k9).is_none()
                || rim_sides.contains(&(pi, spoke))
            {
                continue;
            }
            let Some(rim0) = row_for_side(patches, pi, side, 0) else {
                continue;
            };
            let Some(spoke0) = row_for_side(patches, pi, spoke, 0) else {
                continue;
            };
            let junction = if at_k9 { rim0[9] } else { rim0[0] };
            let from_junction =
                spoke0[0].distance_to(&junction) <= spoke0[9].distance_to(&junction);
            ends.push((
                junction,
                SpokeSide {
                    pi,
                    side: spoke,
                    from_junction,
                },
            ));
        }
    }
    let mut pairs = Vec::new();
    for (i, (pa, a)) in ends.iter().enumerate() {
        for (pb, b) in ends.iter().skip(i + 1) {
            if a.pi == b.pi || pa.distance_to(pb) > SAME_POINT {
                continue;
            }
            let same_row = rows(patches, *a)
                .zip(rows(patches, *b))
                .is_some_and(|((ra, _), (rb, _))| ra[9].distance_to(&rb[9]) <= SAME_POINT);
            if same_row {
                pairs.push((*a, *b));
            }
        }
    }
    pairs
}

/// Write `delta` into a side's depth-1 cell `k` (junction-indexed):
/// the side's own claim in an activated dual cell, the grid otherwise.
fn add_d1(patches: &mut [GregoryPatch], s: SpokeSide, k: usize, delta: Vec3) {
    use super::setback_corner_dual_jet::coord_for;
    let kk = if s.from_junction { k } else { 9 - k };
    let (i, j) = coord_for(s.side, 1, kk);
    if let GregoryPatch::RectangularBiNonicDualJet(r) = &patches[s.pi]
        && let Some((corner, cell)) = dual_jet_cell_of(i, j)
        && r.corners[corner].is_some()
    {
        add_claim_delta(patches, s.pi, s.side, corner, cell, delta);
        return;
    }
    if let Some(grid) = patches[s.pi].binonic_grid_mut() {
        let g = grid[i][j];
        grid[i][j] = Point3::new(g.x + delta.x, g.y + delta.y, g.z + delta.z);
    }
}

/// Re-fit every spoke pair's depth-1 cells so the pair's pre-pass G1
/// relation holds along the rebuilt spoke (see the module docs).
pub(super) fn reconcile_spokes(
    patches: &mut [GregoryPatch],
    before: &[GregoryPatch],
    candidates: &[SurfaceGeometry],
    border_loop: &[Vec<Point3>],
    rim_off_loop_allowance: f64,
) -> SpokeReconcile {
    let rim_sides = super::setback_corner_dual_jet::rim_sides(
        patches,
        candidates,
        border_loop,
        rim_off_loop_allowance,
    );
    let mut report = SpokeReconcile::default();
    for (a, b) in spoke_pairs(patches, &rim_sides) {
        let (Some((a0, a1)), Some((b0, b1))) = (rows(before, a), rows(before, b)) else {
            continue;
        };
        let (Some((na0, na1)), Some((nb0, nb1))) = (rows(patches, a), rows(patches, b)) else {
            continue;
        };
        let mut lsq = VecLsq::new(2 * REACH);
        let mut worst = 0.0_f64;
        let mut rows_used = Vec::with_capacity(SAMPLES);
        for j in 0..SAMPLES {
            let s = j as f64 / (SAMPLES - 1) as f64;
            let (ca, t) = jets(&a0, &a1, s);
            let (cb, _) = jets(&b0, &b1, s);
            let Some((lambda, mu)) = relation(ca, cb, t) else {
                continue;
            };
            let (nca, nt) = jets(&na0, &na1, s);
            let (ncb, _) = jets(&nb0, &nb1, s);
            let residual = nca + ncb * lambda - nt * mu;
            worst = worst.max(residual.length() / nt.length().max(1e-300));
            let mut entries = Vec::with_capacity(2 * REACH);
            for k in 1..=REACH {
                let w = 9.0 * bern(9, k, s);
                entries.push((k - 1, w));
                entries.push((REACH + k - 1, lambda * w));
            }
            lsq.add_vec_row(&entries, residual * -1.0);
            rows_used.push((s, lambda, mu));
        }
        for i in 0..2 * REACH {
            lsq.add_vec_row(&[(i, DAMPING)], Vec3::zero());
        }
        let Some(x) = lsq.solve() else {
            continue;
        };
        for k in 1..=REACH {
            add_d1(patches, a, k, x[k - 1]);
            add_d1(patches, b, k, x[REACH + k - 1]);
        }
        let mut after = 0.0_f64;
        if let (Some((ra0, ra1)), Some((rb0, rb1))) = (rows(patches, a), rows(patches, b)) {
            for &(s, lambda, mu) in &rows_used {
                let (ca, t) = jets(&ra0, &ra1, s);
                let (cb, _) = jets(&rb0, &rb1, s);
                let r = ca + cb * lambda - t * mu;
                after = after.max(r.length() / t.length().max(1e-300));
            }
        }
        report.spokes += 1;
        report.worst_before = report.worst_before.max(worst);
        report.worst_after = report.worst_after.max(after);
    }
    report
}
