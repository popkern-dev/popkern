//! Rim cross-jet reconciliation (PK-212): after the G0 rim solve
//! lands the fan's rim rows on the live supports, project each rim
//! row's CROSS-JET (depth-1 row) into the support's tangent plane at
//! the landed border, so the cap meets the strips G1.
//!
//! Why this exists: the border derivation's tangent-plane
//! prolongation (kp20_paper_border_extract) fixes the dominant
//! chord-tilt term at the four control points of each cubic arc
//! piece, but tangency is needed at every EVALUATED point — the
//! per-CP tangent planes rotate along the arc and their cubic
//! combination leaves a smooth residual (measured 4.5° on the blunt
//! pentagon under G2 sweeps). The committed pieces are bi-9: their
//! ten-CP depth-1 rows have the freedom to null the residual field.
//!
//! Cross-seam discipline: at a shared fan corner the spoke relations
//! tie both pieces' station-adjacent columns, so the three CPs
//! nearest each end take ONE shared per-corner delta (the average of
//! the incident sides' endpoint targets, clustered by distance for
//! skew-trim corners), and the interior blends LINEARLY between the
//! two — zero second difference along the row except at the clamp
//! joints. Every k-varying interior alternative was falsified by
//! measurement: a per-column least-squares fit and a one-DOF smooth
//! bump both wreck the fan's cross-spoke transverse G2 through the
//! twist they introduce (1.98 and full postcondition declines
//! respectively). Depth application is curvature-invariant: row 1
//! moves by δ, row 2 by 2δ, keeping row2 − 2·row1 + row0 exact.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{DualJetCell, GregoryPatch};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;

use super::bezier9_eval;
use super::setback_corner_dual_jet::{
    RIM_NEAR, active_corner_at_end, add_claim_delta, coord_for, dual_end_cells, rim_sides,
    row_for_side, side_is_row,
};

const M: usize = 64;

fn row_of(patches: &[GregoryPatch], patch: usize, side: usize, d: usize) -> Option<[Point3; 10]> {
    row_for_side(patches, patch, side, d)
}

/// The claim cells a side's depth-1 row addresses at a dual end,
/// ordered by distance from the corner: the twist cell, then the
/// cell holding this side's d1 band two steps out.
fn dual_end_d1_cells(side: usize) -> [DualJetCell; 2] {
    if side_is_row(side) {
        [DualJetCell::Cell11, DualJetCell::Cell12]
    } else {
        [DualJetCell::Cell11, DualJetCell::Cell21]
    }
}

/// Support normal at the landed point `p`: nearest candidate surface
/// within the rim band, finite-difference normal at its projection.
fn support_normal(
    candidates: &[SurfaceGeometry],
    p: &Point3,
) -> Option<popkern_core::geometry::Vec3> {
    let mut best: Option<(f64, usize, (f64, f64))> = None;
    for (idx, sfc) in candidates.iter().enumerate() {
        let Some((u, v)) = project_point_to_uv(sfc, p) else {
            continue;
        };
        let d = sfc.point_at(u, v).distance_to(p);
        if d < RIM_NEAR && best.map(|(bd, _, _)| d < bd).unwrap_or(true) {
            best = Some((d, idx, (u, v)));
        }
    }
    let (_, idx, (u, v)) = best?;
    let sfc = &candidates[idx];
    let h = 1e-6;
    let du = sfc.point_at(u + h, v) - sfc.point_at(u - h, v);
    let dv = sfc.point_at(u, v + h) - sfc.point_at(u, v - h);
    let (lu, lv) = (du.length(), dv.length());
    if lu.min(lv) < 1e-3 * lu.max(lv) {
        return None;
    }
    du.cross(&dv).normalized()
}

/// Run the rim cross-jet reconciliation over the committed fan.
pub(super) fn reconcile_rim_jets(
    patches: &mut [GregoryPatch],
    candidates: &[SurfaceGeometry],
    strip_count: usize,
    border_loop: &[Vec<Point3>],
    rim_off_loop_allowance: f64,
) {
    // The tangency targets come from the STRIPS only (the first
    // `strip_count` candidates): the setback arcs lie on the strips,
    // and near the pulled-back corners the LATERAL PLANES can sit
    // nearer in space — mixing their normals into the corner targets
    // pollutes the jet field with the wrong plane.
    let strips = &candidates[..strip_count.min(candidates.len())];
    // Rim-side registry (same gates as the G0 solve: near the
    // supports along the whole row + on the border loop).
    let sides = rim_sides(patches, candidates, border_loop, rim_off_loop_allowance);
    if sides.is_empty() {
        return;
    }

    // Pass 1: per-corner shared endpoint targets, clustered by
    // DISTANCE — on skew-trim fixtures (octagon8_mixed) adjacent
    // sides' corner endpoints are closest-approach midpoints that do
    // NOT coincide bit-wise, and grid-keyed grouping left each side
    // its own corner value, breaking the shared-column contract
    // (measured: cross-spoke transverse G2 1.99).
    let mut cluster_centers: Vec<Point3> = Vec::new();
    let mut corner_acc: Vec<(Point3, usize)> = Vec::new();
    const CORNER_CLUSTER: f64 = 1.0e-3;
    let endpoint_target =
        |patches: &[GregoryPatch], pi: usize, side: usize, t: f64| -> Option<(Point3, Point3)> {
            let row0 = row_of(patches, pi, side, 0)?;
            let row1 = row_of(patches, pi, side, 1)?;
            let p0 = bezier9_eval(&row0, t);
            let p1 = bezier9_eval(&row1, t);
            let n = support_normal(strips, &p0)?;
            let j = p1 - p0;
            let d = n * (-j.dot(&n));
            Some((p0, Point3::new(d.x, d.y, d.z)))
        };
    let cluster_of =
        |p: &Point3, centers: &mut Vec<Point3>, acc: &mut Vec<(Point3, usize)>| -> usize {
            for (i, c) in centers.iter().enumerate() {
                if c.distance_to(p) < CORNER_CLUSTER {
                    return i;
                }
            }
            centers.push(*p);
            acc.push((Point3::new(0.0, 0.0, 0.0), 0));
            centers.len() - 1
        };
    for &(pi, side) in &sides {
        for t in [0.0, 1.0] {
            // A dual-jet corner end never joins the shared-average
            // cluster: its correction is per-side through the claims.
            if active_corner_at_end(patches, pi, side, t > 0.5).is_some() {
                continue;
            }
            if let Some((p0, d)) = endpoint_target(patches, pi, side, t) {
                let i = cluster_of(&p0, &mut cluster_centers, &mut corner_acc);
                let e = &mut corner_acc[i];
                e.0 = Point3::new(e.0.x + d.x, e.0.y + d.y, e.0.z + d.z);
                e.1 += 1;
            }
        }
    }
    let corner_avg = |p: &Point3| -> [f64; 3] {
        cluster_centers
            .iter()
            .enumerate()
            .filter(|(_, c)| c.distance_to(p) < CORNER_CLUSTER)
            .filter(|(i, _)| corner_acc[*i].1 > 0)
            .map(|(i, _)| {
                let (sum, n) = &corner_acc[i];
                let n = *n as f64;
                [sum.x / n, sum.y / n, sum.z / n]
            })
            .next()
            .unwrap_or([0.0; 3])
    };

    // Pass 2: per-side damped LSQ for the depth-1 row deltas, pinned
    // to the shared corner targets at both ends.
    let ts: Vec<f64> = (0..M).map(|m| m as f64 / (M - 1) as f64).collect();
    for &(pi, side) in &sides {
        let Some(row0) = row_of(patches, pi, side, 0) else {
            continue;
        };
        let Some(row1) = row_of(patches, pi, side, 1) else {
            continue;
        };
        let mut field: Vec<[f64; 3]> = Vec::with_capacity(M);
        let mut ok = true;
        for &t in &ts {
            let p0 = bezier9_eval(&row0, t);
            let p1 = bezier9_eval(&row1, t);
            let Some(n) = support_normal(strips, &p0) else {
                ok = false;
                break;
            };
            let j = p1 - p0;
            let d = n * (-j.dot(&n));
            field.push([d.x, d.y, d.z]);
        }
        if !ok {
            continue;
        }
        let start_dual = active_corner_at_end(patches, pi, side, false);
        let end_dual = active_corner_at_end(patches, pi, side, true);
        // At a dual end the side takes its OWN endpoint target (the
        // shared average is exactly the compromise the dual claims
        // exist to remove); elsewhere the shared cluster average.
        let own = |t: f64| -> [f64; 3] {
            endpoint_target(patches, pi, side, t)
                .map(|(_, d)| [d.x, d.y, d.z])
                .unwrap_or([0.0; 3])
        };
        let d_start = if start_dual.is_some() {
            own(0.0)
        } else {
            corner_avg(&bezier9_eval(&row0, 0.0))
        };
        let d_end = if end_dual.is_some() {
            own(1.0)
        } else {
            corner_avg(&bezier9_eval(&row0, 1.0))
        };
        // LINEAR-in-k deltas between the shared corner values: zero
        // second difference along the row (twist-neutral), and the
        // 2δ depth application keeps the cross-second-difference
        // invariant — an interior least-squares fit was measured to
        // explode octagon8_mixed's cross-spoke transverse G2 to
        // 1.98 through exactly the twist it introduced. The residual
        // this leaves is the field's nonlinear part; the tilt field
        // is near-uniform (the chord-tilt mechanism), so the linear
        // component carries it under the commercial bar.
        // Clamped profile: the three columns nearest each end take
        // the SHARED corner delta (the spoke relations tie those
        // columns across the seam — a linear-through-corners profile
        // was measured to crease the cap-cap transverse G2 to 0.55);
        // the interior blends linearly, so its second differences
        // vanish everywhere except the two clamp joints, which sit
        // on interior columns away from any spoke seam.
        let mut deltas = [[0.0_f64; 3]; 10];
        for (k, d) in deltas.iter_mut().enumerate() {
            let f = ((k as f64 - 2.0) / 5.0).clamp(0.0, 1.0);
            for c in 0..3 {
                d[c] = d_start[c] * (1.0 - f) + d_end[c] * f;
            }
        }
        for (k, delta) in deltas.iter().enumerate() {
            // Row 1 carries the jet correction; row 2 moves by 2δ so
            // the cross-second-difference row2 − 2·row1 + row0 — the
            // piece's cross-curvature — is EXACTLY preserved (a
            // row-1-only move exploded octagon8_mixed's transverse
            // seam G2 to 1.78; with the invariant it is untouched).
            //
            // Dual-end routing: the end cell itself is the OTHER rim
            // side's d0/d1 row — fixed (the shared-average write that
            // used to land there was the compromise); one and two
            // steps out, the d1/d2 corrections go into this side's
            // own claims.
            let dual = match k {
                0..=2 => start_dual,
                7..=9 => end_dual,
                _ => None,
            };
            let dv = popkern_core::geometry::Vec3::new(delta[0], delta[1], delta[2]);
            if let Some(corner) = dual {
                let at_end_pos = if k <= 2 { k } else { 9 - k };
                if at_end_pos == 0 {
                    continue;
                }
                if at_end_pos <= 2 {
                    let d1_cell = dual_end_d1_cells(side)[at_end_pos - 1];
                    let d2_cell = dual_end_cells(side)[at_end_pos - 1];
                    add_claim_delta(patches, pi, side, corner, d1_cell, dv);
                    add_claim_delta(patches, pi, side, corner, d2_cell, dv * 2.0);
                    continue;
                }
            }
            let (i1, j1) = coord_for(side, 1, k);
            let (i2, j2) = coord_for(side, 2, k);
            let Some(grid) = patches[pi].binonic_grid_mut() else {
                continue;
            };
            let g1 = grid[i1][j1];
            grid[i1][j1] = Point3::new(g1.x + delta[0], g1.y + delta[1], g1.z + delta[2]);
            let g2 = grid[i2][j2];
            grid[i2][j2] = Point3::new(
                g2.x + 2.0 * delta[0],
                g2.y + 2.0 * delta[1],
                g2.z + 2.0 * delta[2],
            );
        }
        if popkern_core::trace::trace_enabled("setback") {
            // Residual attribution: the linear-in-k profile nulls the
            // tilt field's linear component only — measure what it
            // leaves (the PK-392 cap↔strip class lives here).
            let (mut pre, mut post) = (0.0_f64, 0.0_f64);
            let mut post_at = 0.0_f64;
            if let (Some(row0), Some(row1)) =
                (row_of(patches, pi, side, 0), row_of(patches, pi, side, 1))
            {
                for &t in &ts {
                    let p0 = bezier9_eval(&row0, t);
                    let p1 = bezier9_eval(&row1, t);
                    if let Some(n) = support_normal(strips, &p0) {
                        let j = p1 - p0;
                        let tilt = j.dot(&n).abs() / j.length().max(1e-15);
                        if tilt > post {
                            post = tilt;
                            post_at = t;
                        }
                    }
                }
                for (m, f) in field.iter().enumerate() {
                    let p0 = bezier9_eval(&row0, ts[m]);
                    let p1 = bezier9_eval(&row1, ts[m]);
                    let j = p1 - p0;
                    let mag = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
                    let _ = p0;
                    pre = pre.max(mag / j.length().max(1e-15));
                }
            }
            eprintln!(
                "  [rim-jet] patch {pi} side {side}: corner deltas applied; tilt sin pre={pre:.3e} post={post:.3e} (post peak at t={post_at:.3})"
            );
        }
    }
}
