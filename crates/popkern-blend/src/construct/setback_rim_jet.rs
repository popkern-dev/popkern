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
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;

use super::bezier9_eval;
use super::setback_rim_reconcile::distance_to_polylines;

const M: usize = 64;
const RIM_NEAR: f64 = 2.0e-3;

fn coord_for(side: usize, d: usize, k: usize) -> (usize, usize) {
    match side {
        0 => (d, k),
        1 => (9 - d, k),
        2 => (k, d),
        _ => (k, 9 - d),
    }
}

fn row_of(patches: &[GregoryPatch], patch: usize, side: usize, d: usize) -> Option<[Point3; 10]> {
    let GregoryPatch::RectangularBiNonic(bp) = &patches[patch] else {
        return None;
    };
    Some(core::array::from_fn(|k| {
        let (i, j) = coord_for(side, d, k);
        bp.grid[i][j]
    }))
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
    let mut sides: Vec<(usize, usize)> = Vec::new();
    for (pi, p) in patches.iter().enumerate() {
        let GregoryPatch::RectangularBiNonic(bp) = p else {
            continue;
        };
        for side in 0..4 {
            let row: [Point3; 10] = core::array::from_fn(|k| {
                let (i, j) = coord_for(side, 0, k);
                bp.grid[i][j]
            });
            let on_support = (0..=16).all(|i| {
                let pt = bezier9_eval(&row, i as f64 / 16.0);
                candidates.iter().any(|sfc| {
                    project_point_to_uv(sfc, &pt)
                        .map(|(u, v)| sfc.point_at(u, v).distance_to(&pt) < RIM_NEAR)
                        .unwrap_or(false)
                })
            });
            if !on_support {
                continue;
            }
            let worst_off_loop = (0..=16)
                .map(|i| distance_to_polylines(bezier9_eval(&row, i as f64 / 16.0), border_loop))
                .fold(0.0_f64, f64::max);
            if worst_off_loop > rim_off_loop_allowance.max(RIM_NEAR) {
                continue;
            }
            sides.push((pi, side));
        }
    }
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
        let d_start = corner_avg(&bezier9_eval(&row0, 0.0));
        let d_end = corner_avg(&bezier9_eval(&row0, 1.0));
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
        let GregoryPatch::RectangularBiNonic(bp) = &mut patches[pi] else {
            continue;
        };
        for (k, delta) in deltas.iter().enumerate() {
            // Row 1 carries the jet correction; row 2 moves by 2δ so
            // the cross-second-difference row2 − 2·row1 + row0 — the
            // piece's cross-curvature — is EXACTLY preserved (a
            // row-1-only move exploded octagon8_mixed's transverse
            // seam G2 to 1.78; with the invariant it is untouched).
            let (i1, j1) = coord_for(side, 1, k);
            let g1 = bp.grid[i1][j1];
            bp.grid[i1][j1] = Point3::new(g1.x + delta[0], g1.y + delta[1], g1.z + delta[2]);
            let (i2, j2) = coord_for(side, 2, k);
            let g2 = bp.grid[i2][j2];
            bp.grid[i2][j2] = Point3::new(
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
            if let (Some(row0), Some(row1)) =
                (row_of(patches, pi, side, 0), row_of(patches, pi, side, 1))
            {
                for &t in &ts {
                    let p0 = bezier9_eval(&row0, t);
                    let p1 = bezier9_eval(&row1, t);
                    if let Some(n) = support_normal(strips, &p0) {
                        let j = p1 - p0;
                        let tilt = j.dot(&n).abs() / j.length().max(1e-15);
                        post = post.max(tilt);
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
                "  [rim-jet] patch {pi} side {side}: corner deltas applied;                  tilt sin pre={pre:.3e} post={post:.3e}"
            );
        }
    }
}
