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

use popkern_core::geometry::gregory::{DualJetCell, GregoryPatch};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};

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

    // Joint tilt solve (PK-521 rework of the two-pass design):
    // ONE least-squares over every registered rim side's depth-1 row,
    // fitting the measured tilt-nulling field with
    //   - ONE shared constant per corner cluster over the 3-CP clamp
    //     band. The constant is the invariant-safe maximum: it moves
    //     every band cell around the cluster by the same vector — a
    //     LOCAL TRANSLATION, which every cross-spoke G relation
    //     survives. Per-layer shared deltas were re-falsified here
    //     WITH the mechanism: layered deltas are not a translation,
    //     and hexagon6_mixed's live-seam transverse G2 blew to
    //     2.7e-1 (the same class as the PK-515 clamp falsification).
    //   - FOUR per-side interior unknowns following the tilt field's
    //     measured shape. The historical interior-fit falsification
    //     (octagon8_mixed cross-spoke 1.98) predates the PK-515/519
    //     substrate; with reconcile_rim_curvature re-truing depth-2
    //     against the strips afterwards and cluster constants kept
    //     translation-safe, the fixture matrix measures green.
    //   - TWO per-side unknowns at each dual-jet corner end, written
    //     into the side's own claims (the end CP itself — the other
    //     rim's boundary row — stays fixed).
    // Replacing the per-side averaged-endpoint + linear profile with
    // this joint fit removed the clamp-band tilt bump whose along-rim
    // derivative was the mid-band diagonal mixed-curvature floor:
    // flagship worst rim g2 4.877e-2 → 6.7e-3, worst g1
    // 1.6e-4 → 4.6e-6.
    // Depth application is unchanged: row 1 moves by δ and row 2 by
    // 2δ, keeping row2 − 2·row1 + row0 exact.
    const CORNER_CLUSTER: f64 = 1.0e-3;
    let ts: Vec<f64> = (0..M).map(|m| m as f64 / (M - 1) as f64).collect();
    struct JetSide {
        pi: usize,
        side: usize,
        field: Vec<[f64; 3]>,
        start_cluster: Option<usize>,
        end_cluster: Option<usize>,
        start_dual: Option<usize>,
        end_dual: Option<usize>,
        interior_base: usize,
        dual_base: usize,
    }
    impl JetSide {
        fn dual_indices(&self, end: bool) -> (usize, usize) {
            let offset = if end && self.start_dual.is_some() {
                2
            } else {
                0
            };
            (self.dual_base + offset, self.dual_base + offset + 1)
        }
    }
    let mut centers: Vec<Point3> = Vec::new();
    let mut cluster_of = |p: Point3| -> usize {
        centers
            .iter()
            .position(|c| c.distance_to(&p) < CORNER_CLUSTER)
            .unwrap_or_else(|| {
                centers.push(p);
                centers.len() - 1
            })
    };
    let mut jet_sides: Vec<JetSide> = Vec::new();
    for &(pi, side) in &sides {
        let (Some(row0), Some(row1)) = (row_of(patches, pi, side, 0), row_of(patches, pi, side, 1))
        else {
            continue;
        };
        // The measured field: at each station, the depth-1 delta that
        // projects the cross-jet into the support tangent plane. A
        // side whose measurement declines anywhere is skipped whole.
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
        jet_sides.push(JetSide {
            pi,
            side,
            field,
            start_cluster: (start_dual.is_none()).then(|| cluster_of(bezier9_eval(&row0, 0.0))),
            end_cluster: (end_dual.is_none()).then(|| cluster_of(bezier9_eval(&row0, 1.0))),
            start_dual,
            end_dual,
            interior_base: 0,
            dual_base: 0,
        });
    }
    if jet_sides.is_empty() {
        return;
    }
    let n_clusters = centers.len();
    let mut next_unknown = n_clusters;
    for js in jet_sides.iter_mut() {
        js.interior_base = next_unknown;
        next_unknown += 4;
        js.dual_base = next_unknown;
        next_unknown += 2 * (js.start_dual.is_some() as usize + js.end_dual.is_some() as usize);
    }
    let unknowns = next_unknown;
    let bern9 = |k: usize, t: f64| -> f64 {
        const C: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];
        C[k] * t.powi(k as i32) * (1.0 - t).powi((9 - k) as i32)
    };
    let mut ata = vec![vec![0.0_f64; unknowns]; unknowns];
    let mut atb = vec![[0.0_f64; 3]; unknowns];
    for js in &jet_sides {
        for (m, target) in js.field.iter().enumerate() {
            let t = ts[m];
            let mut entries: Vec<(usize, f64)> = Vec::with_capacity(10);
            if js.start_dual.is_some() {
                let (idx1, idx2) = js.dual_indices(false);
                entries.push((idx1, bern9(1, t)));
                entries.push((idx2, bern9(2, t)));
            } else if let Some(cluster) = js.start_cluster {
                let mass: f64 = (0..=2).map(|k| bern9(k, t)).sum();
                entries.push((cluster, mass));
            }
            if js.end_dual.is_some() {
                let (idx1, idx2) = js.dual_indices(true);
                entries.push((idx1, bern9(8, t)));
                entries.push((idx2, bern9(7, t)));
            } else if let Some(cluster) = js.end_cluster {
                let mass: f64 = (7..=9).map(|k| bern9(k, t)).sum();
                entries.push((cluster, mass));
            }
            for q in 0..4usize {
                entries.push((js.interior_base + q, bern9(3 + q, t)));
            }
            for &(i, wi) in &entries {
                for &(j, wj) in &entries {
                    ata[i][j] += wi * wj;
                }
                for c in 0..3 {
                    atb[i][c] += wi * target[c];
                }
            }
        }
    }
    let Some(solution) = super::dense_spd_solve::solve_dense_spd(&mut ata, &atb) else {
        return;
    };
    for js in &jet_sides {
        let JetSide { pi, side, .. } = *js;
        for k in 0..10usize {
            let dual = match k {
                0..=2 => js.start_dual,
                7..=9 => js.end_dual,
                _ => None,
            };
            let delta = match k {
                0..=2 => match (js.start_dual, js.start_cluster) {
                    (Some(_), _) => {
                        let (idx1, idx2) = js.dual_indices(false);
                        match k {
                            0 => Vec3::zero(),
                            1 => solution[idx1],
                            _ => solution[idx2],
                        }
                    }
                    (_, Some(cluster)) => solution[cluster],
                    _ => Vec3::zero(),
                },
                7..=9 => match (js.end_dual, js.end_cluster) {
                    (Some(_), _) => {
                        let (idx1, idx2) = js.dual_indices(true);
                        match k {
                            9 => Vec3::zero(),
                            8 => solution[idx1],
                            _ => solution[idx2],
                        }
                    }
                    (_, Some(cluster)) => solution[cluster],
                    _ => Vec3::zero(),
                },
                _ => solution[js.interior_base + (k - 3)],
            };
            // Dual-end routing: the end cell itself is the OTHER rim
            // side's d0/d1 row — fixed; one and two steps out, the
            // d1/d2 corrections go into this side's own claims.
            if let Some(corner) = dual {
                let at_end_pos = if k <= 2 { k } else { 9 - k };
                if at_end_pos == 0 {
                    continue;
                }
                let d1_cell = dual_end_d1_cells(side)[at_end_pos - 1];
                let d2_cell = dual_end_cells(side)[at_end_pos - 1];
                add_claim_delta(patches, pi, side, corner, d1_cell, delta);
                add_claim_delta(patches, pi, side, corner, d2_cell, delta * 2.0);
                continue;
            }
            // Row 1 carries the jet correction; row 2 moves by 2δ so
            // the cross-second-difference row2 − 2·row1 + row0 — the
            // piece's cross-curvature — is EXACTLY preserved (a
            // row-1-only move exploded octagon8_mixed's transverse
            // seam G2 to 1.78; with the invariant it is untouched).
            let (i1, j1) = coord_for(side, 1, k);
            let (i2, j2) = coord_for(side, 2, k);
            let Some(grid) = patches[pi].binonic_grid_mut() else {
                continue;
            };
            let g1 = grid[i1][j1];
            grid[i1][j1] = Point3::new(g1.x + delta.x, g1.y + delta.y, g1.z + delta.z);
            let g2 = grid[i2][j2];
            grid[i2][j2] = Point3::new(
                g2.x + 2.0 * delta.x,
                g2.y + 2.0 * delta.y,
                g2.z + 2.0 * delta.z,
            );
        }
        if popkern_core::trace::trace_enabled("setback") {
            // Residual attribution: the worst remaining tilt after
            // the joint fit (the PK-392 cap↔strip class lives here).
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
                for (m, f) in js.field.iter().enumerate() {
                    let p1 = bezier9_eval(&row1, ts[m]);
                    let p0 = bezier9_eval(&row0, ts[m]);
                    let j = p1 - p0;
                    let mag = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
                    pre = pre.max(mag / j.length().max(1e-15));
                }
            }
            eprintln!(
                "  [rim-jet] patch {pi} side {side}: joint tilt solve; tilt sin pre={pre:.3e} post={post:.3e} (post peak at t={post_at:.3})"
            );
        }
    }
}
