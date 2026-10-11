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
use popkern_core::geometry::{Point3, Vec3};

use super::bezier9_eval;
use popkern_core::geometry::uv_projection::project_point_to_uv;

use super::dense_spd_solve::VecLsq;
use super::setback_cluster_field::{
    CLUSTER_FIELD_UNKNOWNS, ClusterFields, fan_side, spoke_dual_end,
};
use super::setback_corner_dual_jet::{
    RIM_NEAR, active_corner_at_end, add_claim_delta, coord_for, dual_end_cells, rim_sides,
    row_for_side, side_is_row,
};
use super::setback_rim_stations::side_support;

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
/// within the rim band, its exact normal at the projection (see
/// `surface_foot`; it was a finite-difference normal, PK-776).
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
    let (_, idx, uv) = best?;
    crate::construct::surface_foot::unit_normal(&candidates[idx], uv)
}

/// The target surfaces for one rim side: its whole-row support when
/// one candidate carries the entire row, else every strip (per-point
/// nearest).
pub(super) fn side_targets<'a>(
    candidates: &'a [SurfaceGeometry],
    strips: &'a [SurfaceGeometry],
    row0: &[Point3; 10],
) -> &'a [SurfaceGeometry] {
    match side_support(candidates, strips.len(), row0) {
        Some(i) => std::slice::from_ref(&candidates[i]),
        None => strips,
    }
}

/// Run the rim cross-jet reconciliation over the committed fan.
pub(super) fn reconcile_rim_jets(
    patches: &mut [GregoryPatch],
    candidates: &[SurfaceGeometry],
    strip_count: usize,
    border_loop: &[Vec<Point3>],
    rim_off_loop_allowance: f64,
) {
    // Tangency targets: each side's WHOLE-ROW support (the candidate
    // carrying the entire row — `side_targets`), else per point from
    // the STRIPS (the first `strip_count` candidates). Per-point
    // nearest-candidate lookup is wrong both ways: near the pulled-back
    // corners a LATERAL PLANE can sit nearer a strip-borne rim in
    // space, and a saddle's bridge rim lies ON its reflex support plane
    // while a strip is within the rim band at every sample. Targeting
    // the bridge against the strip tilted the saddle's rims and spokes
    // to 4e-3..1e-2 (PK-576; whole-row attribution: spokes 1e-11,
    // cap-support 1e-10).
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
        // Targets from the side's whole-row support (a bridge side
        // on the reflex plane takes the plane), else per-point from
        // the nearest strip.
        let targets = side_targets(candidates, strips, &row0);
        let mut field: Vec<[f64; 3]> = Vec::with_capacity(M);
        let mut ok = true;
        for &t in &ts {
            let p0 = bezier9_eval(&row0, t);
            let p1 = bezier9_eval(&row1, t);
            let Some(n) = support_normal(targets, &p0) else {
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
            start_cluster: (start_dual.is_none() || spoke_dual_end(&sides, pi, side, false))
                .then(|| cluster_of(bezier9_eval(&row0, 0.0))),
            end_cluster: (end_dual.is_none() || spoke_dual_end(&sides, pi, side, true))
                .then(|| cluster_of(bezier9_eval(&row0, 1.0))),
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
    // Affine cluster fields (PK-575, `setback_cluster_field`), as in
    // the curvature pass: each band cell moves by the field at ITS OWN
    // position (row 1 by δ(x₁), row 2 by 2·δ(x₂)), which keeps every
    // linear seam relation exact.
    let mut fields_at = ClusterFields::new(centers);
    for js in &jet_sides {
        let (Some(r1), Some(r2)) = (
            row_of(patches, js.pi, js.side, 1),
            row_of(patches, js.pi, js.side, 2),
        ) else {
            continue;
        };
        for (cluster, ks, t) in [
            (js.start_cluster, 0..=2usize, 0.0),
            (js.end_cluster, 7..=9usize, 1.0),
        ] {
            if let Some(c) = cluster {
                fields_at.include_side(c, js.pi, js.side, t > 0.5);
                for k in ks {
                    fields_at.include(c, r1[k]);
                    fields_at.include(c, r2[k]);
                }
            }
        }
    }
    fields_at.classify(patches);
    let mut next_unknown = CLUSTER_FIELD_UNKNOWNS * n_clusters;
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
    let mut lsq = VecLsq::new(unknowns);
    for js in &jet_sides {
        let Some(row1_of) = row_of(patches, js.pi, js.side, 1) else {
            continue;
        };
        for (m, target) in js.field.iter().enumerate() {
            let t = ts[m];
            let mut entries: Vec<(usize, f64)> = Vec::with_capacity(10);
            if js.start_dual.is_some() {
                let (idx1, idx2) = js.dual_indices(false);
                entries.push((idx1, bern9(1, t)));
                entries.push((idx2, bern9(2, t)));
            }
            if let Some(cluster) = js.start_cluster {
                for (k, &x1) in row1_of.iter().enumerate().take(3) {
                    let q = fields_at.coords(cluster, x1);
                    for (i, x) in q.iter().enumerate() {
                        entries.push((CLUSTER_FIELD_UNKNOWNS * cluster + i, bern9(k, t) * x));
                    }
                }
            }
            if js.end_dual.is_some() {
                let (idx1, idx2) = js.dual_indices(true);
                entries.push((idx1, bern9(8, t)));
                entries.push((idx2, bern9(7, t)));
            }
            if let Some(cluster) = js.end_cluster {
                for (k, &x1) in row1_of.iter().enumerate().take(9 + 1).skip(7) {
                    let q = fields_at.coords(cluster, x1);
                    for (i, x) in q.iter().enumerate() {
                        entries.push((CLUSTER_FIELD_UNKNOWNS * cluster + i, bern9(k, t) * x));
                    }
                }
            }
            for q in 0..4usize {
                entries.push((js.interior_base + q, bern9(3 + q, t)));
            }
            lsq.add_vec_row(&entries, Vec3::new(target[0], target[1], target[2]));
        }
    }
    let Some(solution) = lsq.solve() else {
        return;
    };
    for js in &jet_sides {
        let JetSide { pi, side, .. } = *js;
        let (Some(row1_pre), Some(row2_pre)) =
            (row_of(patches, pi, side, 1), row_of(patches, pi, side, 2))
        else {
            continue;
        };
        for k in 0..10usize {
            let dual = match k {
                0..=2 => js.start_dual,
                7..=9 => js.end_dual,
                _ => None,
            };
            let cluster_here = match k {
                0..=2 => js.start_cluster,
                7..=9 => js.end_cluster,
                _ => None,
            };
            // The cluster field at each row's own cell (a spoke
            // cluster's field is a translation).
            let field1 = cluster_here.map(|c| fields_at.delta(&solution, c, row1_pre[k]));
            let field2 = cluster_here.map(|c| fields_at.delta(&solution, c, row2_pre[k]) * 2.0);
            // Dual-end routing: the end cell itself is the OTHER
            // boundary's d0/d1 row. At a setback corner it stays fixed;
            // at a spoke junction it moves with the junction's shared
            // translation. One and two steps out, the side's own d1/d2
            // corrections go into its own claims, and the shared
            // translation moves BOTH claim slots (the spoke keeps its
            // cross vectors).
            if let Some(corner) = dual {
                let at_end_pos = if k <= 2 { k } else { 9 - k };
                let shared1 = field1.unwrap_or(Vec3::zero());
                let shared2 = field2.unwrap_or(Vec3::zero());
                if at_end_pos == 0 {
                    if cluster_here.is_some() {
                        let (i1, j1) = coord_for(side, 1, k);
                        let (i2, j2) = coord_for(side, 2, k);
                        if let Some(grid) = patches[pi].binonic_grid_mut() {
                            let g1 = grid[i1][j1];
                            grid[i1][j1] =
                                Point3::new(g1.x + shared1.x, g1.y + shared1.y, g1.z + shared1.z);
                            let g2 = grid[i2][j2];
                            grid[i2][j2] =
                                Point3::new(g2.x + shared2.x, g2.y + shared2.y, g2.z + shared2.z);
                        }
                    }
                    continue;
                }
                let (idx1, idx2) = js.dual_indices(k >= 7);
                let own = if at_end_pos == 1 {
                    solution[idx1]
                } else {
                    solution[idx2]
                };
                let d1_cell = dual_end_d1_cells(side)[at_end_pos - 1];
                let d2_cell = dual_end_cells(side)[at_end_pos - 1];
                add_claim_delta(patches, pi, side, corner, d1_cell, own + shared1);
                add_claim_delta(patches, pi, side, corner, d2_cell, own * 2.0 + shared2);
                if cluster_here.is_some() {
                    let other = fan_side(side, k >= 7);
                    add_claim_delta(patches, pi, other, corner, d1_cell, shared1);
                    add_claim_delta(patches, pi, other, corner, d2_cell, shared2);
                }
                continue;
            }
            let delta = match field1 {
                Some(f) => f,
                None => solution[js.interior_base + (k - 3)],
            };
            // Row 2 under an affine cluster field: displaced by the
            // field at its OWN position (2·δ(x₂)), so seam relations
            // along row 2 stay exact.
            let delta2 = field2.unwrap_or(delta * 2.0);
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
            grid[i2][j2] = Point3::new(g2.x + delta2.x, g2.y + delta2.y, g2.z + delta2.z);
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
                    if let Some(n) = support_normal(side_targets(candidates, strips, &row0), &p0) {
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
