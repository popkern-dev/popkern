//! Rim transverse-curvature reconciliation (PK-515): after the G0
//! rim landing and the G1 cross-jet reconcile, move each rim side's
//! DEPTH-2 row so the cap's transverse normal curvature along the
//! rim follows the strip's, measured pointwise on both sides.
//!
//! Why this exists: the border derivation encodes the strips' cross
//! curvature into the bi-6 fan through ONE shared Hermite block per
//! corner, and two neighbouring fillet strips are different surfaces
//! whose cross fields agree only to second order at the corner — the
//! block cannot carry both (measured floor 3.4e-4 in row units on the
//! G2 pentagon flagship at ANY objective weight, ≈ 1 in transverse
//! curvature against the 5e-2 commercial gate; re-weighting only
//! relocated the error into jet tilt and corner twist). The committed
//! pieces are bi-9: their ten-CP depth-2 rows have the freedom the
//! bi-6 block lacks, and a depth-2-only move leaves rows 0 and 1 —
//! G0, tangency, and the rim twist — bit-identical.
//!
//! Cross-seam discipline (as the jet reconcile): the three CPs
//! nearest each end take ONE shared per-corner delta (translation-
//! invariant G rows); the interior four CPs fit the measured field.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};

use super::bezier9_eval;
use super::dense_spd_solve::VecLsq;
use super::setback_cluster_field::{
    CLUSTER_FIELD_UNKNOWNS, ClusterFields, fan_side, spoke_dual_end,
};
use super::setback_corner_dual_jet::{
    active_corner_at_end, add_claim_delta, coord_for, dual_end_cells, rim_sides, row_for_side,
};
use super::setback_rim_stations::{Station, side_stations};

const CORNER_CLUSTER: f64 = 1.0e-3;

fn row_of(patches: &[GregoryPatch], patch: usize, side: usize, d: usize) -> Option<[Point3; 10]> {
    row_for_side(patches, patch, side, d)
}

/// Run the rim transverse-curvature reconciliation over the fan.
pub(super) fn reconcile_rim_curvature(
    patches: &mut [GregoryPatch],
    candidates: &[SurfaceGeometry],
    strip_count: usize,
    border_loop: &[Vec<Point3>],
    rim_off_loop_allowance: f64,
) {
    let strips = &candidates[..strip_count.min(candidates.len())];
    let sides = rim_sides(patches, candidates, border_loop, rim_off_loop_allowance);
    if sides.is_empty() {
        return;
    }
    // Per-side station fields; a side whose measurement declines at
    // any station is skipped whole (a partial correction would kink
    // the row against its clamped ends).
    let mut fields: Vec<((usize, usize), Vec<Station>)> = Vec::new();
    for &(pi, side) in &sides {
        let Some(row0) = row_for_side(patches, pi, side, 0) else {
            continue;
        };
        let targets = super::setback_rim_jet::side_targets(candidates, strips, &row0);
        if let Some(stations) = side_stations(patches, targets, pi, side) {
            fields.push(((pi, side), stations));
        }
    }
    // ONE joint least-squares over the whole ring: unknowns are the
    // shared per-corner deltas (one Vec3 per corner cluster - the
    // three CPs nearest each end of every incident side take it) and
    // each side's four interior CPs. Solving jointly picks the corner
    // constant that best serves the clamp-region FIELD of both
    // incident sides - an endpoint-pinned constant left the measured
    // residual peaked at 0.68 exactly inside the clamp region.
    // A side end is either clamped to a shared per-cluster constant
    // (the fallback semantics) or, at an activated dual-jet corner,
    // free through the side's OWN claims of the two cells nearest the
    // corner (the end CP itself — the other rim's G0 curve CP — stays
    // fixed). Dual ends never join a cluster: the whole point of the
    // substrate is that the two incident sides' corrections are
    // independent there.
    struct SideLayout {
        pi: usize,
        side: usize,
        /// Shared-constant cluster per end, `None` at a dual end.
        start_cluster: Option<usize>,
        end_cluster: Option<usize>,
        /// Activated dual-jet corner per end.
        start_dual: Option<usize>,
        end_dual: Option<usize>,
        interior_base: usize,
        /// Unknown base for the dual-end claim deltas: two unknowns
        /// per dual end (1-step cell then 2-step cell), start end
        /// first when both are dual.
        dual_base: usize,
    }
    impl SideLayout {
        /// Unknown indices of a dual end's (1-step, 2-step) cells.
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
    let mut layouts: Vec<SideLayout> = Vec::new();
    for ((pi, side), _) in &fields {
        let Some(row0) = row_of(patches, *pi, *side, 0) else {
            continue;
        };
        let start_dual = active_corner_at_end(patches, *pi, *side, false);
        let end_dual = active_corner_at_end(patches, *pi, *side, true);
        // A cluster joins every end off a SETBACK corner: plain ends
        // (affine field) and spoke-junction dual ends (translation of
        // the junction neighbourhood, per-side claims on top).
        let clustered = |at_k9: bool, dual: Option<usize>| {
            dual.is_none() || spoke_dual_end(&sides, *pi, *side, at_k9)
        };
        layouts.push(SideLayout {
            pi: *pi,
            side: *side,
            start_cluster: clustered(false, start_dual)
                .then(|| cluster_of(bezier9_eval(&row0, 0.0))),
            end_cluster: clustered(true, end_dual).then(|| cluster_of(bezier9_eval(&row0, 1.0))),
            start_dual,
            end_dual,
            interior_base: 0,
            dual_base: 0,
        });
    }
    let n_clusters = centers.len();
    // Affine cluster fields (PK-575, `setback_cluster_field`): each
    // non-dual corner cluster's band cells move by one shared affine
    // displacement field (seam-safe; no freeze needed).
    let mut fields_at = ClusterFields::new(centers);
    let rows2: Vec<Option<[Point3; 10]>> = layouts
        .iter()
        .map(|l| row_of(patches, l.pi, l.side, 2))
        .collect();
    for (layout, row2) in layouts.iter().zip(&rows2) {
        let Some(row2) = row2 else {
            continue;
        };
        for (cluster, ks, t) in [
            (layout.start_cluster, 0..=2usize, 0.0),
            (layout.end_cluster, 7..=9usize, 1.0),
        ] {
            if let Some(cluster) = cluster {
                fields_at.include_side(cluster, layout.pi, layout.side, t > 0.5);
                for k in ks {
                    fields_at.include(cluster, row2[k]);
                }
            }
        }
    }
    fields_at.classify(patches);
    let band_coords = |cluster: usize, li: usize, k: usize| -> [f64; CLUSTER_FIELD_UNKNOWNS] {
        match &rows2[li] {
            Some(row2) => fields_at.coords(cluster, row2[k]),
            None => [1.0, 0.0, 0.0, 0.0],
        }
    };
    let mut next_unknown = CLUSTER_FIELD_UNKNOWNS * n_clusters;
    for layout in layouts.iter_mut() {
        layout.interior_base = next_unknown;
        next_unknown += 4;
        layout.dual_base = next_unknown;
        next_unknown +=
            2 * (layout.start_dual.is_some() as usize + layout.end_dual.is_some() as usize);
    }
    let unknowns = next_unknown;
    let bern9 = |k: usize, t: f64| -> f64 {
        const C: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];
        C[k] * t.powi(k as i32) * (1.0 - t).powi((9 - k) as i32)
    };
    let mut lsq = VecLsq::new(unknowns);
    for (li, (layout, ((_, _), stations))) in layouts.iter().zip(&fields).enumerate() {
        for station in stations {
            let target = station.normal * (station.delta_n_bvv / 72.0);
            // Row: per-end unknowns (dual claims, or clamp mass on a
            // shared cluster constant) + the four interiors.
            let mut entries: Vec<(usize, f64)> = Vec::with_capacity(8);
            if layout.start_dual.is_some() {
                let (idx1, idx2) = layout.dual_indices(false);
                entries.push((idx1, bern9(1, station.t)));
                entries.push((idx2, bern9(2, station.t)));
            }
            if let Some(cluster) = layout.start_cluster {
                for k in 0..=2usize {
                    let w = bern9(k, station.t);
                    for (q, x) in band_coords(cluster, li, k).iter().enumerate() {
                        entries.push((CLUSTER_FIELD_UNKNOWNS * cluster + q, w * x));
                    }
                }
            }
            if layout.end_dual.is_some() {
                let (idx1, idx2) = layout.dual_indices(true);
                entries.push((idx1, bern9(8, station.t)));
                entries.push((idx2, bern9(7, station.t)));
            }
            if let Some(cluster) = layout.end_cluster {
                for k in 7..=9usize {
                    let w = bern9(k, station.t);
                    for (q, x) in band_coords(cluster, li, k).iter().enumerate() {
                        entries.push((CLUSTER_FIELD_UNKNOWNS * cluster + q, w * x));
                    }
                }
            }
            for k in 3..=6_usize {
                entries.push((layout.interior_base + (k - 3), bern9(k, station.t)));
            }
            lsq.add_vec_row(&entries, target);
        }
    }
    let Some(solution) = lsq.solve() else {
        return;
    };
    // Protected cells: any registered rim side's d = 0 (G0) and
    // d = 1 (tangency/twist) rows are untouchable. The CORNER piece
    // carries TWO rim sides meeting at `c`, and each side's d = 2 row
    // crosses the OTHER side's shallow rows there — writing those
    // cells moved a genuine rim boundary row and opened the
    // cap↔strip seam (measured: pyramid5_mixed cap_outer g0 6.1e-6
    // against its 1e-9 bar). Cells claimed by two sides' d = 2 rows
    // (the (7,7)-class corner cell) take the MEAN of the claims.
    let mut protected: std::collections::BTreeSet<(usize, usize, usize)> =
        std::collections::BTreeSet::new();
    for layout in &layouts {
        for d in 0..=1_usize {
            for k in 0..10_usize {
                let (i, j) = coord_for(layout.side, d, k);
                protected.insert((layout.pi, i, j));
            }
        }
    }
    // Per-end delta of the d2-row CP at position `k`, or `None` for
    // a dual end's cells (those route through the side's own claims)
    // and the fixed end CP at a dual corner.
    let cluster_delta = |c: usize, li: usize, k: usize| -> Vec3 {
        rows2[li].map_or(Vec3::zero(), |row2| fields_at.delta(&solution, c, row2[k]))
    };
    let end_delta = |li: usize, layout: &SideLayout, k: usize| -> Option<Vec3> {
        match k {
            0 if layout.start_dual.is_some() => {
                layout.start_cluster.map(|c| cluster_delta(c, li, k))
            }
            0..=2 if layout.start_dual.is_some() => None,
            0..=2 => Some(
                layout
                    .start_cluster
                    .map_or(Vec3::zero(), |c| cluster_delta(c, li, k)),
            ),
            9 if layout.end_dual.is_some() => layout.end_cluster.map(|c| cluster_delta(c, li, k)),
            7..=9 if layout.end_dual.is_some() => None,
            7..=9 => Some(
                layout
                    .end_cluster
                    .map_or(Vec3::zero(), |c| cluster_delta(c, li, k)),
            ),
            _ => Some(solution[layout.interior_base + (k - 3)]),
        }
    };
    let mut writes: std::collections::BTreeMap<(usize, usize, usize), (Vec3, usize)> =
        std::collections::BTreeMap::new();
    for (li, layout) in layouts.iter().enumerate() {
        for k in 0..10_usize {
            let Some(delta) = end_delta(li, layout, k) else {
                continue;
            };
            let (i2, j2) = coord_for(layout.side, 2, k);
            if protected.contains(&(layout.pi, i2, j2)) {
                continue;
            }
            let entry = writes
                .entry((layout.pi, i2, j2))
                .or_insert((Vec3::zero(), 0));
            entry.0 = entry.0 + delta;
            entry.1 += 1;
        }
    }
    for ((pi, i2, j2), (sum, count)) in writes {
        let Some(grid) = patches[pi].binonic_grid_mut() else {
            continue;
        };
        let delta = sum * (1.0 / count.max(1) as f64);
        let g = grid[i2][j2];
        grid[i2][j2] = Point3::new(g.x + delta.x, g.y + delta.y, g.z + delta.z);
    }
    // Dual-end claim writes: each side's correction lands in its OWN
    // claims of the two cells nearest the corner — no shared
    // constant, no doubly-claimed mean, no protection conflict (the
    // other side's G0/G1 rows read their own claims and the plain
    // grid, which these writes never touch).
    let mut worst_claim_delta = 0.0_f64;
    for layout in &layouts {
        for (end, dual) in [(false, layout.start_dual), (true, layout.end_dual)] {
            let Some(corner) = dual else {
                continue;
            };
            let (idx1, idx2) = layout.dual_indices(end);
            let cells = dual_end_cells(layout.side);
            let cluster = if end {
                layout.end_cluster
            } else {
                layout.start_cluster
            };
            let li = layouts
                .iter()
                .position(|l| l.pi == layout.pi && l.side == layout.side)
                .unwrap_or(0);
            for (step, (cell, idx)) in cells.iter().zip([idx1, idx2]).enumerate() {
                // Spoke junction: the junction neighbourhood's shared
                // translation moves BOTH claim slots (the spoke keeps
                // its cross vectors); the side's own correction lands
                // in its slot only.
                let k = if end { 8 - step } else { 1 + step };
                let shared = cluster.map_or(Vec3::zero(), |c| cluster_delta(c, li, k));
                add_claim_delta(
                    patches,
                    layout.pi,
                    layout.side,
                    corner,
                    *cell,
                    solution[idx] + shared,
                );
                if cluster.is_some() {
                    let other = fan_side(layout.side, end);
                    add_claim_delta(patches, layout.pi, other, corner, *cell, shared);
                }
                worst_claim_delta = worst_claim_delta.max(solution[idx].length());
            }
        }
    }
    if popkern_core::trace::trace_enabled("setback") {
        for (li, layout) in layouts.iter().enumerate() {
            let d_max = (0..10_usize)
                .filter_map(|k| end_delta(li, layout, k).map(|d| d.length()))
                .fold(0.0_f64, f64::max);
            let frozen_label = |cluster: Option<usize>, dual: Option<usize>| match (cluster, dual) {
                (_, Some(c)) => format!("dual{c}"),
                (Some(c), _) => format!("field{c}"),
                _ => "-".to_string(),
            };
            eprintln!(
                "  [rim-curvature] applied patch {} side {} |delta|max {:.3e} ends=({},{})",
                layout.pi,
                layout.side,
                d_max,
                frozen_label(layout.start_cluster, layout.start_dual),
                frozen_label(layout.end_cluster, layout.end_dual)
            );
        }
        let worst = fields
            .iter()
            .flat_map(|(_, s)| s.iter())
            .map(|s| (s.delta_n_bvv / 72.0).abs())
            .fold(0.0_f64, f64::max);
        let dual_ends: usize = layouts
            .iter()
            .map(|l| l.start_dual.is_some() as usize + l.end_dual.is_some() as usize)
            .sum();
        eprintln!(
            "  [rim-curvature] {} sides, {} corner clusters, {} dual ends (worst claim delta {worst_claim_delta:.3e}): worst |drow2| target {worst:.3e} applied (joint solve)",
            layouts.len(),
            n_clusters,
            dual_ends
        );
    }
}
