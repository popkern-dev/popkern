//! Face simplification: merge contiguous co-surfacial faces.
//!
//! After redundant edge removal, adjacent faces that share the same
//! analytical surface can be merged into a single face. This reduces
//! face count and simplifies downstream operations.
//!
//! Only merges faces on analytical surfaces (Plane, Cylinder, Cone,
//! Sphere, Torus) where co-surfacial detection is exact. NURBS faces
//! are not merged because approximate comparison is unreliable.

use popkern_core::geometry::intersection::coincident_analytic::are_surfaces_coincident_analytic;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Result of face simplification analysis.
#[derive(Debug, Clone, Default)]
pub struct FaceSimplifyReport {
    /// Groups of faces that could be merged (each group shares a surface).
    pub merge_groups: Vec<Vec<FaceHandle>>,
    /// Total faces analyzed.
    pub faces_analyzed: usize,
}

impl FaceSimplifyReport {
    /// Number of faces that could be eliminated by merging.
    pub fn potential_reduction(&self) -> usize {
        self.merge_groups
            .iter()
            .map(|g| if g.len() > 1 { g.len() - 1 } else { 0 })
            .sum()
    }
}

/// Find groups of contiguous faces that share the same analytic surface.
///
/// Returns groups of face handles where each group could be merged into
/// a single face. Only analytic surfaces are considered.
///
/// Note: actual merging (Kill-Edge-Face Euler operations) is not performed
/// here — this function only identifies the candidates.
pub fn find_merge_candidates(brep: &BRep, tolerance: &Tolerance) -> FaceSimplifyReport {
    let mut report = FaceSimplifyReport {
        faces_analyzed: brep.faces.len(),
        ..Default::default()
    };

    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    let mut visited = std::collections::HashSet::new();
    // Each face's edge-neighbours (the other side of each of its edges),
    // in edge order.
    let mut adjacency: std::collections::HashMap<FaceHandle, Vec<FaceHandle>> =
        std::collections::HashMap::new();
    for (_eh, edge) in &brep.edges {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let (Some(la), Some(lb)) = (brep.loops.get(he.loop_), brep.loops.get(twin.loop_)) else {
            continue;
        };
        adjacency.entry(la.face).or_default().push(lb.face);
        if lb.face != la.face {
            adjacency.entry(lb.face).or_default().push(la.face);
        }
    }

    for &fh in &face_handles {
        if visited.contains(&fh) {
            continue;
        }

        let surf = match brep.faces.get(fh).and_then(|f| f.surface.as_ref()) {
            Some(s) => s,
            None => continue,
        };

        // Skip NURBS surfaces
        if matches!(
            surf.geometry,
            popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_)
        ) {
            continue;
        }

        // Find all adjacent faces with the same surface
        let mut group = vec![fh];
        visited.insert(fh);

        // Chart-span budget for DOUBLY-periodic surfaces (torus):
        // merging every co-surfacial patch of a torus produced one
        // full-wrap mega-face whose loop cannot close in any UV
        // branch (measured: 31-use loop, UV gap 0.67, STEP demoted
        // to OPEN_SHELL, OCC NULL shape). A merge is admitted only
        // while the group's combined chart cell stays strictly
        // inside one period per axis — the same cell-resolvability
        // contract as clipping and classification. Singly-periodic
        // surfaces (cylinder walls) keep the long-tuned unguarded
        // behavior.
        let (u_period, v_period) = surf.geometry.periodicity();
        let doubly_periodic = u_period.is_some() && v_period.is_some();
        let mut group_span: Option<(f64, f64, f64, f64)> = None;
        if doubly_periodic {
            match crate::clip::chart_cell::face_chart_cell(brep, fh, tolerance) {
                Some(cell) => {
                    group_span = Some((
                        cell.u_range.0,
                        cell.u_range.1,
                        cell.v_range.0,
                        cell.v_range.1,
                    ));
                }
                // A non-cell-resolvable periodic face (already
                // winding) must not seed further growth.
                None => continue,
            }
        }

        let mut queue = std::collections::VecDeque::new();
        queue.push_back(fh);

        while let Some(current) = queue.pop_front() {
            // Adjacent faces via shared edges, in edge order (the
            // adjacency is built once: scanning every edge per visited
            // face was quadratic in body size — PK-804).
            for &neighbor in adjacency.get(&current).map(Vec::as_slice).unwrap_or(&[]) {
                if visited.contains(&neighbor) {
                    continue;
                }
                let neighbor_surf = match brep.faces.get(neighbor).and_then(|f| f.surface.as_ref())
                {
                    Some(s) => s,
                    None => continue,
                };

                // Same surface, same orientation
                if surf.orientation == neighbor_surf.orientation
                    && are_surfaces_coincident_analytic(
                        &surf.geometry,
                        &neighbor_surf.geometry,
                        tolerance,
                    )
                {
                    if doubly_periodic {
                        let Some((lo_u, hi_u, lo_v, hi_v)) = group_span else {
                            continue;
                        };
                        let Some(ncell) =
                            crate::clip::chart_cell::face_chart_cell(brep, neighbor, tolerance)
                        else {
                            continue;
                        };
                        // Align the neighbour's cell to the group's
                        // branch (shift by whole periods toward the
                        // group centre), then check the combined span.
                        let align = |lo: f64, hi: f64, glo: f64, ghi: f64, period: f64| {
                            let centre = 0.5 * (lo + hi);
                            let gcentre = 0.5 * (glo + ghi);
                            let shift = ((gcentre - centre) / period).round() * period;
                            (lo + shift, hi + shift)
                        };
                        let (p_u, p_v) = match (u_period, v_period) {
                            (Some(pu), Some(pv)) => (pu, pv),
                            _ => continue,
                        };
                        let (nlo_u, nhi_u) =
                            align(ncell.u_range.0, ncell.u_range.1, lo_u, hi_u, p_u);
                        let (nlo_v, nhi_v) =
                            align(ncell.v_range.0, ncell.v_range.1, lo_v, hi_v, p_v);
                        let clo_u = lo_u.min(nlo_u);
                        let chi_u = hi_u.max(nhi_u);
                        let clo_v = lo_v.min(nlo_v);
                        let chi_v = hi_v.max(nhi_v);
                        if chi_u - clo_u >= p_u - 1.0e-9 || chi_v - clo_v >= p_v - 1.0e-9 {
                            continue;
                        }
                        group_span = Some((clo_u, chi_u, clo_v, chi_v));
                    }
                    group.push(neighbor);
                    visited.insert(neighbor);
                    queue.push_back(neighbor);
                }
            }
        }

        if group.len() > 1 {
            report.merge_groups.push(group);
        }
    }

    report
}

/// Apply face simplification: merge co-surfacial adjacent faces via KEF.
///
/// Iterates until no more merges are possible. After each individual
/// KEF, re-runs `find_merge_candidates` so the next iteration sees
/// the post-merge topology — this avoids the stale-handle bug where
/// `shared_edges` was computed once but the BFS group used stale
/// face references after the first merge.
///
/// Validates after each KEF; on validation failure, rolls back the
/// merge by restoring from a snapshot and stops the pass. Returns
/// the number of faces eliminated.
pub fn apply_face_simplification(brep: &mut BRep, tolerance: &Tolerance) -> usize {
    let mut eliminated = 0;
    // Hard cap so a pathological case can't run forever. Each merge
    // strictly decreases the face count, so the actual cap is the
    // initial face count anyway.
    let max_iterations = brep.faces.len() + 1;

    for _ in 0..max_iterations {
        let report = find_merge_candidates(brep, tolerance);
        if report.merge_groups.is_empty() {
            break;
        }

        let mut merged_this_iter = false;
        for group in &report.merge_groups {
            if group.len() < 2 {
                continue;
            }
            let group_set: std::collections::HashSet<FaceHandle> = group.iter().copied().collect();

            // Find one shared edge between any two distinct faces in
            // the group. Picking just one — and re-running find_merge_candidates
            // afterwards — keeps the loop self-consistent with the
            // mutating topology.
            //
            // KEF is only sound when the pair shares exactly ONE edge:
            // dissolving one edge of a multi-edge shared boundary
            // leaves the surviving shared edges traversed twice by the
            // merged loop — a zero-width slit that passes Euler and
            // manifold validation but cannot be triangulated
            // (measured on a flush coplanar union: an arm wall and the
            // gusset cap shared both the slope record and the corner
            // column; the merged face lost 70% of the swept volume to
            // mesh failures). Pairs sharing several edges stay split —
            // a valid, meshable result.
            let mut pair_shares: Vec<(
                (FaceHandle, FaceHandle),
                popkern_core::topology::HalfEdgeHandle,
                usize,
            )> = Vec::new();
            for (_, edge) in brep.edges.iter() {
                let Some(he) = brep.half_edges.get(edge.half_edge) else {
                    continue;
                };
                let Some(twin) = brep.half_edges.get(he.twin) else {
                    continue;
                };
                let (Some(la), Some(lb)) = (brep.loops.get(he.loop_), brep.loops.get(twin.loop_))
                else {
                    continue;
                };
                let (fa, fb) = (la.face, lb.face);
                if fa == fb || !group_set.contains(&fa) || !group_set.contains(&fb) {
                    continue;
                }
                let key = if fa < fb { (fa, fb) } else { (fb, fa) };
                match pair_shares.iter_mut().find(|(k, _, _)| *k == key) {
                    Some((_, _, count)) => *count += 1,
                    None => pair_shares.push((key, edge.half_edge, 1)),
                }
            }
            let shared = pair_shares
                .iter()
                .find(|(_, _, count)| *count == 1)
                .map(|(_, he_h, _)| *he_h);

            let Some(he_h) = shared else {
                continue;
            };

            // Snapshot for rollback on validation regression.
            let snapshot = brep.clone();
            if brep.kef(he_h).is_ok() && brep.validate().is_ok() {
                eliminated += 1;
                merged_this_iter = true;
                // Restart the outer loop so find_merge_candidates
                // sees the post-merge topology.
                break;
            } else {
                *brep = snapshot;
            }
        }

        if !merged_this_iter {
            break;
        }
    }

    eliminated
}
