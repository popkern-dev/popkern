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
    find_merge_candidates_in(brep, tolerance, &crate::serve_region::ServeRegion::Whole)
}

/// [`find_merge_candidates`] among the faces `region` admits: groups
/// neither start nor grow outside it.
fn find_merge_candidates_in(
    brep: &BRep,
    tolerance: &Tolerance,
    region: &crate::serve_region::ServeRegion,
) -> FaceSimplifyReport {
    let mut report = FaceSimplifyReport {
        faces_analyzed: brep.faces.len(),
        ..Default::default()
    };

    let face_handles: Vec<FaceHandle> = brep
        .faces
        .keys()
        .filter(|&f| region.admits_face(f))
        .collect();
    // Dense on the body's face handles: the pass runs on every
    // boolean's whole result (PK-815).
    let mut visited: slotmap::SecondaryMap<FaceHandle, ()> =
        slotmap::SecondaryMap::with_capacity(brep.faces.len());
    // Each face's edge-neighbours (the other side of each of its edges),
    // in edge order.
    let mut adjacency: slotmap::SecondaryMap<FaceHandle, Vec<FaceHandle>> =
        slotmap::SecondaryMap::with_capacity(brep.faces.len());
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
        if !(region.admits_face(la.face) && region.admits_face(lb.face)) {
            continue;
        }
        if let Some(e) = adjacency.entry(la.face) {
            e.or_default().push(lb.face);
        }
        if lb.face != la.face
            && let Some(e) = adjacency.entry(lb.face)
        {
            e.or_default().push(la.face);
        }
    }

    for &fh in &face_handles {
        if visited.contains_key(fh) {
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
        visited.insert(fh, ());

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
        //
        // A face that wraps ONE axis through its own seam (a torus
        // band round the tube, or either half of a cut-domain walk
        // split at the u-seam, PK-922) has no chart cell; its span on
        // that axis is the whole period and the other axis is bounded
        // by its pcurve window. Such faces group with neighbours
        // wrapping the same axis while the bounded axis stays strictly
        // inside one period.
        let (u_period, v_period) = surf.geometry.periodicity();
        let periods = u_period.zip(v_period);
        let mut group_span: Option<ChartSpan> = None;
        if let Some(periods) = periods {
            match doubly_periodic_span(brep, fh, tolerance, periods) {
                Some(span) => group_span = Some(span),
                // A face wrapping both axes, or one with no readable
                // span, must not seed further growth.
                None => continue,
            }
        }

        let mut queue = std::collections::VecDeque::new();
        queue.push_back(fh);

        while let Some(current) = queue.pop_front() {
            // Adjacent faces via shared edges, in edge order (the
            // adjacency is built once: scanning every edge per visited
            // face was quadratic in body size — PK-804).
            for &neighbor in adjacency.get(current).map(Vec::as_slice).unwrap_or(&[]) {
                if visited.contains_key(neighbor) {
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
                    if let Some(periods) = periods {
                        let Some(span) = group_span else {
                            continue;
                        };
                        let Some(nspan) = doubly_periodic_span(brep, neighbor, tolerance, periods)
                        else {
                            continue;
                        };
                        let Some(combined) = combine_spans(span, nspan, periods) else {
                            continue;
                        };
                        group_span = Some(combined);
                    }
                    group.push(neighbor);
                    visited.insert(neighbor, ());
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

/// A doubly-periodic face's chart span per axis, `None` on an axis the
/// face wraps through its own seam.
type ChartSpan = [Option<(f64, f64)>; 2];

/// The chart span of `face` on a surface with periods `(pu, pv)`: its
/// chart cell when it has one; else, when it wraps exactly one axis by
/// its own seam (`chart_seam_axes`), its pcurve window on the other,
/// which must lie strictly inside one period. `None` for a face
/// wrapping both axes or with no readable span.
fn doubly_periodic_span(
    brep: &BRep,
    face: FaceHandle,
    tolerance: &Tolerance,
    (pu, pv): (f64, f64),
) -> Option<ChartSpan> {
    if let Some(cell) = crate::clip::chart_cell::face_chart_cell(brep, face, tolerance) {
        return Some([Some(cell.u_range), Some(cell.v_range)]);
    }
    let seams = crate::merge::walk::chart_seam_axes::face_seam_axes(brep, face, (pu, pv));
    let ((u0, u1), (v0, v1)) = brep.face_trim_uv_window(face)?;
    match (seams.u.is_some(), seams.v.is_some()) {
        (true, false) if v1 - v0 < pv - 1.0e-9 => Some([None, Some((v0, v1))]),
        (false, true) if u1 - u0 < pu - 1.0e-9 => Some([Some((u0, u1)), None]),
        _ => None,
    }
}

/// The span of two adjacent faces taken as one, or `None` when they
/// wrap different axes or the merged span would reach a full period on
/// an axis neither wraps. The neighbour's window is shifted by whole
/// periods toward the group's before combining.
fn combine_spans(group: ChartSpan, neighbor: ChartSpan, periods: (f64, f64)) -> Option<ChartSpan> {
    let mut out = group;
    for (k, period) in [periods.0, periods.1].into_iter().enumerate() {
        out[k] = match (group[k], neighbor[k]) {
            (None, None) => None,
            (Some((glo, ghi)), Some((lo, hi))) => {
                let shift = ((0.5 * (glo + ghi) - 0.5 * (lo + hi)) / period).round() * period;
                let (clo, chi) = (glo.min(lo + shift), ghi.max(hi + shift));
                if chi - clo >= period - 1.0e-9 {
                    return None;
                }
                Some((clo, chi))
            }
            _ => return None,
        };
    }
    Some(out)
}

/// After a KEF, make the survivor's loops chart-continuous when its
/// surface is doubly periodic: each half kept its pcurves on its own
/// branch, so the merged loop jumped a period at the old shared edge.
/// Shifts are whole periods, so the pcurves stay exact. Other surfaces
/// keep their pcurves.
fn rechart_doubly_periodic(brep: &mut BRep, surviving: popkern_core::topology::LoopHandle) -> bool {
    let Some(face) = brep.loops.get(surviving).map(|l| l.face) else {
        return false;
    };
    let doubly = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .is_some_and(|s| matches!(s.geometry.periodicity(), (Some(_), Some(_))));
    if doubly {
        brep.normalize_loop_chart_continuity(&[face]);
    }
    true
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
    apply_face_simplification_in(brep, tolerance, &crate::serve_region::ServeRegion::Whole)
}

/// [`apply_face_simplification`] among the faces `region` admits (a
/// merge kills a face and creates none, so the region stays exact).
pub(crate) fn apply_face_simplification_in(
    brep: &mut BRep,
    tolerance: &Tolerance,
    region: &crate::serve_region::ServeRegion,
) -> usize {
    let mut eliminated = 0;
    // Hard cap so a pathological case can't run forever. Each merge
    // strictly decreases the face count, so the actual cap is the
    // initial face count anyway.
    let max_iterations = brep.faces.len() + 1;

    for _ in 0..max_iterations {
        let report = find_merge_candidates_in(brep, tolerance, region);
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
            let merged = brep
                .kef(he_h)
                .ok()
                .is_some_and(|k| rechart_doubly_periodic(brep, k.surviving_loop));
            if merged && brep.validate().is_ok() {
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
