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

        let mut queue = std::collections::VecDeque::new();
        queue.push_back(fh);

        while let Some(current) = queue.pop_front() {
            // Find adjacent faces via shared edges
            for (_eh, edge) in &brep.edges {
                let he = match brep.half_edges.get(edge.half_edge) {
                    Some(h) => h,
                    None => continue,
                };
                let twin = match brep.half_edges.get(he.twin) {
                    Some(t) => t,
                    None => continue,
                };

                let loop_a = match brep.loops.get(he.loop_) {
                    Some(l) => l,
                    None => continue,
                };
                let loop_b = match brep.loops.get(twin.loop_) {
                    Some(l) => l,
                    None => continue,
                };

                let neighbor = if loop_a.face == current && !visited.contains(&loop_b.face) {
                    loop_b.face
                } else if loop_b.face == current && !visited.contains(&loop_a.face) {
                    loop_a.face
                } else {
                    continue;
                };

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
