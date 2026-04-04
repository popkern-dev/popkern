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

/// Apply face simplification: merge faces within each group by removing
/// shared edges via KEF. Validates after each removal; stops a group
/// if validation fails.
///
/// Returns the number of faces eliminated.
pub fn apply_face_simplification(brep: &mut BRep, tolerance: &Tolerance) -> usize {
    let report = find_merge_candidates(brep, tolerance);
    let mut eliminated = 0;

    for group in &report.merge_groups {
        if group.len() < 2 {
            continue;
        }

        // Find shared edges between faces in this group
        let group_set: std::collections::HashSet<FaceHandle> = group.iter().copied().collect();

        let shared_edges: Vec<_> = brep
            .edges
            .keys()
            .filter(|&eh| {
                let edge = &brep.edges[eh];
                let he = match brep.half_edges.get(edge.half_edge) {
                    Some(h) => h,
                    None => return false,
                };
                let twin = match brep.half_edges.get(he.twin) {
                    Some(t) => t,
                    None => return false,
                };
                let fa = match brep.loops.get(he.loop_) {
                    Some(l) => l.face,
                    None => return false,
                };
                let fb = match brep.loops.get(twin.loop_) {
                    Some(l) => l.face,
                    None => return false,
                };
                fa != fb && group_set.contains(&fa) && group_set.contains(&fb)
            })
            .collect();

        // Remove shared edges one at a time with validation
        for eh in shared_edges {
            let he = match brep.edges.get(eh) {
                Some(e) => e.half_edge,
                None => continue,
            };

            if brep.kef(he).is_ok() {
                if brep.validate().is_err() {
                    // Topology corrupted — stop this group
                    break;
                }
                eliminated += 1;
            }
        }
    }

    eliminated
}
