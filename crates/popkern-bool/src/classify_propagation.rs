//! Topological propagation for face classification.
//!
//! Instead of ray-casting every face independently, this classifier finds
//! connected components of faces sharing non-intersection edges, classifies
//! one seed face per component, then propagates the classification to all
//! faces in the component via BFS.
//!
//! This reduces the number of expensive ray-cast calls from O(n_faces) to
//! O(n_components), which is typically a 50-90% reduction.

use std::collections::{HashMap, HashSet, VecDeque};

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::classify::point_in_solid;
use crate::error::{BoolError, Classification, PointClassification};

/// Classify faces using topological propagation.
///
/// 1. Build a face-adjacency graph (faces sharing non-intersection edges).
/// 2. Find connected components via BFS.
/// 3. Classify one seed face per component via ray casting.
/// 4. Propagate the seed's classification to all faces in its component.
///
/// Returns the per-face classification map plus the number of ray casts used.
///
/// # Arguments
///
/// * `brep` — the BRep whose faces are being classified
/// * `other_brep` — the solid to test containment against
/// * `other_aabbs` — precomputed face AABBs for the other solid
/// * `intersection_edges` — edges lying on the intersection (DO NOT propagate across these)
/// * `tolerance` — geometric tolerance
pub fn classify_with_propagation(
    brep: &BRep,
    other_brep: &BRep,
    other_aabbs: &[(FaceHandle, Aabb)],
    intersection_edges: &HashSet<EdgeHandle>,
    tolerance: &Tolerance,
) -> Result<(HashMap<FaceHandle, Classification>, usize), BoolError> {
    // Step 1: Build face adjacency (exclude intersection edges)
    let adjacency = build_face_adjacency(brep, intersection_edges);

    // Step 2: Find connected components
    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    let components = find_components(&face_handles, &adjacency);

    // Step 3+4: Classify seed, propagate
    let mut result = HashMap::new();
    let mut ray_cast_count = 0;

    for component in &components {
        if component.is_empty() {
            continue;
        }

        // Pick the seed face (first in component)
        let seed = component[0];
        let sample = crate::classify::pick_sample_point(brep, seed);
        let pc = point_in_solid(&sample, other_brep, other_aabbs, tolerance)?;
        ray_cast_count += 1;

        let class = match pc {
            PointClassification::Inside => Classification::Inside,
            PointClassification::Outside => Classification::Outside,
            PointClassification::OnBoundary => Classification::OnBoundaryCooriented,
        };

        // Propagate to all faces in the component
        for &face in component {
            result.insert(face, class);
        }
    }

    Ok((result, ray_cast_count))
}

/// Build a face-adjacency map: for each face, list the adjacent faces
/// reachable via non-intersection edges.
fn build_face_adjacency(
    brep: &BRep,
    intersection_edges: &HashSet<EdgeHandle>,
) -> HashMap<FaceHandle, Vec<FaceHandle>> {
    let mut adjacency: HashMap<FaceHandle, Vec<FaceHandle>> = HashMap::new();

    for (eh, edge) in &brep.edges {
        // Skip intersection edges — classification may differ across these
        if intersection_edges.contains(&eh) {
            continue;
        }

        // Get the two faces incident to this edge
        let he = match brep.half_edges.get(edge.half_edge) {
            Some(h) => h,
            None => continue,
        };
        let twin = match brep.half_edges.get(he.twin) {
            Some(t) => t,
            None => continue,
        };

        let face_a = match brep.loops.get(he.loop_) {
            Some(l) => l.face,
            None => continue,
        };
        let face_b = match brep.loops.get(twin.loop_) {
            Some(l) => l.face,
            None => continue,
        };

        if face_a != face_b {
            adjacency.entry(face_a).or_default().push(face_b);
            adjacency.entry(face_b).or_default().push(face_a);
        }
    }

    adjacency
}

/// Find connected components using BFS.
fn find_components(
    faces: &[FaceHandle],
    adjacency: &HashMap<FaceHandle, Vec<FaceHandle>>,
) -> Vec<Vec<FaceHandle>> {
    let mut visited = HashSet::new();
    let mut components = Vec::new();

    for &face in faces {
        if visited.contains(&face) {
            continue;
        }

        let mut component = Vec::new();
        let mut queue = VecDeque::new();
        queue.push_back(face);
        visited.insert(face);

        while let Some(current) = queue.pop_front() {
            component.push(current);

            if let Some(neighbors) = adjacency.get(&current) {
                for &neighbor in neighbors {
                    if visited.insert(neighbor) {
                        queue.push_back(neighbor);
                    }
                }
            }
        }

        components.push(component);
    }

    components
}
