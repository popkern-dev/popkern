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
    let adjacency = build_face_adjacency(brep, intersection_edges, tolerance);

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
        if std::env::var_os("PK_BOOL_TRACE").is_some() {
            eprintln!(
                "  [propagate] seed {seed:?} (component of {}) sample ({:.3},{:.3},{:.3}) -> {pc:?}",
                component.len(),
                sample.x,
                sample.y,
                sample.z
            );
        }

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
    tolerance: &Tolerance,
) -> HashMap<FaceHandle, Vec<FaceHandle>> {
    let mut adjacency: HashMap<FaceHandle, Vec<FaceHandle>> = HashMap::new();

    let prop_trace = std::env::var_os("PK_PROP_TRACE").is_some();
    for (eh, edge) in &brep.edges {
        // Skip intersection edges — classification may differ across these
        if intersection_edges.contains(&eh) {
            if prop_trace {
                eprintln!("  [prop-barrier] {eh:?}");
            }
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

        // A sub-tolerance edge is a POINT contact, not shared
        // boundary: classification may legitimately differ across a
        // point (two regions of a split face meet at the junction
        // vertex where the intersection loop crosses a grid corner).
        // Propagating across one bridged inside/outside components
        // (measured: three ~1e-13 stub edges merged the torus disc
        // into the rest and the union kept every fragment). Skipping
        // only ever DISCONNECTS — the cost is one extra ray cast.
        let degenerate = {
            let a = brep.vertices.get(he.origin).map(|v| v.position);
            let b = brep
                .half_edges
                .get(he.twin)
                .and_then(|t| brep.vertices.get(t.origin))
                .map(|v| v.position);
            match (a, b) {
                (Some(a), Some(b)) if tolerance.points_coincident(&a, &b) => {
                    // Coincident endpoints alone do not make an edge
                    // degenerate — every topologically CLOSED edge
                    // (cylinder rim, seam circle) closes on one
                    // vertex (PK-296). The curve arbitrates: a real
                    // closed edge sweeps away from the vertex at its
                    // midpoint; residue has no curve or never leaves
                    // the point.
                    match edge.curve.as_ref() {
                        None => true,
                        Some(c) => {
                            use popkern_core::geometry::curve::CurveEval;
                            let mid = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
                            tolerance.points_coincident(&a, &mid)
                        }
                    }
                }
                _ => false,
            }
        };
        if degenerate {
            if prop_trace {
                eprintln!("  [prop-degen] {eh:?} skipped (point contact)");
            }
            continue;
        }
        if face_a != face_b {
            if prop_trace {
                let mid = brep
                    .vertices
                    .get(brep.half_edges[edge.half_edge].origin)
                    .map(|v| v.position);
                let mid2 = brep
                    .half_edges
                    .get(he.twin)
                    .and_then(|t| brep.vertices.get(t.origin))
                    .map(|v| v.position);
                if let (Some(a), Some(b)) = (mid, mid2) {
                    eprintln!(
                        "  [prop-adj] {face_a:?} <-> {face_b:?} via {eh:?} ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) len {:.3e}",
                        a.x,
                        a.y,
                        a.z,
                        b.x,
                        b.y,
                        b.z,
                        a.distance_to(&b)
                    );
                } else {
                    eprintln!("  [prop-adj] {face_a:?} <-> {face_b:?} via {eh:?}");
                }
            }
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
