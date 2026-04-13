//! Topology-surgery helpers for [`crate::unfillet`].
//!
//! Holds the post-KEF cleanup primitives that the unfillet pipeline
//! needs after the Euler operators run: vertex-pair detection,
//! degenerate self-loop edge removal, curve-geometry rebuild after
//! vertex repositioning, and loop traversal helpers. Extracted from
//! `unfillet.rs` to keep each source file under the 500-line
//! structural cap.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::error::BlendError;

/// Find the two vertices on the blend face's outer loop that are
/// closest to the given 3D position. Returns them sorted by
/// distance (closest first). Errors if fewer than 2 distinct
/// vertices are reachable.
pub(crate) fn corner_vertex_pair_near(
    brep: &BRep,
    blend_face: FaceHandle,
    target: &Point3,
) -> Result<Vec<VertexHandle>, BlendError> {
    let face_data = brep
        .faces
        .get(blend_face)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "blend_face".to_string(),
            value: format!("{blend_face:?}"),
            reason: "face handle does not exist".to_string(),
        })?;
    let loop_data =
        brep.loops
            .get(face_data.outer_loop)
            .ok_or_else(|| BlendError::InvalidParameter {
                parameter: "blend_face.outer_loop".to_string(),
                value: format!("{:?}", face_data.outer_loop),
                reason: "loop handle does not exist".to_string(),
            })?;
    let start = loop_data.half_edge;

    let mut candidates: Vec<(VertexHandle, f64)> = Vec::new();
    let mut he = start;
    let mut safety = 0;
    loop {
        safety += 1;
        if safety > 10_000 {
            break;
        }
        let Some(he_data) = brep.half_edges.get(he) else {
            break;
        };
        let v = he_data.origin;
        if !candidates.iter().any(|(c, _)| *c == v) {
            let pos = brep.vertices[v].position;
            let d2 = (pos.x - target.x).powi(2)
                + (pos.y - target.y).powi(2)
                + (pos.z - target.z).powi(2);
            candidates.push((v, d2));
        }
        he = he_data.next;
        if he == start {
            break;
        }
    }
    candidates.sort_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    if candidates.len() < 2 {
        return Err(BlendError::InvalidParameter {
            parameter: "blend_face".to_string(),
            value: format!("{blend_face:?}"),
            reason: format!(
                "blend face has only {} distinct boundary vertices; need ≥2 for unfillet",
                candidates.len()
            ),
        });
    }
    Ok(vec![candidates[0].0, candidates[1].0])
}

/// Walk every edge in the BRep and, if the edge touches one of the
/// welded vertices, rebuild the curve to span the edge's *current*
/// vertex positions as a straight line. Closes the gap that vertex
/// repositioning leaves in loop validation: edges whose endpoints
/// moved but whose curve geometry kept the old start/end points.
///
/// Replaces *any* curve type with a fresh straight line — for
/// Round 5c-1's plane-plane scope all the affected edges are
/// straight cube edges or their fillet-corner replacements, so a
/// linear curve is exact.
pub(crate) fn rebuild_curves_touching(brep: &mut BRep, welded: &[VertexHandle]) {
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};

    let mut to_rebuild: Vec<EdgeHandle> = Vec::new();
    for (eh, edge) in brep.edges.iter() {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let touches_welded = welded.contains(&he.origin) || welded.contains(&twin.origin);
        if touches_welded {
            to_rebuild.push(eh);
        }
    }

    for eh in to_rebuild {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let he_handle = edge.half_edge;
        let Some(he) = brep.half_edges.get(he_handle).cloned() else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin).cloned() else {
            continue;
        };
        let Some(start_v) = brep.vertices.get(he.origin) else {
            continue;
        };
        let Some(end_v) = brep.vertices.get(twin.origin) else {
            continue;
        };
        let start_pos = start_v.position;
        let end_pos = end_v.position;
        let direction = Vec3::new(
            end_pos.x - start_pos.x,
            end_pos.y - start_pos.y,
            end_pos.z - start_pos.z,
        );
        let length = direction.length();
        if length < 1e-12 {
            // Degenerate (should have been removed already). Skip.
            continue;
        }

        brep.edges[eh].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: start_pos,
                direction,
            },
            t_start: 0.0,
            t_end: 1.0,
        });
    }
}

/// Remove every edge in the BRep whose two half-edges have the
/// same `origin` vertex (a "self-loop" edge). After
/// [`BRep::weld_vertex_pair`] collapses two vertices into one,
/// every edge that previously connected those two vertices becomes
/// a self-loop and contributes garbage to the V−E+F invariant.
///
/// For each such edge:
/// 1. Splice both half-edges out of their respective loops
/// 2. Update the loop's `half_edge` pointer if it was one of the
///    removed half-edges
/// 3. Delete the edge and both half-edges from the slotmaps
pub(crate) fn remove_degenerate_self_loop_edges(brep: &mut BRep) {
    // Collect edges to remove first to avoid mutating during iteration.
    let mut to_remove: Vec<EdgeHandle> = Vec::new();
    for (eh, edge) in brep.edges.iter() {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        if he.origin == twin.origin {
            to_remove.push(eh);
        }
    }

    for eh in to_remove {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let he_handle = edge.half_edge;
        let Some(he) = brep.half_edges.get(he_handle).cloned() else {
            continue;
        };
        let twin_handle = he.twin;
        let Some(twin) = brep.half_edges.get(twin_handle).cloned() else {
            continue;
        };

        // Splice he out of its loop.
        if he.next != he_handle {
            brep.half_edges[he.prev].next = he.next;
            brep.half_edges[he.next].prev = he.prev;
        }
        // Splice twin out of its loop.
        if twin.next != twin_handle {
            brep.half_edges[twin.prev].next = twin.next;
            brep.half_edges[twin.next].prev = twin.prev;
        }
        // Update each loop's `half_edge` pointer if it pointed at
        // one of the removed half-edges.
        if let Some(loop_data) = brep.loops.get(he.loop_)
            && loop_data.half_edge == he_handle
            && he.next != he_handle
        {
            brep.loops[he.loop_].half_edge = he.next;
        }
        if let Some(loop_data) = brep.loops.get(twin.loop_)
            && loop_data.half_edge == twin_handle
            && twin.next != twin_handle
        {
            brep.loops[twin.loop_].half_edge = twin.next;
        }

        brep.half_edges.remove(he_handle);
        brep.half_edges.remove(twin_handle);
        brep.edges.remove(eh);
    }
}

/// Walk a face's outer loop and return the first half-edge whose
/// `origin` is the given vertex.
pub(crate) fn find_he_in_face_with_origin(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
) -> Option<HalfEdgeHandle> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    let mut he = start;
    let mut safety = 0;
    loop {
        safety += 1;
        if safety > 10_000 {
            return None;
        }
        let he_data = brep.half_edges.get(he)?;
        if he_data.origin == vertex {
            return Some(he);
        }
        he = he_data.next;
        if he == start {
            return None;
        }
    }
}
