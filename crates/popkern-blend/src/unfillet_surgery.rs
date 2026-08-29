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

/// Find ALL blend-face loop vertices that are "near" `target` —
/// the cube corner vertex (often at distance 0 because the
/// construct phase leaves it in the loop) plus the perpendicular-
/// edge split vertices on each support face (at distance ≈ R from
/// the cube corner). Used by [`crate::unfillet::unfillet_face`] to
/// weld the entire cluster of vertices at one fillet endpoint into
/// a single canonical vertex, eliminating leftover perp-split
/// vertices that the old `corner_vertex_pair_near`-based pipeline
/// missed.
///
/// Threshold is `2 *` the second-closest distance, which catches
/// vertices roughly at the same fillet-radius distance from
/// `target` while excluding the far-side endpoint and unrelated
/// vertices.
pub(crate) fn corners_near_endpoint(
    brep: &BRep,
    blend_face: FaceHandle,
    target: &Point3,
) -> Result<Vec<VertexHandle>, BlendError> {
    let candidates = sorted_loop_candidates(brep, blend_face, target)?;
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
    // Use the second-closest distance² as the radius² scale. The
    // closest is typically a cube corner at distance 0; the second
    // is one of the perp-split vertices at distance ≈ R, so ≈ R².
    // Anything within 2× that distance is a sibling perp-split
    // vertex on the same fillet endpoint.
    let r2 = candidates[1].1;
    let threshold = (r2 * 4.0).max(1e-12);
    Ok(candidates
        .into_iter()
        .filter(|(_, d2)| *d2 <= threshold)
        .map(|(v, _)| v)
        .collect())
}

/// Walk the blend face's outer loop and return every distinct
/// vertex with its squared distance to `target`, sorted ascending.
fn sorted_loop_candidates(
    brep: &BRep,
    blend_face: FaceHandle,
    target: &Point3,
) -> Result<Vec<(VertexHandle, f64)>, BlendError> {
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
    Ok(candidates)
}

/// Walk every edge in the BRep and, if the edge touches one of the
/// welded vertices, rebuild the curve to span the edge's *current*
/// vertex positions as a straight line. Closes the gap that vertex
/// repositioning leaves in loop validation: edges whose endpoints
/// moved but whose curve geometry kept the old start/end points.
///
/// Replaces *any* curve type with a fresh straight line — for
/// the initial plane-plane scope all the affected edges are
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
    // One edge per pass, re-scanning after each removal: removing a
    // degenerate edge rewrites its neighbours' next/prev pointers,
    // and a batch collected up front dereferences stale handles when
    // degenerate edges form CHAINS (the collapsed cap ring of a
    // corner-group unfillet panicked on an already-removed prev).
    while let Some(eh) = brep.edges.iter().find_map(|(eh, edge)| {
        let he = brep.half_edges.get(edge.half_edge)?;
        let twin = brep.half_edges.get(he.twin)?;
        (he.origin == twin.origin).then_some(eh)
    }) {
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

        if he.next == twin_handle && twin.next == he_handle {
            // Isolated 2-cycle: the pair IS its whole loop. Nothing
            // outside references the half-edges; drop them and leave
            // the (empty) loop for the caller's loop-level cleanup.
            brep.remove_half_edge_reanchored(he_handle);
            brep.remove_half_edge_reanchored(twin_handle);
            brep.edges.remove(eh);
            continue;
        }
        if he.next == twin_handle || twin.next == he_handle {
            // Zero-length SPUR: the pair sits consecutively in one
            // loop (a collapsed cap-ring arc after the corner weld).
            // Splice the pair out as a unit — the generic
            // one-at-a-time splice below writes through the partner's
            // stale pointers and dangles the neighbours (measured:
            // invalid-SlotMap-key panic on the corner-group restore).
            let (first, second) = if he.next == twin_handle {
                (he_handle, twin_handle)
            } else {
                (twin_handle, he_handle)
            };
            let before = brep.half_edges[first].prev;
            let after = brep.half_edges[second].next;
            brep.half_edges[before].next = after;
            brep.half_edges[after].prev = before;
            for h in [first, second] {
                let lp = brep.half_edges[h].loop_;
                if let Some(loop_data) = brep.loops.get(lp)
                    && loop_data.half_edge == h
                {
                    brep.loops[lp].half_edge = after;
                }
            }
            brep.remove_half_edge_reanchored(he_handle);
            brep.remove_half_edge_reanchored(twin_handle);
            brep.edges.remove(eh);
            continue;
        }

        // Separated pair: splice each side out of its loop.
        brep.half_edges[he.prev].next = he.next;
        brep.half_edges[he.next].prev = he.prev;
        brep.half_edges[twin.prev].next = twin.next;
        brep.half_edges[twin.next].prev = twin.prev;
        // Update each loop's `half_edge` pointer if it pointed at
        // one of the removed half-edges.
        if let Some(loop_data) = brep.loops.get(he.loop_)
            && loop_data.half_edge == he_handle
        {
            brep.loops[he.loop_].half_edge = he.next;
        }
        if let Some(loop_data) = brep.loops.get(twin.loop_)
            && loop_data.half_edge == twin_handle
        {
            brep.loops[twin.loop_].half_edge = twin.next;
        }

        brep.remove_half_edge_reanchored(he_handle);
        brep.remove_half_edge_reanchored(twin_handle);
        brep.edges.remove(eh);
    }
    // Boundary (twin-less) edges collapsed to zero length by the
    // welds: splice the lone half-edge out and retire the edge.
    // Twin-less edges only exist after tolerant residual removal
    // (corner-group unfillet); paired-edge paths never enter here.
    while let Some((eh, he_handle)) = brep.edges.iter().find_map(|(eh, edge)| {
        let he = brep.half_edges.get(edge.half_edge)?;
        if brep.half_edges.contains_key(he.twin) {
            return None;
        }
        let dest = brep.half_edges.get(he.next)?.origin;
        (he.origin == dest).then_some((eh, edge.half_edge))
    }) {
        let he = brep.half_edges[he_handle].clone();
        if he.next != he_handle {
            brep.half_edges[he.prev].next = he.next;
            brep.half_edges[he.next].prev = he.prev;
            if let Some(loop_data) = brep.loops.get(he.loop_)
                && loop_data.half_edge == he_handle
            {
                brep.loops[he.loop_].half_edge = he.next;
            }
        }
        brep.remove_half_edge_reanchored(he_handle);
        brep.edges.remove(eh);
        brep.boundary_edges.remove(&eh);
    }
}

/// Remove residual zero-area slit spurs: an edge whose two
/// half-edges sit CONSECUTIVELY in one loop (`he.next == twin` —
/// the walk goes out and straight back along the same span, a dead
/// end bounding no area regardless of length). The corner-group
/// merge strands one per redundant contact edge on the merged
/// region — the mixed-radius DirectApex corner class (PK-255 serve)
/// carries more group contact edges than the uniform 6-ring, and
/// the leftover slit survived the mef re-split as a spur embedded
/// in a restored support plane's outer loop, geometrically
/// coincident with a real restored edge (measured: restored "cube"
/// with 13 edges). A spur built from ONE edge's own half-edge pair
/// is always zero-area, so removal never changes face geometry; a
/// bridge to an inner island has `he.next != twin` and is never
/// touched. Zero-length self-loop pairs remain
/// [`remove_degenerate_self_loop_edges`]'s job.
pub(crate) fn remove_slit_loop_edges(brep: &mut BRep) {
    while let Some((eh, first, second)) = brep.edges.iter().find_map(|(eh, edge)| {
        let he = brep.half_edges.get(edge.half_edge)?;
        let twin = brep.half_edges.get(he.twin)?;
        if he.loop_ != twin.loop_ || he.origin == twin.origin {
            return None;
        }
        if he.next == he.twin {
            Some((eh, edge.half_edge, he.twin))
        } else if twin.next == edge.half_edge {
            Some((eh, he.twin, edge.half_edge))
        } else {
            None
        }
    }) {
        let loop_handle = brep.half_edges[first].loop_;
        let before = brep.half_edges[first].prev;
        let after = brep.half_edges[second].next;
        if before == second && after == first {
            // The pair IS the whole loop: an isolated slit loop.
            // Retire the loop record and its face registration.
            let face = brep.loops[loop_handle].face;
            if let Some(face_data) = brep.faces.get_mut(face) {
                face_data
                    .inner_loops
                    .retain(|&candidate| candidate != loop_handle);
            }
            brep.loops.remove(loop_handle);
        } else {
            brep.half_edges[before].next = after;
            brep.half_edges[after].prev = before;
            if let Some(loop_data) = brep.loops.get(loop_handle)
                && (loop_data.half_edge == first || loop_data.half_edge == second)
            {
                brep.loops[loop_handle].half_edge = after;
            }
        }
        for handle in [first, second] {
            let vertex = brep.half_edges[handle].origin;
            if brep
                .vertices
                .get(vertex)
                .is_some_and(|data| data.outgoing == handle)
                && let Some(replacement) = brep
                    .half_edges
                    .iter()
                    .find(|(candidate, data)| {
                        data.origin == vertex && *candidate != first && *candidate != second
                    })
                    .map(|(candidate, _)| candidate)
            {
                brep.vertices[vertex].outgoing = replacement;
            }
        }
        brep.remove_half_edge_reanchored(first);
        brep.remove_half_edge_reanchored(second);
        brep.edges.remove(eh);
        brep.boundary_edges.remove(&eh);
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
