//! §8.3 blend cap **insertion** — post-construct pass that decides
//! what to do with the cap geometry already computed in
//! `EdgeGeometry::start_cap` / `end_cap`.
//!
//! For each fillet edge with a non-`Open` cap kind, this pass
//! determines whether the endpoint is **closed** (already has an
//! adjacent solid face that serves as the natural cap) or **free**
//! (open shell — no neighbouring face). On closed endpoints the
//! cap is implicit; on free endpoints the cap surface is a
//! candidate for topology insertion as a new face.
//!
//! Today the insertion path is a stub (no open-shell tests exist
//! yet because popkern has no sheet-body constructor), but the
//! detection + warning pipeline is fully wired so callers can
//! observe the §8.3 semantics end-to-end on existing solid
//! configurations and so tests can pin the expected behaviour.
//!
//! Emits [`BlendWarning::CapApplied`] once per `(edge, endpoint)`
//! pair whenever `cap_kind != Open`, with `topology_inserted` set
//! to reflect the chosen path.

use popkern_core::geometry::surface::SurfaceOnFace;
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

use crate::construct::EdgeGeometry;
use crate::error::BlendWarning;

/// Post-construct cap pass.
///
/// Walks every edge in `edge_geometries`; for each with a non-None
/// `start_cap` or `end_cap`, decides "closed" vs "free", and emits
/// the corresponding warning. Currently no topology mutation
/// occurs on either path — closed-body cap semantics are already
/// correct via the adjacent face, and open-shell insertion is a
/// tracked follow-up gated on sheet-body constructor support.
///
/// This is the final §8.3 pipeline hook and should be called after
/// the consumed-merge stage (face merging), vertex-blend handling, and surface-fix
/// correction in `apply_blend_geometry`.
pub(crate) fn insert_caps(
    brep: &mut BRep,
    original_brep: &BRep,
    edge_geometries: &[EdgeGeometry],
    warnings: &mut Vec<(popkern_core::topology::EdgeHandle, BlendWarning)>,
) {
    for eg in edge_geometries {
        process_endpoint(
            brep,
            original_brep,
            eg,
            /* at_start = */ true,
            warnings,
        );
        process_endpoint(
            brep,
            original_brep,
            eg,
            /* at_start = */ false,
            warnings,
        );
    }
}

fn process_endpoint(
    brep: &mut BRep,
    original_brep: &BRep,
    eg: &EdgeGeometry,
    at_start: bool,
    warnings: &mut Vec<(popkern_core::topology::EdgeHandle, BlendWarning)>,
) {
    let cap = if at_start {
        eg.start_cap.as_ref()
    } else {
        eg.end_cap.as_ref()
    };
    let Some(cap) = cap else {
        return;
    };
    // `Open` shouldn't reach here (commit 2 drops it in
    // `build_edge_caps`), and `UserDefined` is treated as open, but
    // belt-and-braces: skip anything that somehow still carries a
    // `None` surface through.
    if cap.surface.is_none() {
        return;
    }

    let vertex = if at_start {
        eg.analysis.vertex_start
    } else {
        eg.analysis.vertex_end
    };

    // Free-boundary detection: an endpoint is "closed" when the
    // vertex is incident to at least one face OTHER than the two
    // support faces, which is the common case for closed solid
    // bodies (a cube corner has 3 faces; two of them are the
    // supports and the third implicitly caps the blend). On an
    // open shell — e.g. a sheet body where the vertex has only
    // the two supports — the endpoint is free.
    let closed = endpoint_is_closed(
        original_brep,
        vertex,
        eg.analysis.face_1,
        eg.analysis.face_2,
    );

    let topology_inserted = if closed {
        false
    } else {
        insert_cap_face(brep, vertex, cap)
    };

    warnings.push((
        eg.analysis.edge,
        BlendWarning::CapApplied {
            edge: eg.analysis.edge,
            at_start,
            kind: cap.kind,
            topology_inserted,
        },
    ));
}

/// Returns `true` if `vertex` has at least one incident face that
/// is NOT `support_1` or `support_2`.
///
/// The implementation walks every half-edge in the BRep whose
/// `origin` equals `vertex` and collects the set of incident
/// faces via `loop_.face`. As soon as a face other than the two
/// supports shows up, the endpoint is closed. An endpoint that
/// sees only the two supports (or fewer) is free.
fn endpoint_is_closed(
    brep: &BRep,
    vertex: VertexHandle,
    support_1: popkern_core::topology::FaceHandle,
    support_2: popkern_core::topology::FaceHandle,
) -> bool {
    for he in brep.half_edges.values() {
        if he.origin != vertex {
            continue;
        }
        let face = brep.loops[he.loop_].face;
        if face == support_1 || face == support_2 {
            continue;
        }
        // Faces with no surface geometry are "virtual" backside
        // faces (e.g. sheet bodies, or a solid face whose
        // surface has been cleared to simulate an open shell).
        // They don't count as a real implicit cap.
        if brep.faces[face].surface.is_none() {
            continue;
        }
        return true;
    }
    false
}

/// Insert a cap face at a free-boundary blend endpoint.
///
/// Finds the two BRep vertices nearest to the cap's first and
/// last boundary points (the contact vertices on the two
/// support faces), locates a common face loop containing both,
/// and calls `mef()` to split that face along a new edge
/// between the contact vertices. The new face is assigned the
/// cap surface. Returns `true` on success.
///
/// The search is confined to faces incident to `vertex` (the
/// original edge endpoint), which keeps the nearest-vertex
/// lookup local. If neither contact can be matched, or if no
/// common face is found, the function returns `false` without
/// modifying topology — the warning will report
/// `topology_inserted = false`.
fn insert_cap_face(brep: &mut BRep, vertex: VertexHandle, cap: &crate::cap::CapGeometry) -> bool {
    let Some(surface) = &cap.surface else {
        return false;
    };

    // Find the two vertices adjacent to `vertex` on a face with
    // a surface (the fillet/torus face). These are the trim
    // curve endpoints — the points where the fillet touches
    // each support face at this endpoint. We derive them from
    // the BRep topology (which is correct after the consumed-merge stage) rather
    // than from `cap.boundary_points` (which can overshoot due
    // to footprint/center-curve extrapolation).
    // Find the 2 nearest BRep vertices to the endpoint that sit
    // on a face with surface geometry. The fillet strip (torus
    // face) has its trim endpoints at these positions — they're
    // the cap's contact points. Direct adjacency from the
    // endpoint vertex doesn't work because the boundary-edge
    // removal + fillet trim can leave the vertex with only 1
    // outgoing half-edge.
    let ep_pos = brep.vertices[vertex].position;
    let mut candidates: Vec<(VertexHandle, f64)> = brep
        .vertices
        .iter()
        .filter_map(|(vh, v)| {
            if vh == vertex {
                return None;
            }
            // Must be on at least one face with surface geometry.
            let on_surfaced = brep.half_edges.values().any(|he| {
                he.origin == vh && brep.faces[brep.loops[he.loop_].face].surface.is_some()
            });
            if !on_surfaced {
                return None;
            }
            let d2 = (v.position.x - ep_pos.x).powi(2)
                + (v.position.y - ep_pos.y).powi(2)
                + (v.position.z - ep_pos.z).powi(2);
            Some((vh, d2))
        })
        .collect();
    candidates.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    if candidates.len() < 2 {
        return false;
    }
    let v1 = candidates[0].0;
    let v2 = candidates[1].0;
    if v1 == v2 {
        return false;
    }

    let Some((he1, he2)) = find_common_loop_pair(brep, v1, v2) else {
        return false;
    };

    let Ok(mef_r) = brep.mef(he1, he2) else {
        return false;
    };

    brep.faces[mef_r.face].surface = Some(SurfaceOnFace {
        geometry: surface.clone(),
        orientation: true,
    });

    true
}

/// Find a pair of half-edges in the same face loop with origins
/// `v1` and `v2`, suitable for passing to `mef()`.
fn find_common_loop_pair(
    brep: &BRep,
    v1: VertexHandle,
    v2: VertexHandle,
) -> Option<(HalfEdgeHandle, HalfEdgeHandle)> {
    // Collect all half-edges with origin v1.
    let hes_v1: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.origin == v1)
        .map(|(h, _)| h)
        .collect();
    // For each, walk its loop to find one with origin v2.
    for he1 in hes_v1 {
        let loop_h = brep.half_edges[he1].loop_;
        let start = brep.loops[loop_h].half_edge;
        let mut cursor = start;
        loop {
            if brep.half_edges[cursor].origin == v2 && cursor != he1 {
                return Some((he1, cursor));
            }
            cursor = brep.half_edges[cursor].next;
            if cursor == start {
                break;
            }
        }
    }
    None
}

#[cfg(test)]
#[path = "cap_insert_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "cap_insert_integration_tests.rs"]
mod integration_tests;
