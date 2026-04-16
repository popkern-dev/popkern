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
/// Phase C (face merging), vertex-blend handling, and surface-fix
/// correction in `apply_blend_geometry`.
pub(crate) fn insert_caps(
    brep: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    warnings: &mut Vec<(popkern_core::topology::EdgeHandle, BlendWarning)>,
) {
    for eg in edge_geometries {
        process_endpoint(brep, eg, /* at_start = */ true, warnings);
        process_endpoint(brep, eg, /* at_start = */ false, warnings);
    }
}

fn process_endpoint(
    brep: &mut BRep,
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
    let closed = endpoint_is_closed(brep, vertex, eg.analysis.face_1, eg.analysis.face_2);

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
        if face != support_1 && face != support_2 {
            return true;
        }
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
    if cap.boundary_points.len() < 2 {
        return false;
    }
    let target_1 = cap.boundary_points[0];
    let target_2 = cap.boundary_points[cap.boundary_points.len() - 1];

    // Find BRep vertices nearest to each contact, searching only
    // among vertices incident to `vertex`'s neighbouring faces.
    let tol_sq = 1.0e-4_f64; // generous tolerance²
    let v1 = find_nearest_vertex(brep, vertex, target_1, tol_sq);
    let v2 = find_nearest_vertex(brep, vertex, target_2, tol_sq);
    let (Some(v1), Some(v2)) = (v1, v2) else {
        return false;
    };
    if v1 == v2 {
        return false;
    }

    // Find a pair of half-edges in the SAME loop, one with origin
    // v1 and the other with origin v2.  `mef(he1, he2)` creates a
    // new edge from he1.origin to he2.origin.
    let Some((he1, he2)) = find_common_loop_pair(brep, v1, v2) else {
        return false;
    };

    // mef creates the new edge + face.
    let Ok(mef_r) = brep.mef(he1, he2) else {
        return false;
    };

    // Assign the cap surface to the new face.
    brep.faces[mef_r.face].surface = Some(SurfaceOnFace {
        geometry: surface.clone(),
        orientation: true,
    });

    true
}

/// Walk faces around `vertex` and return the vertex whose
/// position is closest to `target` (within `tol_sq` distance²).
fn find_nearest_vertex(
    brep: &BRep,
    vertex: VertexHandle,
    target: popkern_core::geometry::Point3,
    tol_sq: f64,
) -> Option<VertexHandle> {
    let mut best: Option<(VertexHandle, f64)> = None;
    for (vh, v) in brep.vertices.iter() {
        if vh == vertex {
            continue;
        }
        let dx = v.position.x - target.x;
        let dy = v.position.y - target.y;
        let dz = v.position.z - target.z;
        let d2 = dx * dx + dy * dy + dz * dz;
        if d2 < tol_sq && best.is_none_or(|(_, prev)| d2 < prev) {
            best = Some((vh, d2));
        }
    }
    best.map(|(vh, _)| vh)
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
