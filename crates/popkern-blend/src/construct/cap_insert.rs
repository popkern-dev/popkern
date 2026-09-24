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
//! On a free end, `insert_cap_face` splits the strip's end off with
//! `mef` between the two contacts, gives the corner side the cap
//! surface, and `cap_insert_boundary` authors the exact boundary.
//! Flat and Rounded caps insert. A cap that cannot be made exact
//! (TangentMatch's bilinear patch, a rim outside a flat cap's plane)
//! is rolled back, reported as `CapInsertionDeclined`, and the end is
//! served open (PK-532).
//!
//! Emits [`BlendWarning::CapApplied`] once per `(edge, endpoint)`
//! pair whenever `cap_kind != Open`, with `topology_inserted` set
//! to reflect the chosen path.

use popkern_core::geometry::surface::SurfaceOnFace;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

use crate::construct::EdgeGeometry;
use crate::error::BlendWarning;

/// Post-construct cap pass.
///
/// Walks every edge in `edge_geometries`; for each with a non-None
/// `start_cap` or `end_cap`, decides "closed" vs "free", and emits
/// the corresponding warning. Closed ends are capped by their
/// adjacent face already; free ends get the cap inserted with an
/// exact boundary, or a `CapInsertionDeclined` fallback.
///
/// This is the final §8.3 pipeline hook and should be called after
/// the consumed-merge stage (face merging), vertex-blend handling, and surface-fix
/// correction in `apply_blend_geometry`.
pub(crate) fn insert_caps(
    brep: &mut BRep,
    original_brep: &BRep,
    edge_geometries: &[EdgeGeometry],
    warnings: &mut Vec<(popkern_core::topology::EdgeHandle, BlendWarning)>,
    tolerance: &Tolerance,
) {
    for eg in edge_geometries {
        process_endpoint(
            brep,
            original_brep,
            eg,
            /* at_start = */ true,
            warnings,
            tolerance,
        );
        process_endpoint(
            brep,
            original_brep,
            eg,
            /* at_start = */ false,
            warnings,
            tolerance,
        );
    }
}

fn process_endpoint(
    brep: &mut BRep,
    original_brep: &BRep,
    eg: &EdgeGeometry,
    at_start: bool,
    warnings: &mut Vec<(popkern_core::topology::EdgeHandle, BlendWarning)>,
    tolerance: &Tolerance,
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
        match insert_cap_face(brep, vertex, cap, tolerance) {
            Ok(inserted) => inserted,
            Err(reason) => {
                warnings.push((
                    eg.analysis.edge,
                    BlendWarning::CapInsertionDeclined {
                        edge: eg.analysis.edge,
                        at_start,
                        kind: cap.kind,
                        reason,
                    },
                ));
                false
            }
        }
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
fn insert_cap_face(
    brep: &mut BRep,
    vertex: VertexHandle,
    cap: &crate::cap::CapGeometry,
    tolerance: &Tolerance,
) -> Result<bool, String> {
    let Some(surface) = &cap.surface else {
        return Ok(false);
    };
    if cap.kind == crate::types::CapKind::TangentMatch {
        // Its bilinear patch has a straight contact-to-contact
        // boundary: it cannot contain the strip's end section arc
        // (PK-532).
        return Err("the tangent-match patch cannot contain the strip's end section".into());
    }

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
        return Ok(false);
    }
    let v1 = candidates[0].0;
    let v2 = candidates[1].0;
    if v1 == v2 {
        return Ok(false);
    }

    let Some((he1, he2)) = find_common_loop_pair(brep, v1, v2) else {
        return Ok(false);
    };

    let snapshot = brep.clone();
    let Ok(mef_r) = brep.mef(he1, he2) else {
        return Ok(false);
    };
    // The cap is the side holding the original corner (the rim L);
    // the other side keeps the strip surface both inherited.
    let corner_side = |b: &BRep, he: HalfEdgeHandle| -> bool {
        let Some(start) = b.half_edges.get(he) else {
            return false;
        };
        let mut cursor = start.next;
        for _ in 0..b.half_edges.len() {
            if cursor == he {
                return false;
            }
            match b.half_edges.get(cursor) {
                Some(h) if h.origin == vertex => return true,
                Some(h) => cursor = h.next,
                None => return false,
            }
        }
        false
    };
    let cap_half_edge = if corner_side(brep, mef_r.half_edge) {
        mef_r.half_edge
    } else if corner_side(brep, mef_r.half_edge_twin) {
        mef_r.half_edge_twin
    } else {
        *brep = snapshot;
        return Err("neither side of the cap split holds the corner".into());
    };
    let cap_face = brep.loops[brep.half_edges[cap_half_edge].loop_].face;
    brep.faces[cap_face].surface = Some(SurfaceOnFace {
        geometry: surface.clone(),
        orientation: true,
    });
    if let Err(reason) =
        super::cap_insert_boundary::author_cap_boundary(brep, cap_half_edge, vertex, cap, tolerance)
    {
        *brep = snapshot;
        return Err(reason);
    }
    Ok(true)
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
