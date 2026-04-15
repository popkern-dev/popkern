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

use popkern_core::topology::{BRep, VertexHandle};

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
    brep: &BRep,
    edge_geometries: &[EdgeGeometry],
    warnings: &mut Vec<(popkern_core::topology::EdgeHandle, BlendWarning)>,
) {
    for eg in edge_geometries {
        process_endpoint(brep, eg, /* at_start = */ true, warnings);
        process_endpoint(brep, eg, /* at_start = */ false, warnings);
    }
}

fn process_endpoint(
    brep: &BRep,
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
        // The adjacent face serves as the natural cap; no
        // topology change is needed or desirable — modifying a
        // closed solid's existing faces would create a
        // non-manifold configuration.
        false
    } else {
        // Open-shell / free-boundary endpoint. Topology insertion
        // would create a new face from the `CapGeometry.surface`
        // and stitch its boundary to the existing trim edges at
        // the two contact points. Implementing that surgery on
        // open shells is a tracked follow-up; for now we record
        // the intent and leave topology unchanged.
        insert_cap_face_stub(cap)
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

/// Placeholder topology insertion for open-shell / free-boundary
/// cap faces.
///
/// **Status.** Sheet-body constructors landed (`BRep::sheet_quad`,
/// `sheet_triangle`, `sheet_polygon`) so the INPUTS a real
/// implementation needs are now buildable. The remaining blocker
/// is that `popkern-fillet`'s edge analyzer currently rejects
/// fillets on planar sheet bodies — the front and backside faces
/// of a `sheet_quad` are coincident planes with opposite normals,
/// which the convexity classifier treats as degenerate. That's
/// the §5.7 "2D fillet / sheet fillet" scope on the requirements
/// gap list; opening it up is the next step before this function
/// can be turned from a no-op into the real surgery.
///
/// The signature already takes a full [`crate::cap::CapGeometry`],
/// so the eventual implementation will slot in without any
/// caller-visible changes: it will create a new [`Face`] with
/// `cap.surface`, stitch boundary half-edges between the two
/// contact vertices using `cap.boundary_points`, and splice the
/// result into the blend face's outer loop.
///
/// [`Face`]: popkern_core::topology::Face
fn insert_cap_face_stub(_cap: &crate::cap::CapGeometry) -> bool {
    // Insertion is a follow-up. The function takes `&cap` so the
    // signature is forward-compatible with the real implementation.
    false
}

#[cfg(test)]
#[path = "cap_insert_tests.rs"]
mod tests;
