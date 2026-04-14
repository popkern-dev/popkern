//! General post-construction sliver-face absorption pass.
//!
//! **Step 2 scope**: infrastructure only. This module subsumes the
//! previous specific-case `blend_cleanup::merge_corner_fragments_into_blends`
//! absorber (Class D — "≤4-vert Plane → Sphere vertex blend") under
//! a general classify-apply-reclassify framework that can be
//! extended to additional sliver classes in step 3 without further
//! plumbing changes.
//!
//! **Not yet in this module (step 3 work):**
//! - Class A: zero-area face detection + absorption
//! - Class B: thin (high aspect ratio) face detection + absorption
//! - Class C: redundant-edge detection (two adjacent faces sharing
//!   the same underlying surface) + absorption
//! - Emission of `SliverAbsorbWarning::Retained` when a classified
//!   sliver cannot be safely absorbed (topology would be violated)
//!
//! ## Architectural invariants
//!
//! **Restart-after-mutation** is mandatory: `BRep::kef` invalidates
//! `brep.faces.keys()` iteration. The general loop therefore
//! follows the exact pattern proved in the previous absorber —
//! classify, find ONE target, apply ONE `kef`, reclassify from
//! scratch. Never batch kef targets.
//!
//! **Directional kef**: when merging a plane fragment into a
//! sphere blend, we pass the sphere-side half-edge so the sphere
//! face *survives* and the plane fragment *dies*. Otherwise the
//! merged face would inherit the plane's surface and lose the
//! sphere labelling. Every absorber that touches surface-typed
//! neighbours must respect this convention.
//!
//! **Iteration cap** is `brep.faces.len() * 2 + 16`. Generous
//! safety net — cascade absorption rarely needs more than a few
//! passes in practice, and each successful absorption reduces the
//! face count, so the pass is guaranteed to terminate well below
//! the cap.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// A classified sliver face that the absorber decided to retain
/// (instead of absorbing) because removing it would violate
/// topology or no suitable neighbour exists.
///
/// Callers consume these via `BlendWarning` plumbing — the
/// `apply_blend_geometry` wiring converts them into
/// `BlendWarning::SliverRetained` (or similar) entries on the
/// result. **Step 2 ships the framework but does not yet emit any
/// retained slivers**: the only classified case today (Class D
/// corner fragments) is always absorbable, so the return list is
/// always empty. Step 3 fills in the remaining cases.
#[derive(Debug, Clone)]
pub(crate) enum SliverAbsorbWarning {
    /// A sliver was classified but could not be absorbed.
    #[allow(dead_code)] // step 3 will emit these
    Retained {
        /// The face that remains in the output BRep.
        face: FaceHandle,
        /// Which sliver category the classifier flagged.
        kind: SliverKind,
    },
}

/// Classification of a sliver face. Each variant carries the
/// measurement that drove the classification, for diagnostic
/// rendering when the sliver ends up retained.
#[derive(Debug, Clone)]
#[allow(dead_code)] // ZeroArea / Thin / RedundantEdge ship in step 3
pub(crate) enum SliverKind {
    /// Polygonal outer-loop area below `tolerance.linear()²`.
    ZeroArea {
        /// The measured polygonal area.
        area: f64,
    },
    /// Face with extreme aspect ratio (one dimension ≪ the other).
    Thin {
        /// Longer dimension.
        length: f64,
        /// Shorter dimension.
        width: f64,
    },
    /// Face shares an edge with a neighbour having the same
    /// underlying surface (duplicate subdivision artifact).
    RedundantEdge,
    /// ≤4-vert face whose twin face via some edge is a vertex-
    /// blend sphere — the Round 5r corner-fragment case.
    CornerFragment {
        /// Number of vertices on the outer loop (≤ 4 when flagged).
        vertex_count: usize,
    },
}

/// Run the general sliver absorption pass on a fillet output BRep.
///
/// **Step 2 behaviour**: only Class D (corner fragment → sphere
/// blend) is classified and absorbed. Classes A/B/C ship in step 3.
/// The iteration pattern, kef dispatch, and return type are final.
///
/// Classify → apply one `kef` → reclassify. Cap at
/// `brep.faces.len() * 2 + 16` iterations. On any `kef` failure
/// (stale handle, same-loop violation), the pass breaks out of the
/// loop rather than retrying — this matches the previous absorber's
/// behaviour and avoids infinite loops if the classifier drifts
/// from the surgery primitive's preconditions.
#[allow(dead_code)] // tolerance is reserved for step 3
pub(crate) fn absorb_slivers(brep: &mut BRep, tolerance: &Tolerance) -> Vec<SliverAbsorbWarning> {
    let _ = tolerance; // step 3 will use this for Class A/B thresholds
    let cap = brep.faces.len() * 2 + 16;
    for _ in 0..cap {
        let Some(he) = find_next_absorbable(brep) else {
            break;
        };
        if brep.kef(he).is_err() {
            break;
        }
    }
    // Step 3: classify any remaining slivers and return
    // `Retained` warnings. Step 2 returns an empty list because
    // the only handled class (corner fragments) is always
    // absorbable when classified.
    Vec::new()
}

/// Walk the BRep faces looking for the first absorbable sliver.
/// Returns the half-edge whose face should SURVIVE the `kef` merge
/// — the caller passes it directly to `brep.kef(he)`.
///
/// Priority order (step 2):
/// 1. Corner fragment (Class D) — ≤4-vert Plane whose twin face
///    via some edge is a Sphere. Returns the sphere-side half-edge
///    so the sphere labelling survives.
///
/// Step 3 will prepend Classes A (zero-area), B (thin), C
/// (redundant-edge) in that order. The priority order is
/// "corner-fragments first" because they are the only known-
/// clean case with full test coverage.
fn find_next_absorbable(brep: &BRep) -> Option<HalfEdgeHandle> {
    find_corner_fragment_target(brep)
}

/// Class D: `≤4-vert Plane face whose outer loop contains a half-edge
/// whose twin face is a Sphere`. Returns the *sphere-side* half-edge
/// so a subsequent `kef` preserves the sphere labelling.
///
/// Exact port of the previous `merge_corner_fragments_into_blends`
/// logic, repackaged as a step of the general classifier.
fn find_corner_fragment_target(brep: &BRep) -> Option<HalfEdgeHandle> {
    for fh in brep.faces.keys() {
        let face = &brep.faces[fh];
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        if !matches!(sof.geometry, SurfaceGeometry::Plane { .. }) {
            continue;
        }
        let positions = brep.face_outer_loop_positions(fh);
        if positions.len() > 4 {
            continue;
        }
        let start_he = brep.loops[face.outer_loop].half_edge;
        let mut cursor = start_he;
        // Safety cap on the inner walk: a malformed half-edge
        // chain shouldn't be able to spin us forever.
        let walk_cap = brep.half_edges.len() + 1;
        for _ in 0..walk_cap {
            let he = &brep.half_edges[cursor];
            let twin_h = he.twin;
            let twin = &brep.half_edges[twin_h];
            let twin_face = brep.loops[twin.loop_].face;
            if twin_face != fh
                && let Some(twin_sof) = brep.faces[twin_face].surface.as_ref()
                && matches!(twin_sof.geometry, SurfaceGeometry::Sphere { .. })
            {
                return Some(twin_h);
            }
            cursor = he.next;
            if cursor == start_he {
                break;
            }
        }
    }
    None
}
