//! General post-construction sliver-face absorption pass.
//!
//! **Step 2 scope**: infrastructure only. This module subsumes the
//! previous specific-case `blend_cleanup::merge_corner_fragments_into_blends`
//! absorber (Class D — "≤4-vert Plane → Sphere/Gregory vertex blend") under
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

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

// Classifier helpers live in a sibling module to keep this file
// under the 500-line structural cap. Re-imported here so the
// `#[path]`-mounted test module can still reach them via
// `super::*`.
use crate::construct::sliver_absorb_classifiers::{
    find_corner_fragment_target, find_redundant_edge_target, find_short_edge_target,
    find_zero_area_target,
};
// `surfaces_equal` is consumed transitively via
// `find_redundant_edge_target` in non-test builds, but the
// `#[path]`-mounted test module reaches it through `super::*`.
#[cfg(test)]
#[allow(unused_imports)]
use crate::construct::sliver_absorb_classifiers::surfaces_equal;

/// Run the sliver absorber and route any retained-sliver warnings
/// into the caller's `BlendWarning` vector. Kept here rather than
/// inline in `construct::mod` so the dispatch logic stays colocated
/// with the classifier and the parent file stays under the 500-line
/// structural cap.
#[allow(dead_code)]
pub(crate) fn run_and_route_warnings(
    brep: &mut BRep,
    tolerance: &Tolerance,
    edge_geometries: &[crate::construct::EdgeGeometry],
    warnings: &mut Vec<(
        popkern_core::topology::EdgeHandle,
        crate::error::BlendWarning,
    )>,
) {
    run_and_route_warnings_with_protected(brep, tolerance, edge_geometries, &[], warnings)
}

/// As [`run_and_route_warnings`], plus an opt-in set of cap faces
/// that the absorber must NOT touch.  Used by Phase 6c-D Stage 4.5:
/// when sphere-fit boundary substitution + bi-9 K-P injection have
/// fired, each cap face's outer loop is the canonical 3-edge
/// per-fillet topology that the injected `hex_boundaries[k]`
/// describes.  Letting `sliver_absorb` merge adjacent Plane
/// fragments into those caps expands the cap-face outer loop to
/// span 2-3 fillets — the Gregory cap surface (row 0 = ONE
/// FilletEnd) cannot reproduce that boundary, and `cap_outer_g0`
/// stalls at 1-6 cm regardless of how well the L-M solver fits the
/// interior.  Skipping absorption for these caps preserves the
/// per-fillet hex topology.
pub(crate) fn run_and_route_warnings_with_protected(
    brep: &mut BRep,
    tolerance: &Tolerance,
    edge_geometries: &[crate::construct::EdgeGeometry],
    protected: &[FaceHandle],
    warnings: &mut Vec<(
        popkern_core::topology::EdgeHandle,
        crate::error::BlendWarning,
    )>,
) {
    let protected_set: std::collections::BTreeSet<FaceHandle> = protected.iter().copied().collect();
    let sliver_warnings = absorb_slivers_with_protected(brep, tolerance, &protected_set);
    let fallback = edge_geometries
        .first()
        .map(|eg| eg.analysis.edge)
        .unwrap_or_default();
    for SliverAbsorbWarning::Retained { face, kind } in sliver_warnings {
        warnings.push((
            fallback,
            crate::error::BlendWarning::SliverRetained {
                face,
                kind: kind.label(),
            },
        ));
    }
}

#[cfg(test)]
#[path = "sliver_absorb_tests.rs"]
mod tests;

/// A classified sliver face that the absorber could not remove.
///
/// Returned by [`absorb_slivers`] for any face that a classifier
/// flagged but the corresponding absorber refused to merge — typically
/// because no suitable absorbable neighbour exists or the merge would
/// violate topology. Callers consume these via `BlendWarning` plumbing:
/// `apply_blend_geometry` converts them into
/// [`crate::error::BlendWarning::SliverRetained`] entries on the
/// per-edge warning list.
#[derive(Debug, Clone)]
pub(crate) enum SliverAbsorbWarning {
    /// A sliver was classified but could not be absorbed.
    Retained {
        /// The face that remains in the output BRep.
        face: FaceHandle,
        /// Which sliver category the classifier flagged.
        kind: SliverKind,
    },
}

/// Classification of a sliver face. Each variant carries the
/// measurement that drove the classification, for diagnostic
/// rendering via `Debug` when the sliver ends up retained. The
/// inner fields are not read by non-`Debug` code paths today —
/// the `#[allow(dead_code)]` is intentional.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub(crate) enum SliverKind {
    /// Polygonal outer-loop area below `tolerance.linear()²`.
    ZeroArea {
        /// The measured polygonal area.
        area: f64,
    },
    /// Face with an outer-loop edge shorter than `tolerance.linear()`
    /// (a "thin strip" produced by near-parallel trim curves).
    ShortEdge {
        /// Length of the shortest edge.
        min_length: f64,
    },
    /// Face shares an edge with a neighbour whose underlying
    /// surface is structurally equal (duplicate-subdivision artifact).
    RedundantEdge,
    /// ≤4-vert face whose twin face via some edge is a vertex-
    /// blend sphere — the Round 5r corner-fragment case.
    CornerFragment {
        /// Number of vertices on the outer loop (≤ 4 when flagged).
        vertex_count: usize,
    },
}

impl SliverKind {
    /// Short static-string label for the `BlendWarning::SliverRetained`
    /// variant's `kind` field. Callers match on this string rather
    /// than on a cloned enum to keep the public `BlendWarning` small.
    pub(crate) fn label(&self) -> &'static str {
        match self {
            SliverKind::ZeroArea { .. } => "zero-area",
            SliverKind::ShortEdge { .. } => "short-edge",
            SliverKind::RedundantEdge => "redundant-edge",
            SliverKind::CornerFragment { .. } => "corner-fragment",
        }
    }
}

/// Run the sliver absorption pass. Classify → apply one `kef` →
/// reclassify, capped at `faces.len() * 2 + 16` iterations.
/// Priority: Class D (corner fragment → sphere, test-covered) →
/// A (zero-area) → C (redundant edge) → B (short edge, roughest
/// heuristic, last). Failed kefs mark the face as retained so the
/// pass moves on rather than spinning.
#[allow(dead_code)]
pub(crate) fn absorb_slivers(brep: &mut BRep, tolerance: &Tolerance) -> Vec<SliverAbsorbWarning> {
    absorb_slivers_with_protected(brep, tolerance, &std::collections::BTreeSet::new())
}

/// Variant of [`absorb_slivers`] that refuses any kef whose surviving
/// face is in `protected`.  Used by Phase 6c-D Stage 4.5 to keep
/// substituted-boundary cap faces from absorbing adjacent Plane
/// fragments and growing across multiple fillets' boundaries.
fn absorb_slivers_with_protected(
    brep: &mut BRep,
    tolerance: &Tolerance,
    protected: &std::collections::BTreeSet<FaceHandle>,
) -> Vec<SliverAbsorbWarning> {
    let mut warnings: Vec<SliverAbsorbWarning> = Vec::new();
    let cap = brep.faces.len() * 2 + 16;
    let mut retained: std::collections::HashSet<FaceHandle> = std::collections::HashSet::new();
    for _ in 0..cap {
        let candidate = find_corner_fragment_target(brep, &retained)
            .map(|(face, he)| {
                let verts = brep.face_outer_loop_positions(face).len();
                (
                    face,
                    SliverKind::CornerFragment {
                        vertex_count: verts,
                    },
                    he,
                )
            })
            .or_else(|| {
                find_zero_area_target(brep, tolerance, &retained)
                    .map(|(face, he, area)| (face, SliverKind::ZeroArea { area }, he))
            })
            .or_else(|| {
                find_redundant_edge_target(brep, tolerance, &retained)
                    .map(|(face, he)| (face, SliverKind::RedundantEdge, he))
            })
            .or_else(|| {
                find_short_edge_target(brep, tolerance, &retained)
                    .map(|(face, he, min_length)| (face, SliverKind::ShortEdge { min_length }, he))
            });
        let Some((face, kind, he)) = candidate else {
            break;
        };
        // Phase 6c-D Stage 4.5: if the surviving face (the one
        // containing `he`) is a substituted-boundary cap, skip the
        // kef.  The caller has injected hex-boundary CPs that
        // describe the cap's CURRENT outer loop; absorbing a
        // neighbour would expand the loop across multiple fillets'
        // boundaries and the Gregory cap surface (row 0 = one
        // FilletEnd) can no longer represent it.  Mark the sliver
        // as retained so the loop moves on rather than retrying.
        let surviving_face = brep.half_edges.get(he).map(|h| brep.loops[h.loop_].face);
        if let Some(sf) = surviving_face
            && protected.contains(&sf)
        {
            warnings.push(SliverAbsorbWarning::Retained { face, kind });
            retained.insert(face);
            continue;
        }
        if brep.kef(he).is_err() {
            warnings.push(SliverAbsorbWarning::Retained { face, kind });
            retained.insert(face);
        }
    }
    warnings
}
