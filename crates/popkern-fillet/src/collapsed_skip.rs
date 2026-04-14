//! §13.1(5) collapsed-edge soft-skip pre-pass.
//!
//! Spec text: *"Collapsed edges (length < tol_geom): treat as
//! vertex blend or reject."* The repair strategy shipped here is the
//! "treat as vertex blend" arm: when a fillet spec targets an edge
//! whose two endpoint vertices are separated by less than
//! `tolerance.linear()`, drop the spec from the batch with a
//! [`BlendWarning::CollapsedEdgeSkipped`] rather than hard-rejecting
//! the whole operation. The pipeline then runs on the healthy
//! spec subset, and the neighbouring edges' vertex-blend phase
//! naturally handles the geometry around the collapsed point.
//!
//! ## Why "soft-skip" matches the spec
//!
//! The commercial CAD interpretation of "treat as vertex blend" on a
//! standalone collapsed edge is: **the user asked to fillet a
//! point-like edge; filleting a point is meaningless, but the point
//! itself is still part of the model's vertex set and gets whatever
//! blending the neighbouring edges provide.** The pipeline already
//! has full vertex-blend machinery
//! ([`crate::types::VertexMode::Smooth`] / `Mitred`) that joins
//! multiple fillets at a shared vertex — the skip path delegates to
//! it implicitly by not competing for the same vertex.
//!
//! Standalone collapsed edges with no filletable neighbours still
//! skip silently (warning, not failure). This is arguably less
//! useful than a hard reject for those cases, but it matches the
//! spec's *"or reject"* phrasing less strictly than the spec's
//! *"treat as vertex blend"* phrasing, which we read as the
//! preferred arm when it applies.
//!
//! ## Ordering
//!
//! The skip runs *before* the cusp-subdivide pre-pass in
//! [`crate::fillet_edges_internal`]: collapsed edges are typically
//! in the original BRep (lossy imports, numerical drift from
//! boolean ops), not created by cusp subdivision, so there is
//! nothing to gain by running in the other order. Running first
//! also means the cusp pre-pass never wastes work scanning a
//! collapsed edge for interior cusps.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BlendWarning;
use crate::types::FilletSpec;

/// Detect and soft-skip collapsed edges in a fillet spec list.
///
/// Returns `(kept_specs, warnings)` where:
///
/// - `kept_specs` is the input spec list with every collapsed-edge
///   spec removed. Length ≤ input length.
/// - `warnings` contains one [`BlendWarning::CollapsedEdgeSkipped`]
///   per dropped spec, reporting the edge handle and measured
///   length.
///
/// **Non-collapsed-edge behaviour:**
/// - A **closed edge** (`v_start == v_end`, e.g. a cylinder seam)
///   is **not** a collapsed edge — it is the intentional
///   closed-rim case, preserved unchanged in `kept_specs`.
/// - An edge whose handle is stale, or whose half-edge/vertex
///   chain is inconsistent, is left in the spec list for the
///   downstream detectors to report the real problem.
///
/// Intended to be called at the top of
/// [`crate::fillet_edges_internal`], immediately after Phase 0
/// validation and before the cusp-subdivide pre-pass.
pub fn skip_collapsed_edges(
    brep: &BRep,
    specs: &[FilletSpec],
    tolerance: &Tolerance,
) -> (Vec<FilletSpec>, Vec<BlendWarning>) {
    let mut kept: Vec<FilletSpec> = Vec::with_capacity(specs.len());
    let mut warnings: Vec<BlendWarning> = Vec::new();
    let tol = tolerance.linear();

    for spec in specs {
        let Some(edge_data) = brep.edges.get(spec.edge) else {
            kept.push(spec.clone());
            continue;
        };
        let Some(he) = brep.half_edges.get(edge_data.half_edge) else {
            kept.push(spec.clone());
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            kept.push(spec.clone());
            continue;
        };
        let v_start = he.origin;
        let v_end = twin.origin;

        // Intentional closed edges — skip the length check.
        if v_start == v_end {
            kept.push(spec.clone());
            continue;
        }

        let Some(start_v) = brep.vertices.get(v_start) else {
            kept.push(spec.clone());
            continue;
        };
        let Some(end_v) = brep.vertices.get(v_end) else {
            kept.push(spec.clone());
            continue;
        };
        let length = (end_v.position - start_v.position).length();
        if length < tol {
            warnings.push(BlendWarning::CollapsedEdgeSkipped {
                edge: spec.edge,
                length,
            });
        } else {
            kept.push(spec.clone());
        }
    }

    (kept, warnings)
}

#[cfg(test)]
#[path = "collapsed_skip_tests.rs"]
mod tests;
