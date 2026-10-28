//! Chamfer declines at mixed corners (PK-994).
//!
//! A vertex where a CONCAVE chamfer (a fill) and a CONVEX chamfer (a
//! cut) both end has no construction: the exact corner depends on how
//! the added and removed wedges are meant to meet, and the per-edge
//! pipeline's residue there fails its own serve gates (measured: an
//! L-prism's re-entrant vertical chamfered with a top edge, a pocket
//! wall's vertical with its rim). Chamfering the concave and the convex
//! edges in separate operations is well defined (each sees the other's
//! result) and serves.
//!
//! When a chamfer request that contains such a corner declines, the
//! per-edge failures at that corner are restated as
//! [`BlendFailure::VertexBlendFailed`] naming the corner, with the
//! original reason kept, so the remediation points at the split instead
//! of at the input body. Requests without a mixed corner, and failures
//! away from one, pass through unchanged.

use crate::analysis::{Convexity, analyze_edges};
use crate::error::{BlendError, BlendFailure};
use crate::types::ChamferSpec;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};
use std::collections::BTreeMap;

/// The reason prefix the remediation mapping recognises.
pub(crate) const MIXED_CORNER: &str = "mixed chamfer corner";

/// Restate a declined chamfer's failures at mixed corners.
pub(crate) fn explain_mixed_corners(
    brep: &BRep,
    specs: &[ChamferSpec],
    error: BlendError,
) -> BlendError {
    let BlendError::PartialFailure {
        total,
        partial_result,
        failures,
        ..
    } = error
    else {
        return error;
    };
    let mixed = mixed_corners(brep, specs);
    if mixed.is_empty() {
        return BlendError::partial_failure(total, *partial_result, failures);
    }
    let failures = failures
        .into_iter()
        .map(
            |(edge, failure)| match mixed.iter().find(|(_, edges)| edges.contains(&edge)) {
                Some((vertex, _)) => (edge, restate(brep, *vertex, failure)),
                None => (edge, failure),
            },
        )
        .collect();
    BlendError::partial_failure(total, *partial_result, failures)
}

fn restate(brep: &BRep, vertex: VertexHandle, failure: BlendFailure) -> BlendFailure {
    let location = brep.vertices.get(vertex).map(|v| v.position);
    BlendFailure::VertexBlendFailed {
        vertex,
        reason: format!(
            "{MIXED_CORNER} (construction declined: {})",
            failure_reason(&failure)
        ),
        location,
    }
}

fn failure_reason(failure: &BlendFailure) -> String {
    match failure {
        BlendFailure::InvalidTopology { reason, .. }
        | BlendFailure::InconsistentGeometry { reason, .. }
        | BlendFailure::VertexBlendFailed { reason, .. }
        | BlendFailure::IntersectionFailed { reason, .. }
        | BlendFailure::OffsetFailed { reason, .. } => reason.clone(),
        other => format!("{other:?}"),
    }
}

/// Vertices where requested chamfers of both convexities end, with the
/// requested edges at each.
fn mixed_corners(brep: &BRep, specs: &[ChamferSpec]) -> Vec<(VertexHandle, Vec<EdgeHandle>)> {
    let edges: Vec<EdgeHandle> = specs.iter().map(|s| s.edge).collect();
    let Ok(analyses) = analyze_edges(brep, &edges, false) else {
        return Vec::new();
    };
    // vertex -> (edges, has convex, has concave)
    let mut at: BTreeMap<VertexHandle, (Vec<EdgeHandle>, bool, bool)> = BTreeMap::new();
    for a in &analyses {
        for v in [a.vertex_start, a.vertex_end] {
            let entry = at.entry(v).or_default();
            entry.0.push(a.edge);
            entry.1 |= a.convexity == Convexity::Convex;
            entry.2 |= a.convexity == Convexity::Concave;
        }
    }
    at.into_iter()
        .filter(|(_, (_, convex, concave))| *convex && *concave)
        .map(|(v, (edges, _, _))| (v, edges))
        .collect()
}
