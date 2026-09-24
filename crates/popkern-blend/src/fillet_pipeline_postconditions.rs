//! Post-blend structural + corner-resolution gates and fallback-edge
//! attribution helpers
//! for [`super::fillet_edges_attempt`] (extracted for the 500-line
//! structural file budget).

use popkern_core::topology::{BRep, EdgeHandle};

use crate::construct::{ApplyBlendResult, EdgeGeometry};
use crate::error::BlendFailure;
use crate::types::FilletSpec;

/// The edge failures are attributed to when no more specific edge is
/// known: the first spec's edge, else the first computed geometry's.
pub(super) fn fallback_edge_for_specs(
    specs: &[FilletSpec],
    edge_geometries: &[EdgeGeometry],
) -> EdgeHandle {
    specs
        .first()
        .map(|s| s.edge)
        .or_else(|| edge_geometries.first().map(|g| g.analysis.edge))
        .unwrap_or_default()
}

/// Append the overlap and §13.1 pre-pass warnings to the pipeline's
/// warning list: collapsed-edge skips name their own edge, the rest
/// attribute to `fallback_edge`.
pub(super) fn attribute_pipeline_warnings(
    all_warnings: &mut Vec<(EdgeHandle, crate::error::BlendWarning)>,
    fallback_edge: EdgeHandle,
    overlap_warnings: Vec<crate::error::BlendWarning>,
    prepass_warnings: Vec<crate::error::BlendWarning>,
) {
    for warn in overlap_warnings {
        all_warnings.push((fallback_edge, warn));
    }
    for warn in prepass_warnings {
        let edge = match &warn {
            crate::error::BlendWarning::CollapsedEdgeSkipped { edge, .. } => *edge,
            _ => fallback_edge,
        };
        all_warnings.push((edge, warn));
    }
}

/// The served body's first transactional-gate failure, if any: a hard
/// structural defect (attributed to `fallback_edge`), else an
/// unresolved corner whose territory fails check 9.
pub(super) fn served_body_failure(
    served: &ApplyBlendResult,
    fallback_edge: EdgeHandle,
) -> Option<(EdgeHandle, BlendFailure)> {
    if let Some(reason) = first_structural_validation_error(&served.brep) {
        return Some((
            fallback_edge,
            BlendFailure::InvalidTopology {
                reason,
                location: None,
            },
        ));
    }
    corner_resolution_failure(served)
}

/// First hard structural defect on the served body, if any — the
/// post-blend transactional gate (broken cycles/twins, inconsistent
/// vertices, non-manifold edges, open loops).
fn first_structural_validation_error(brep: &BRep) -> Option<String> {
    [
        brep.validate_half_edge_pairing().err(),
        brep.validate_half_edge_cycles().err(),
        brep.validate_vertex_consistency().err(),
        brep.validate_manifold_edges().err(),
        brep.validate_loop_closure().err(),
    ]
    .into_iter()
    .flatten()
    .find(|err| {
        matches!(
            err,
            popkern_core::error::ValidationError::BrokenCycle { .. }
                | popkern_core::error::ValidationError::BrokenTwin { .. }
                | popkern_core::error::ValidationError::InconsistentVertex { .. }
                | popkern_core::error::ValidationError::NonManifoldEdge { .. }
                | popkern_core::error::ValidationError::OpenLoop { .. }
        )
    })
    .map(|err| format!("post-blend structural validation failed: {err:?}"))
}

/// PK-516 corner-resolution gate on the served body: the first corner
/// left without a corner blend whose territory fails check 9, as the
/// transactional `VertexBlendFailed` decline (see
/// [`crate::construct::corner_resolution`]). Runs after every tail
/// pass, so corners a later pass closes exactly serve.
fn corner_resolution_failure(served: &ApplyBlendResult) -> Option<(EdgeHandle, BlendFailure)> {
    let violation = crate::construct::corner_resolution::first_corner_resolution_violation(
        &served.brep,
        &served.unresolved_corners,
    )?;
    Some((
        violation.corner.edge,
        BlendFailure::VertexBlendFailed {
            vertex: violation.corner.vertex,
            reason: violation.reason(),
            location: Some(violation.corner.position),
        },
    ))
}
