//! Post-blend structural + corner-resolution gates and fallback-edge
//! attribution helpers
//! for [`super::fillet_edges_attempt`] (extracted for the 500-line
//! structural file budget).

use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

use crate::construct::{ApplyBlendResult, EdgeGeometry};
use crate::error::{BlendError, BlendFailure};
use crate::types::{ApexConstructionRoute, FilletSpec};

/// The transactional served-body gate of one pipeline attempt: the
/// first [`served_body_failure`] as the attempt's decline. When it
/// declines, `gated_direct_apex_corners` receives every corner the
/// attempt sent through direct-apex, served or fallen back on its
/// anchor convention, with the gate's reason. The serving retry
/// re-runs exactly those corners through the vertex-blend dispatch
/// (see [`crate::fillet_serving_retry`]).
pub(super) fn served_body_gate(
    source: &BRep,
    spec_count: usize,
    served: &ApplyBlendResult,
    fallback_edge: EdgeHandle,
    gated_direct_apex_corners: &mut Vec<(VertexHandle, String)>,
) -> Result<(), BlendError> {
    let Some(failure) = served_body_failure(source, served, fallback_edge) else {
        return Ok(());
    };
    *gated_direct_apex_corners = served
        .repair_report
        .apex_routes
        .iter()
        .filter(|route| match route.route {
            ApexConstructionRoute::DirectApex => !route.built_faces.is_empty(),
            ApexConstructionRoute::DirectApexFallbackToVertexBlend => true,
            _ => false,
        })
        .map(|route| {
            let fallback = route
                .fallback_reason
                .as_ref()
                .map(|reason| format!(" (direct-apex fallback: {reason})"))
                .unwrap_or_default();
            (
                route.vertex,
                format!(
                    "served-body gate declined the {:?} serve{fallback}: {:?}",
                    route.route, failure.1
                ),
            )
        })
        .collect();
    Err(BlendError::partial_failure(
        spec_count,
        source.clone(),
        vec![failure],
    ))
}

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
/// unresolved corner whose territory fails check 9, else any
/// inconsistency in the whole blend territory (PK-530).
fn served_body_failure(
    source: &BRep,
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
    corner_resolution_failure(served).or_else(|| {
        crate::construct::territory_consistency::territory_inconsistency(source, &served.brep)
            .map(|failure| (fallback_edge, failure))
    })
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
                | popkern_core::error::ValidationError::LoopJunctionMismatch { .. }
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
