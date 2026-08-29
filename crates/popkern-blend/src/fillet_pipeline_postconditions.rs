//! Post-blend structural gate + fallback-edge attribution helpers
//! for [`super::fillet_edges_attempt`] (extracted for the 500-line
//! structural file budget).

use popkern_core::topology::{BRep, EdgeHandle};

use crate::construct::EdgeGeometry;
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

/// First hard structural defect on the served body, if any — the
/// post-blend transactional gate (broken cycles/twins, inconsistent
/// vertices, non-manifold edges, open loops).
pub(super) fn first_structural_validation_error(brep: &BRep) -> Option<String> {
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
