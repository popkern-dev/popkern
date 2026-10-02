//! The corner-cap preference in the apply phase (PK-603): which
//! corners a stated setback / fitted preference keeps from direct-apex
//! promotion, and how a strict decline becomes a fillet failure.
//! Design: `docs/architecture/fillet-corner-cap-preference.md`.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};
use std::collections::BTreeSet;

use super::EdgeGeometry;
use crate::error::BlendFailure;
use crate::types::{ApexConstructionDiagnostic, CornerCapPreference, CornerSetback, VertexMode};

/// The corners the setback rung is eligible for when the preference
/// names a rung construction, or a user setback (PK-623) applies there:
/// a uniform distance claims every eligible corner, per-end overrides
/// only the corners they name. Direct apex must not claim them: it runs
/// before the rung and would pre-empt the stated preference, and its
/// strip ends (the mitre complex) are not a border the bridged carve can
/// cut. Corners outside the rung's reach keep their direct-apex route.
/// Eligibility is the rung's own gate (`setback_min_n_for`, classified
/// on the pre-modification body as the dispatch loop does).
pub(super) fn rung_claimed_corners(
    brep: &BRep,
    vertex_mode: &VertexMode,
    edge_geometries: &[EdgeGeometry],
    (corner_cap, setback): (CornerCapPreference, &CornerSetback),
    tolerance: &Tolerance,
    explicit_subpatch_family: bool,
) -> BTreeSet<VertexHandle> {
    let every_corner = corner_cap.names_setback_rung() || setback.distance.is_some();
    if !(every_corner || !setback.ends.is_empty())
        || *vertex_mode != VertexMode::Smooth
        || edge_geometries.len() < 2
        || !super::setback_vertex_topology::setback_topology_commit_enabled()
    {
        return BTreeSet::new();
    }
    crate::analysis::vertex_incidences(
        edge_geometries
            .iter()
            .map(|eg| (eg.analysis.vertex_start, eg.analysis.vertex_end)),
    )
    .into_iter()
    .filter(|(vertex, _)| every_corner || setback.names_vertex(*vertex))
    .filter(|(vertex, edge_indices)| {
        // A user setback claims any corner of three or more rounds.
        if setback.distance.is_some() || setback.names_vertex(*vertex) {
            return edge_indices.len() >= 3;
        }
        let edges: Vec<_> = edge_indices
            .iter()
            .filter_map(|&i| edge_geometries.get(i).map(|g| g.analysis.edge))
            .collect();
        let classification =
            crate::vertex_classifier::classify_vertex_blend(brep, *vertex, &edges, tolerance)
                .classification;
        let fillet_surfaces =
            super::blend_topo::build_fillet_surface_list(edge_indices, edge_geometries);
        edge_indices.len()
            >= super::blend_insert::setback_min_n_for(
                edge_indices,
                edge_geometries,
                &fillet_surfaces,
                classification,
                tolerance,
                explicit_subpatch_family,
            )
    })
    .map(|(vertex, _)| vertex)
    .collect()
}

/// The failure a strict preference's decline at `vertex` becomes, if
/// the rung reported one (a user setback's refusals included, PK-623).
pub(super) fn strict_decline_failure(
    brep: &BRep,
    diagnostics: &[ApexConstructionDiagnostic],
    edge: EdgeHandle,
) -> Option<(EdgeHandle, BlendFailure)> {
    diagnostics.iter().find_map(|diagnostic| match diagnostic {
        ApexConstructionDiagnostic::CornerCapUnavailable {
            vertex,
            construction,
            reason,
        } => Some((
            edge,
            BlendFailure::CornerCapUnavailable {
                vertex: *vertex,
                construction: *construction,
                reason: reason.clone(),
                location: brep.vertices.get(*vertex).map(|v| v.position),
            },
        )),
        // PK-623: a user setback outside the corner's range.
        ApexConstructionDiagnostic::SetbackOutOfRange {
            vertex,
            edge,
            requested,
            minimum,
            maximum,
        } => Some((
            *edge,
            BlendFailure::SetbackOutOfRange {
                vertex: *vertex,
                edge: *edge,
                requested: *requested,
                minimum: *minimum,
                maximum: *maximum,
                location: brep.vertices.get(*vertex).map(|v| v.position),
            },
        )),
        _ => None,
    })
}
