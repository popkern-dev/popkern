use popkern_core::topology::{FaceHandle, VertexHandle};

use super::EdgeGeometry;

pub(super) struct ApexRouteDecisionInput<'a> {
    pub vertex: VertexHandle,
    pub edge_indices: &'a [usize],
    pub edge_geometries: &'a [EdgeGeometry],
    pub route: crate::types::ApexConstructionRoute,
    pub built_faces: Vec<FaceHandle>,
    pub corner_patch_count: usize,
    pub exact_boundary_sector_count: usize,
    pub exact_boundary_prepared_curve_count: usize,
    pub exact_boundary_inner_spoke_count: usize,
    pub exact_boundary_ribbon_count: usize,
    pub exact_boundary_solved_surface_count: usize,
    pub exact_boundary_pcurve_pair_count: usize,
    pub exact_boundary_performance: crate::types::ExactBoundaryCandidatePerformanceReport,
    pub exact_boundary_report: crate::types::ExactBoundaryCapReport,
    pub fallback_reason: Option<String>,
}

pub(super) fn apex_route_decision(
    input: ApexRouteDecisionInput<'_>,
) -> crate::types::ApexConstructionRouteDecision {
    crate::types::ApexConstructionRouteDecision {
        vertex: input.vertex,
        incident_edges: input
            .edge_indices
            .iter()
            .filter_map(|&index| {
                input
                    .edge_geometries
                    .get(index)
                    .map(|geometry| geometry.analysis.edge)
            })
            .collect(),
        route: input.route,
        built_faces: input.built_faces,
        corner_patch_count: input.corner_patch_count,
        exact_boundary_sector_count: input.exact_boundary_sector_count,
        exact_boundary_prepared_curve_count: input.exact_boundary_prepared_curve_count,
        exact_boundary_inner_spoke_count: input.exact_boundary_inner_spoke_count,
        exact_boundary_ribbon_count: input.exact_boundary_ribbon_count,
        exact_boundary_solved_surface_count: input.exact_boundary_solved_surface_count,
        exact_boundary_pcurve_pair_count: input.exact_boundary_pcurve_pair_count,
        exact_boundary_performance: input.exact_boundary_performance,
        exact_boundary_report: input.exact_boundary_report,
        fallback_reason: input.fallback_reason,
    }
}
