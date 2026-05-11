//! Boundary-side separator bridge preflight.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::euler::{
    BridgeDiagnosticOperationRequest, BridgeDiagnosticReport, BridgeDiagnosticRequest,
    DiagnosticHandle, DiagnosticOperationPurpose, DiagnosticResourceKind,
    FaceFamilyBridgeCycleRequest, ResourceConflictKind,
};
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle};

/// Aggregate bridge request requiring one already-boundary separator side.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundarySideSeparatorBridgeRequest {
    /// Candidate bridge cycles in semantic request order.
    pub cycles: Vec<FaceFamilyBridgeCycleRequest>,
}

/// Validated boundary-side separator bridge plan.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundarySideSeparatorBridgePlan {
    /// Candidate bridge cycles in semantic request order.
    pub cycles: Vec<FaceFamilyBridgeCycleRequest>,
    /// Boundary edge that needs one additional allocation.
    pub separator_source_edge: EdgeHandle,
    /// Boundary half-edge side that is currently selected twice.
    pub separator_source_half_edge: HalfEdgeHandle,
    /// Expected topology delta for separator allocation plus bridge faces.
    pub expected_delta: BoundarySideSeparatorBridgeTopologyDelta,
}

/// Expected topology count delta for a boundary-side separator bridge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BoundarySideSeparatorBridgeTopologyDelta {
    /// Vertex count delta.
    pub vertices: i64,
    /// Edge count delta.
    pub edges: i64,
    /// Face count delta.
    pub faces: i64,
    /// Loop count delta.
    pub loops: i64,
}

impl BRep {
    /// Validate an aggregate bridge shape that needs one boundary-side separator.
    pub fn preflight_boundary_side_separator_bridge(
        &self,
        request: BoundarySideSeparatorBridgeRequest,
    ) -> Result<BoundarySideSeparatorBridgePlan, CoreError> {
        if request.cycles.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "boundary-side separator bridge requires at least one cycle".into(),
            ));
        }
        let report = self.bridge_diagnostic_report(BridgeDiagnosticRequest {
            operations: request
                .cycles
                .iter()
                .map(|cycle| BridgeDiagnosticOperationRequest {
                    purpose: DiagnosticOperationPurpose::BridgeCycle,
                    first_loop: cycle.first_loop,
                    first_interval: cycle.first_interval.clone(),
                    second_loop: cycle.second_loop,
                    second_interval: cycle.second_interval.clone(),
                })
                .collect(),
        })?;
        require_boundary_side_separator_candidate_gate(&report)?;
        let summary = report.summary;
        if summary.already_boundary_overlaps != 1
            || summary.duplicate_boundary_edge_consumptions != 1
            || summary.duplicate_half_edge_consumptions != 1
            || summary.missing_operators != 1
        {
            return Err(CoreError::PreconditionFailed(format!(
                "boundary-side separator bridge requires exactly one already-boundary duplicate side; {}",
                summary.diagnostic_label()
            )));
        }
        let separator_source_edge =
            unique_conflict_edge(&report.graph, ResourceConflictKind::AlreadyBoundaryOverlap)?;
        let duplicate_edge = unique_conflict_edge(
            &report.graph,
            ResourceConflictKind::DuplicateBoundaryEdgeConsumption,
        )?;
        if separator_source_edge != duplicate_edge {
            return Err(CoreError::PreconditionFailed(
                "boundary-side separator conflict edge does not match duplicate boundary edge"
                    .into(),
            ));
        }
        let separator_source_half_edge = unique_conflict_half_edge(
            &report.graph,
            ResourceConflictKind::DuplicateHalfEdgeConsumption,
        )?;
        let half_edge = self
            .half_edges
            .get(separator_source_half_edge)
            .ok_or_else(|| {
                CoreError::InvalidHandle(
                    "boundary-side separator duplicate half-edge is stale".into(),
                )
            })?;
        if half_edge.edge != separator_source_edge {
            return Err(CoreError::PreconditionFailed(
                "boundary-side separator duplicate half-edge does not own separator edge".into(),
            ));
        }
        if !self.boundary_edges.contains(&separator_source_edge)
            || self.half_edges.contains_key(half_edge.twin)
        {
            return Err(CoreError::PreconditionFailed(
                "boundary-side separator source edge is not an unowned boundary side".into(),
            ));
        }
        require_one_same_loop_already_boundary_overlap(
            self,
            &request.cycles,
            separator_source_edge,
        )?;
        let cycle_count = request.cycles.len() as i64;
        Ok(BoundarySideSeparatorBridgePlan {
            cycles: request.cycles,
            separator_source_edge,
            separator_source_half_edge,
            expected_delta: BoundarySideSeparatorBridgeTopologyDelta {
                vertices: 0,
                edges: 1 + cycle_count * 2,
                faces: cycle_count,
                loops: cycle_count,
            },
        })
    }
}

fn require_boundary_side_separator_candidate_gate(
    report: &BridgeDiagnosticReport,
) -> Result<(), CoreError> {
    let Some(candidate) = report.graph.candidates.iter().find(|candidate| {
        candidate
            .feasibility_reason
            .contains("boundary-side-separator")
    }) else {
        return Err(CoreError::PreconditionFailed(
            "boundary-side separator bridge requires a typed construction candidate".into(),
        ));
    };
    candidate.validate_mutation_gate().map_err(|err| {
        CoreError::PreconditionFailed(format!(
            "boundary-side separator bridge candidate gate rejected {}",
            err.missing_facet.diagnostic_label()
        ))
    })?;
    Ok(())
}

fn unique_conflict_edge(
    graph: &crate::euler::TopologyDiagnosticGraph,
    kind: ResourceConflictKind,
) -> Result<EdgeHandle, CoreError> {
    let edges = graph
        .conflicts
        .iter()
        .filter(|conflict| conflict.kind == kind)
        .flat_map(|conflict| conflict.resources.iter())
        .filter_map(|resource_id| graph.resources.get(resource_id.0))
        .filter(|resource| resource.kind == DiagnosticResourceKind::EdgeIdentity)
        .filter_map(|resource| match resource.handle {
            DiagnosticHandle::Edge(edge) => Some(edge),
            _ => None,
        })
        .collect::<BTreeSet<_>>();
    if edges.len() == 1 {
        return edges.iter().next().copied().ok_or_else(|| {
            CoreError::PreconditionFailed(
                "boundary-side separator expected one edge but none was retained".into(),
            )
        });
    }
    Err(CoreError::PreconditionFailed(format!(
        "boundary-side separator expected one {kind:?} edge, got {}",
        edges.len()
    )))
}

fn unique_conflict_half_edge(
    graph: &crate::euler::TopologyDiagnosticGraph,
    kind: ResourceConflictKind,
) -> Result<HalfEdgeHandle, CoreError> {
    let half_edges = graph
        .conflicts
        .iter()
        .filter(|conflict| conflict.kind == kind)
        .flat_map(|conflict| conflict.resources.iter())
        .filter_map(|resource_id| graph.resources.get(resource_id.0))
        .filter(|resource| resource.kind == DiagnosticResourceKind::HalfEdgeSide)
        .filter_map(|resource| match resource.handle {
            DiagnosticHandle::HalfEdge(half_edge) => Some(half_edge),
            _ => None,
        })
        .collect::<BTreeSet<_>>();
    if half_edges.len() == 1 {
        return half_edges.iter().next().copied().ok_or_else(|| {
            CoreError::PreconditionFailed(
                "boundary-side separator expected one half-edge but none was retained".into(),
            )
        });
    }
    Err(CoreError::PreconditionFailed(format!(
        "boundary-side separator expected one {kind:?} half-edge, got {}",
        half_edges.len()
    )))
}

fn require_one_same_loop_already_boundary_overlap(
    brep: &BRep,
    cycles: &[FaceFamilyBridgeCycleRequest],
    separator_edge: EdgeHandle,
) -> Result<(), CoreError> {
    let mut matching_same_loop_overlaps = 0;
    for cycle in cycles {
        if cycle.first_loop != cycle.second_loop {
            continue;
        }
        let first_edges = interval_edges(brep, &cycle.first_interval)?;
        let second_edges = interval_edges(brep, &cycle.second_interval)?;
        let overlap = first_edges
            .intersection(&second_edges)
            .copied()
            .collect::<BTreeSet<_>>();
        if overlap.contains(&separator_edge) {
            matching_same_loop_overlaps += 1;
        }
    }
    if matching_same_loop_overlaps == 1 {
        return Ok(());
    }
    Err(CoreError::PreconditionFailed(format!(
        "boundary-side separator expected one same-loop overlap for separator edge, got {matching_same_loop_overlaps}"
    )))
}

fn interval_edges(
    brep: &BRep,
    interval: &[HalfEdgeHandle],
) -> Result<BTreeSet<EdgeHandle>, CoreError> {
    if interval.len() < 2 {
        return Err(CoreError::PreconditionFailed(
            "boundary-side separator interval must include an after-last half-edge".into(),
        ));
    }
    interval[..interval.len() - 1]
        .iter()
        .map(|half_edge| {
            brep.half_edges
                .get(*half_edge)
                .map(|data| data.edge)
                .ok_or_else(|| {
                    CoreError::InvalidHandle(
                        "boundary-side separator interval has stale half-edge".into(),
                    )
                })
        })
        .collect()
}

#[cfg(test)]
#[path = "boundary_side_separator_bridge_tests.rs"]
mod tests;
