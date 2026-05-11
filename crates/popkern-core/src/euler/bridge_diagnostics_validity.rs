//! Shape validity and repair-readiness facts for bridge diagnostics.

use std::collections::BTreeSet;

use crate::euler::{
    BridgeDiagnosticOperationRequest, ConstructionFeasibility, ValidityRepairFact,
    ValidityRepairReadiness, ValidityRepairStatus,
};
use crate::geometry::surface::SurfaceEval;
use crate::topology::{BRep, HalfEdgeHandle, LoopHandle};

pub(super) fn validity_repair_facts_for_candidate(
    brep: &BRep,
    operation: &BridgeDiagnosticOperationRequest,
    feasibility: ConstructionFeasibility,
) -> Vec<ValidityRepairFact> {
    if let Some(fact) = structural_repair_fact(brep, operation) {
        return vec![fact];
    }
    if let Some(fact) = small_feature_repair_fact(brep, operation) {
        return vec![fact];
    }
    if let Some(fact) = trim_readiness_fact(brep, operation) {
        return vec![fact];
    }
    vec![fact(
        readiness_for_feasibility(feasibility),
        reason_for_feasibility(feasibility),
        status_for_feasibility(feasibility),
        false,
        false,
    )]
}

fn trim_readiness_fact(
    brep: &BRep,
    operation: &BridgeDiagnosticOperationRequest,
) -> Option<ValidityRepairFact> {
    for half_edge in operation_half_edges(operation) {
        let data = brep.half_edges.get(half_edge)?;
        let loop_data = brep.loops.get(data.loop_)?;
        let face = brep.faces.get(loop_data.face)?;
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        let Some(pcurve) = data.pcurve.as_ref() else {
            return Some(fact(
                ValidityRepairReadiness::ValidAfterNonMutatingNormalization,
                "bridge-half-edge-pcurve-missing-recomputable",
                ValidityRepairStatus::NonMutatingNormalizationAvailable,
                false,
                false,
            ));
        };
        let tolerance = brep
            .face_tolerance(loop_data.face)
            .max(brep.edge_tolerance(data.edge));
        let (t0, t1) = pcurve.parameter_range();
        let start_uv = pcurve.point_at(t0);
        let end_uv = pcurve.point_at(t1);
        let start_surface = surface.geometry.point_at(start_uv.x, start_uv.y);
        let end_surface = surface.geometry.point_at(end_uv.x, end_uv.y);
        let origin = brep.vertices.get(data.origin)?.position;
        let next = brep.half_edges.get(data.next)?;
        let destination = brep.vertices.get(next.origin)?.position;
        if start_surface.distance_to(&origin) > tolerance
            || end_surface.distance_to(&destination) > tolerance
        {
            return Some(fact(
                ValidityRepairReadiness::RequiresGeometryRepair,
                "bridge-pcurve-endpoint-same-parameter-deviation",
                ValidityRepairStatus::RepairRequired,
                false,
                false,
            ));
        }
        if brep.measure_edge_deviation(data.edge) > tolerance {
            return Some(fact(
                ValidityRepairReadiness::RequiresGeometryRepair,
                "bridge-edge-pcurve-same-parameter-deviation",
                ValidityRepairStatus::RepairRequired,
                false,
                false,
            ));
        }
    }
    None
}

fn structural_repair_fact(
    brep: &BRep,
    operation: &BridgeDiagnosticOperationRequest,
) -> Option<ValidityRepairFact> {
    for loop_ in [operation.first_loop, operation.second_loop] {
        if !brep.loops.contains_key(loop_) {
            return Some(requires_topology_repair("bridge-source-loop-missing"));
        }
        let face = brep.loops[loop_].face;
        let Some(face_data) = brep.faces.get(face) else {
            return Some(requires_topology_repair("bridge-source-face-missing"));
        };
        if face_data.outer_loop == LoopHandle::default()
            || !brep.loops.contains_key(face_data.outer_loop)
        {
            return Some(requires_topology_repair(
                "bridge-source-face-boundary-missing",
            ));
        }
        if !brep.shells.contains_key(face_data.shell) {
            return Some(requires_topology_repair("bridge-source-shell-missing"));
        }
        if !loop_walk_closes(brep, loop_) {
            return Some(requires_topology_repair(
                "bridge-source-loop-connectivity-invalid",
            ));
        }
    }
    for interval in [&operation.first_interval, &operation.second_interval] {
        if directed_interval_invalid(brep, interval) {
            return Some(requires_topology_repair(
                "bridge-interval-wire-connectivity-invalid",
            ));
        }
    }
    None
}

fn small_feature_repair_fact(
    brep: &BRep,
    operation: &BridgeDiagnosticOperationRequest,
) -> Option<ValidityRepairFact> {
    let tolerance = brep.tolerance.linear();
    let shortest = [&operation.first_interval, &operation.second_interval]
        .into_iter()
        .filter(|interval| interval.len() >= 2)
        .flat_map(|interval| &interval[..interval.len() - 1])
        .filter_map(|half_edge| edge_length(brep, *half_edge))
        .min_by(|lhs, rhs| lhs.total_cmp(rhs));
    if shortest.is_some_and(|length| length <= tolerance) {
        return Some(fact(
            ValidityRepairReadiness::RequiresGeometryRepair,
            "bridge-interval-has-small-edge-shape-risk",
            ValidityRepairStatus::RepairRequired,
            false,
            false,
        ));
    }
    None
}

fn readiness_for_feasibility(feasibility: ConstructionFeasibility) -> ValidityRepairReadiness {
    match feasibility {
        ConstructionFeasibility::ExecutableByExistingOperator
        | ConstructionFeasibility::ExecutableAfterKnownNormalization
        | ConstructionFeasibility::RequiresSeparatorSubcomplex
        | ConstructionFeasibility::RequiresCustomAggregateAuthoring
        | ConstructionFeasibility::RequiresCycleRegrouping => ValidityRepairReadiness::ValidAsInput,
        ConstructionFeasibility::InvalidNonManifold => {
            ValidityRepairReadiness::RequiresTopologyRepair
        }
        ConstructionFeasibility::UnsupportedUnknownFamily | ConstructionFeasibility::Unknown => {
            ValidityRepairReadiness::Unknown
        }
    }
}

fn reason_for_feasibility(feasibility: ConstructionFeasibility) -> &'static str {
    match feasibility {
        ConstructionFeasibility::InvalidNonManifold => {
            "bridge-candidate-invalid-non-manifold-requires-topology-repair"
        }
        ConstructionFeasibility::UnsupportedUnknownFamily | ConstructionFeasibility::Unknown => {
            "bridge-candidate-validity-repair-not-classified"
        }
        _ => "bridge-candidate-validity-assumed-from-preflight",
    }
}

fn status_for_feasibility(feasibility: ConstructionFeasibility) -> ValidityRepairStatus {
    match feasibility {
        ConstructionFeasibility::InvalidNonManifold => ValidityRepairStatus::RepairRequired,
        ConstructionFeasibility::UnsupportedUnknownFamily | ConstructionFeasibility::Unknown => {
            ValidityRepairStatus::Unknown
        }
        _ => ValidityRepairStatus::NoRepairRequired,
    }
}

fn requires_topology_repair(reason: &'static str) -> ValidityRepairFact {
    fact(
        ValidityRepairReadiness::RequiresTopologyRepair,
        reason,
        ValidityRepairStatus::RepairRequired,
        false,
        false,
    )
}

fn fact(
    readiness: ValidityRepairReadiness,
    reason: &'static str,
    status: ValidityRepairStatus,
    changes_semantic_ownership: bool,
    changes_user_visible_intent: bool,
) -> ValidityRepairFact {
    ValidityRepairFact {
        readiness,
        reason,
        status,
        changes_semantic_ownership,
        changes_user_visible_intent,
    }
}

fn loop_walk_closes(brep: &BRep, loop_: LoopHandle) -> bool {
    let Some(loop_data) = brep.loops.get(loop_) else {
        return false;
    };
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return false;
    }
    let mut current = start;
    let mut visited = BTreeSet::new();
    loop {
        if !visited.insert(current) {
            return current == start;
        }
        let Some(half_edge) = brep.half_edges.get(current) else {
            return false;
        };
        if half_edge.loop_ != loop_ {
            return false;
        }
        current = half_edge.next;
        if current == start {
            return true;
        }
        if visited.len() > brep.half_edges.len() {
            return false;
        }
    }
}

fn directed_interval_invalid(brep: &BRep, interval: &[HalfEdgeHandle]) -> bool {
    if interval.len() < 2 {
        return true;
    }
    interval.windows(2).any(|pair| {
        let Some(first) = brep.half_edges.get(pair[0]) else {
            return true;
        };
        !brep.half_edges.contains_key(pair[1]) || first.next != pair[1]
    })
}

fn operation_half_edges(
    operation: &BridgeDiagnosticOperationRequest,
) -> impl Iterator<Item = HalfEdgeHandle> + '_ {
    [&operation.first_interval, &operation.second_interval]
        .into_iter()
        .filter(|interval| interval.len() >= 2)
        .flat_map(|interval| interval[..interval.len() - 1].iter().copied())
}

fn edge_length(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<f64> {
    let data = brep.half_edges.get(half_edge)?;
    let origin = brep.vertices.get(data.origin)?.position;
    let next = brep.half_edges.get(data.next)?;
    let destination = brep.vertices.get(next.origin)?.position;
    Some(origin.distance_to(&destination))
}
