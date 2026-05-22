//! Transactional coordinated vertex-star repair for edge/surface compatibility.

use std::collections::HashSet;

use crate::topology::brep_edge_surface_compatibility::{
    BrepEdgeSurfaceCompatibilityOptions, BrepEdgeSurfaceCompatibilityRepairDecision,
    BrepEdgeSurfaceCompatibilityRepairOutcome, BrepEdgeSurfaceCompatibilityRepairReport,
    analyze_brep_edge_surface_compatibility,
};
use crate::topology::brep_edge_surface_compatibility_helpers::{
    edge_vertices, move_vertex_with_tolerance,
};
use crate::topology::brep_edge_surface_compatibility_repair::repair_edge;
use crate::topology::brep_edge_surface_compatibility_vertex_star::BrepEdgeSurfaceVertexStarOutcome;
use crate::topology::{BRep, BrepRepairKnowledgeGraph, EdgeHandle};

pub(super) fn repair_coherent_vertex_stars(
    brep: &mut BRep,
    knowledge: &BrepRepairKnowledgeGraph,
    report: &mut BrepEdgeSurfaceCompatibilityRepairReport,
    options: BrepEdgeSurfaceCompatibilityOptions,
) {
    let coherent_stars: Vec<_> = report
        .vertex_stars
        .decisions
        .iter()
        .filter(|star| {
            matches!(
                star.outcome,
                BrepEdgeSurfaceVertexStarOutcome::CoherentCandidate
            )
        })
        .cloned()
        .collect();
    if coherent_stars.is_empty() {
        return;
    }
    let before_stars = brep.clone();
    let coherent_vertices: HashSet<_> = coherent_stars.iter().map(|star| star.vertex).collect();
    let mut affected_edges = Vec::new();
    for star in &coherent_stars {
        let Some(proposed) = star.proposed_position else {
            *brep = before_stars;
            return;
        };
        if !brep.vertices.contains_key(star.vertex) {
            *brep = before_stars;
            return;
        }
        let current = brep.vertices[star.vertex].position;
        move_vertex_with_tolerance(brep, star.vertex, proposed, current.distance_to(&proposed));
        affected_edges.extend(star.incompatible_edges.iter().copied());
    }
    affected_edges.sort_by_key(|edge| format!("{edge:?}"));
    affected_edges.dedup();
    affected_edges.retain(|edge| {
        edge_vertices(brep, *edge).is_some_and(|(start, end)| {
            coherent_vertices.contains(&start) && coherent_vertices.contains(&end)
        })
    });
    if affected_edges.is_empty() {
        *brep = before_stars;
        return;
    }

    let mut local_decisions = Vec::new();
    for edge in &affected_edges {
        if edge_is_currently_incompatible(brep, *edge, options) {
            local_decisions.push(repair_edge(brep, knowledge, *edge, options));
        }
    }
    let local_report = analyze_brep_edge_surface_compatibility(brep, options);
    let local_clean = affected_edges.iter().all(|edge| {
        !local_report
            .findings
            .iter()
            .any(|finding| finding.edge == *edge)
    });
    let validates = !options.validate_output || brep.validate().is_ok();
    if local_clean && validates {
        for decision in local_decisions {
            record_repair_decision(report, decision);
        }
    } else {
        *brep = before_stars;
    }
}

pub(super) fn edge_is_currently_incompatible(
    brep: &BRep,
    edge: EdgeHandle,
    options: BrepEdgeSurfaceCompatibilityOptions,
) -> bool {
    analyze_brep_edge_surface_compatibility(brep, options)
        .findings
        .iter()
        .any(|finding| finding.edge == edge)
}

pub(super) fn record_repair_decision(
    report: &mut BrepEdgeSurfaceCompatibilityRepairReport,
    decision: BrepEdgeSurfaceCompatibilityRepairDecision,
) {
    report.attempted_edges += 1;
    match decision.outcome {
        BrepEdgeSurfaceCompatibilityRepairOutcome::RepairedFromSurfaceIntersection
        | BrepEdgeSurfaceCompatibilityRepairOutcome::RepairedFromSurfaceBoundary => {
            report.repaired_edges += 1;
            report.max_endpoint_gap = report.max_endpoint_gap.max(decision.endpoint_gap);
        }
        BrepEdgeSurfaceCompatibilityRepairOutcome::Declined(_) => {
            report.declined_edges += 1;
        }
    }
    report.decisions.push(decision);
}
