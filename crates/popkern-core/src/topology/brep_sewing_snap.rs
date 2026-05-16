//! Transactional tolerant snap-and-sew for classified BRep edge pairs.

use crate::geometry::Point3;
use crate::topology::brep_sewing_facts::{sewing_pair_fact_ids, sewing_pair_knowledge};
use crate::topology::{
    BRep, BrepClassifiedEdgePair, BrepEdgePairClassification, BrepRepairFactId,
    BrepSewingDeclineReason, BrepSewingEdgePair, BrepSewingOptions, BrepSewingReport, EdgeHandle,
    VertexHandle, sew_clean_boundary_edge_pairs,
};

/// Policy for transactional tolerant snap sewing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepSnapSewingOptions {
    /// Maximum allowed movement for any snapped vertex.
    pub max_vertex_motion: f64,
    /// Options passed to the final clean sewing operation.
    pub clean_sewing: BrepSewingOptions,
}

/// Successful transactional snap sewing result.
#[derive(Debug, Clone)]
pub struct BrepSnapSewingResult {
    /// Sewn clone of the input BRep.
    pub brep: BRep,
    /// Structured snap and final sewing diagnostics.
    pub report: BrepSnapSewingReport,
}

/// Aggregate diagnostics for tolerant snap sewing.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepSnapSewingReport {
    /// Number of classified pairs considered.
    pub classified_pairs: usize,
    /// Number of exact-clean pairs forwarded without vertex motion.
    pub exact_pairs: usize,
    /// Number of near-clean pairs whose vertices were snapped.
    pub snapped_pairs: usize,
    /// Pair-level declines before clean sewing.
    pub declined_pairs: Vec<BrepSnapSewingDecline>,
    /// Maximum vertex movement applied to the cloned BRep.
    pub max_vertex_motion: f64,
    /// Final clean sewing report when that stage completes.
    pub clean_report: Option<BrepSewingReport>,
    /// Shared repair knowledge facts referenced by this snap-sewing attempt.
    pub shared_fact_ids: Vec<BrepRepairFactId>,
}

/// Pair-level snap-sewing decline.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepSnapSewingDecline {
    /// Candidate pair that declined.
    pub pair: BrepSewingEdgePair,
    /// Reason the pair was not snapped or forwarded.
    pub reason: BrepSnapSewingDeclineReason,
    /// Shared repair knowledge facts used to decline this pair.
    pub shared_fact_ids: Vec<BrepRepairFactId>,
}

/// Reason a classified pair was not accepted for snap sewing.
#[derive(Debug, Clone, PartialEq)]
pub enum BrepSnapSewingDeclineReason {
    /// Classification already rejected the pair.
    RejectedClassification(BrepSewingDeclineReason),
    /// Edge handles were missing or incomplete.
    MissingEdge,
    /// Pair would move a vertex farther than the configured limit.
    VertexMotionExceeded {
        /// Largest proposed vertex motion.
        motion: f64,
        /// Configured maximum allowed vertex motion.
        limit: f64,
    },
}

/// Transaction-level failure from final clean sewing validation.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepSnapSewingFailure {
    /// Failure reason.
    pub reason: BrepSnapSewingFailureReason,
    /// Report accumulated before the transaction failed.
    pub report: BrepSnapSewingReport,
}

/// Reason the snap-sewing transaction failed.
#[derive(Debug, Clone, PartialEq)]
pub enum BrepSnapSewingFailureReason {
    /// The final clean sewing stage failed.
    CleanSewingFailed(String),
}

/// Snap classified near-clean pairs on a clone, then run clean edge sewing.
pub fn snap_and_sew_classified_edge_pairs(
    brep: &BRep,
    classified_pairs: &[BrepClassifiedEdgePair],
    options: BrepSnapSewingOptions,
) -> Result<BrepSnapSewingResult, Box<BrepSnapSewingFailure>> {
    let pairs = classified_pairs
        .iter()
        .map(|classified| classified.pair.clone())
        .collect::<Vec<_>>();
    let knowledge = sewing_pair_knowledge(
        brep,
        &pairs,
        options
            .clean_sewing
            .tolerance
            .max(options.max_vertex_motion),
        options.clean_sewing.samples,
    );
    if !knowledge.validate_fresh(brep) {
        return Err(Box::new(BrepSnapSewingFailure {
            reason: BrepSnapSewingFailureReason::CleanSewingFailed(
                "stale snap-sewing repair knowledge graph".to_string(),
            ),
            report: empty_report(classified_pairs),
        }));
    }
    let mut candidate = brep.clone();
    let mut clean_pairs = Vec::new();
    let mut report = BrepSnapSewingReport {
        classified_pairs: classified_pairs.len(),
        exact_pairs: 0,
        snapped_pairs: 0,
        declined_pairs: Vec::new(),
        max_vertex_motion: 0.0,
        clean_report: None,
        shared_fact_ids: Vec::new(),
    };

    for classified in classified_pairs {
        let shared_fact_ids = sewing_pair_fact_ids(&knowledge, brep, &classified.pair);
        match &classified.classification {
            BrepEdgePairClassification::ExactClean => {
                clean_pairs.push(classified.pair.clone());
                report.exact_pairs += 1;
            }
            BrepEdgePairClassification::NearClean => {
                match snap_pair_vertices(&mut candidate, &classified.pair, options) {
                    Ok(max_motion) => {
                        report.max_vertex_motion = report.max_vertex_motion.max(max_motion);
                        report.snapped_pairs += 1;
                        clean_pairs.push(classified.pair.clone());
                    }
                    Err(reason) => report.declined_pairs.push(BrepSnapSewingDecline {
                        pair: classified.pair.clone(),
                        reason,
                        shared_fact_ids: shared_fact_ids.clone(),
                    }),
                }
            }
            BrepEdgePairClassification::Rejected(reason) => {
                report.declined_pairs.push(BrepSnapSewingDecline {
                    pair: classified.pair.clone(),
                    reason: BrepSnapSewingDeclineReason::RejectedClassification(reason.clone()),
                    shared_fact_ids: shared_fact_ids.clone(),
                });
            }
        }
        report.shared_fact_ids.extend(shared_fact_ids);
    }

    let sewn = sew_clean_boundary_edge_pairs(&candidate, &clean_pairs, options.clean_sewing)
        .map_err(|failure| {
            Box::new(BrepSnapSewingFailure {
                reason: BrepSnapSewingFailureReason::CleanSewingFailed(format!(
                    "{:?}",
                    failure.reason
                )),
                report: BrepSnapSewingReport {
                    clean_report: Some(failure.report),
                    ..report.clone()
                },
            })
        })?;
    report
        .shared_fact_ids
        .extend(sewn.report.shared_fact_ids.iter().copied());
    report.shared_fact_ids.sort();
    report.shared_fact_ids.dedup();
    report.clean_report = Some(sewn.report);
    Ok(BrepSnapSewingResult {
        brep: sewn.brep,
        report,
    })
}

fn empty_report(classified_pairs: &[BrepClassifiedEdgePair]) -> BrepSnapSewingReport {
    BrepSnapSewingReport {
        classified_pairs: classified_pairs.len(),
        exact_pairs: 0,
        snapped_pairs: 0,
        declined_pairs: Vec::new(),
        max_vertex_motion: 0.0,
        clean_report: None,
        shared_fact_ids: Vec::new(),
    }
}

fn snap_pair_vertices(
    brep: &mut BRep,
    pair: &BrepSewingEdgePair,
    options: BrepSnapSewingOptions,
) -> Result<f64, BrepSnapSewingDeclineReason> {
    let Some(first) = edge_endpoint_handles(brep, pair.first) else {
        return Err(BrepSnapSewingDeclineReason::MissingEdge);
    };
    let Some(second) = edge_endpoint_handles(brep, pair.second) else {
        return Err(BrepSnapSewingDeclineReason::MissingEdge);
    };
    let first_start = brep.vertices[first.start].position;
    let first_end = brep.vertices[first.end].position;
    let second_start = brep.vertices[second.start].position;
    let second_end = brep.vertices[second.end].position;
    let start_target = midpoint(first_start, second_end);
    let end_target = midpoint(first_end, second_start);
    let max_motion = first_start
        .distance_to(&start_target)
        .max(second_end.distance_to(&start_target))
        .max(first_end.distance_to(&end_target))
        .max(second_start.distance_to(&end_target));
    if max_motion > options.max_vertex_motion {
        return Err(BrepSnapSewingDeclineReason::VertexMotionExceeded {
            motion: max_motion,
            limit: options.max_vertex_motion,
        });
    }
    move_vertex_with_tolerance(brep, first.start, start_target, max_motion);
    move_vertex_with_tolerance(brep, second.end, start_target, max_motion);
    move_vertex_with_tolerance(brep, first.end, end_target, max_motion);
    move_vertex_with_tolerance(brep, second.start, end_target, max_motion);
    Ok(max_motion)
}

#[derive(Debug, Clone, Copy)]
struct EdgeEndpoints {
    start: VertexHandle,
    end: VertexHandle,
}

fn edge_endpoint_handles(brep: &BRep, edge: EdgeHandle) -> Option<EdgeEndpoints> {
    let half_edge = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let next = brep.half_edges.get(half_edge.next)?;
    Some(EdgeEndpoints {
        start: half_edge.origin,
        end: next.origin,
    })
}

fn midpoint(a: Point3, b: Point3) -> Point3 {
    Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z))
}

fn move_vertex_with_tolerance(brep: &mut BRep, vertex: VertexHandle, target: Point3, motion: f64) {
    let default = brep.tolerance.linear();
    let current = brep.vertices[vertex].tolerance.unwrap_or(default);
    brep.vertices[vertex].position = target;
    brep.vertices[vertex].tolerance = Some(current.max(motion));
    if motion > brep.max_imprecision {
        brep.max_imprecision = (motion * 1.01).max(motion + brep.tolerance.numerical_floor());
    }
}

#[cfg(test)]
#[path = "brep_sewing_snap_tests.rs"]
mod tests;
