//! Read-only near-edge classification for future BRep sewing repair.

use crate::geometry::{Point3, SurfaceEval};
use crate::topology::brep_sewing_facts::{
    sewing_pair_fact_ids, sewing_pair_knowledge, shared_fact_ids_from_pairs,
};
use crate::topology::{
    BRep, BrepEdgePairPolicyContext, BrepEdgePairPolicyDecision, BrepEdgePairPolicyDeclineReason,
    BrepEdgePairPolicyOutcome, BrepRepairFactId, BrepSewingDeclineReason, BrepSewingEdgePair,
    BrepTopologyRepairPolicy, EdgeHandle, assess_edge_pair_topology_policy, edge_segment,
};

use std::collections::HashSet;

/// Local tolerances and policy for read-only sewing classification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepSewingClassificationOptions {
    /// Tolerance under which a pair is already exact-clean.
    pub exact_tolerance: f64,
    /// Tolerance under which a pair is a possible future snap/merge candidate.
    pub snap_tolerance: f64,
    /// Number of samples for optional pcurve/surface agreement checks.
    pub samples: usize,
    /// Whether pcurve/surface sample agreement is required for classification.
    pub require_surface_agreement: bool,
    /// Whether existing twin/radial topology is allowed.
    pub allow_nonmanifold_edges: bool,
}

/// Aggregate read-only classification report.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepSewingClassificationReport {
    /// Number of caller-supplied pairs.
    pub candidate_pairs: usize,
    /// Per-pair classifications.
    pub pairs: Vec<BrepClassifiedEdgePair>,
    /// Number classified as exact-clean.
    pub exact_clean_pairs: usize,
    /// Number classified as near-clean and eligible for later snap review.
    pub near_clean_pairs: usize,
    /// Number rejected as unsafe for exact or tolerant sewing.
    pub rejected_pairs: usize,
    /// Maximum reversed endpoint gap observed.
    pub max_endpoint_gap: f64,
    /// Maximum sampled surface gap observed.
    pub max_surface_gap: f64,
    /// Stable bucket counts for GUI/API diagnostics and performance triage.
    pub buckets: Vec<BrepSewingClassificationBucketCount>,
    /// Shared repair knowledge facts referenced by this classification report.
    pub shared_fact_ids: Vec<BrepRepairFactId>,
}

/// Count for one stable sewing-classification bucket.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BrepSewingClassificationBucketCount {
    /// Stable bucket name.
    pub bucket: BrepSewingClassificationBucket,
    /// Number of candidate pairs in this bucket.
    pub count: usize,
}

/// Stable high-level classification buckets for candidate sewing pairs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepSewingClassificationBucket {
    /// Pair can be welded without vertex motion.
    ExactClean,
    /// Pair is eligible for snap review.
    NearClean,
    /// Candidate referenced the same edge twice.
    SameEdge,
    /// Required edge topology or geometry was missing.
    MissingEdge,
    /// Candidate edge had already been consumed.
    EdgeAlreadyUsed,
    /// Candidate edge was not a boundary edge.
    NotBoundaryEdge,
    /// Candidate would create unsupported non-manifold topology.
    NonManifoldConflict,
    /// Candidate touched an already joined edge.
    AlreadyJoined,
    /// Candidate topology could not be classified unambiguously.
    AmbiguousTopology,
    /// Candidate edges appear to have the same orientation.
    ForwardOrientation,
    /// Candidate endpoint gap exceeded snap tolerance.
    EndpointGapExceeded,
    /// Surface sample data needed for classification was missing.
    MissingSurfaceSampleData,
    /// Surface sample gap exceeded snap tolerance.
    SurfaceGapExceeded,
}

/// Classification for one candidate edge pair.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepClassifiedEdgePair {
    /// Original caller-supplied pair.
    pub pair: BrepSewingEdgePair,
    /// Outcome for the pair.
    pub classification: BrepEdgePairClassification,
    /// Geometric measurements used to classify the pair.
    pub measurement: BrepEdgePairMeasurement,
    /// Detailed topology policy decision when topology policy declined the pair.
    pub policy: Option<BrepEdgePairPolicyDecision>,
    /// Shared repair knowledge facts used to classify or reject the pair.
    pub shared_fact_ids: Vec<BrepRepairFactId>,
}

/// Read-only sewing classification outcome.
#[derive(Debug, Clone, PartialEq)]
pub enum BrepEdgePairClassification {
    /// Pair can be welded by exact/clean sewing.
    ExactClean,
    /// Pair is close enough for a later snap/merge policy to consider.
    NearClean,
    /// Pair is unsafe for this sewing family.
    Rejected(BrepSewingDeclineReason),
}

/// Orientation and gap measurements for one edge pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepEdgePairMeasurement {
    /// Gap if the two edges are treated as opposite directions.
    pub reversed_endpoint_gap: f64,
    /// Gap if the two edges are treated as the same direction.
    pub forward_endpoint_gap: f64,
    /// Optional reversed pcurve/surface sample gap.
    pub reversed_surface_gap: Option<f64>,
}

/// Classify caller-supplied boundary-edge pairs without mutating the BRep.
pub fn classify_boundary_edge_pairs(
    brep: &BRep,
    pairs: &[BrepSewingEdgePair],
    options: BrepSewingClassificationOptions,
) -> BrepSewingClassificationReport {
    let knowledge = sewing_pair_knowledge(brep, pairs, options.snap_tolerance, options.samples);
    let mut report = BrepSewingClassificationReport {
        candidate_pairs: pairs.len(),
        pairs: Vec::with_capacity(pairs.len()),
        exact_clean_pairs: 0,
        near_clean_pairs: 0,
        rejected_pairs: 0,
        max_endpoint_gap: 0.0,
        max_surface_gap: 0.0,
        buckets: Vec::new(),
        shared_fact_ids: Vec::new(),
    };

    for pair in pairs {
        let mut classified = classify_pair(brep, pair, options);
        classified.shared_fact_ids = sewing_pair_fact_ids(&knowledge, brep, pair);
        update_report_counts(&mut report, &classified);
        report.pairs.push(classified);
    }
    report.shared_fact_ids =
        shared_fact_ids_from_pairs(report.pairs.iter().map(|pair| pair.shared_fact_ids.clone()));
    report
}

fn classify_pair(
    brep: &BRep,
    pair: &BrepSewingEdgePair,
    options: BrepSewingClassificationOptions,
) -> BrepClassifiedEdgePair {
    let mut measurement = BrepEdgePairMeasurement {
        reversed_endpoint_gap: f64::INFINITY,
        forward_endpoint_gap: f64::INFINITY,
        reversed_surface_gap: None,
    };
    let policy_decision = topology_policy_decision(brep, pair, options);
    let policy = (!policy_decision.is_accepted()).then_some(policy_decision);
    let classification = if let Some(policy) = &policy {
        BrepEdgePairClassification::Rejected(map_policy_decline(policy))
    } else {
        match measure_pair(brep, pair, options, &mut measurement) {
            Ok(()) => classify_measured_pair(measurement, options),
            Err(reason) => BrepEdgePairClassification::Rejected(reason),
        }
    };
    BrepClassifiedEdgePair {
        pair: pair.clone(),
        classification,
        measurement,
        policy,
        shared_fact_ids: Vec::new(),
    }
}

fn measure_pair(
    brep: &BRep,
    pair: &BrepSewingEdgePair,
    options: BrepSewingClassificationOptions,
    measurement: &mut BrepEdgePairMeasurement,
) -> Result<(), BrepSewingDeclineReason> {
    let empty_used = HashSet::new();
    let policy_decision = assess_edge_pair_topology_policy(
        brep,
        pair,
        BrepTopologyRepairPolicy {
            allow_nonmanifold_edges: options.allow_nonmanifold_edges,
        },
        BrepEdgePairPolicyContext {
            used_edges: &empty_used,
        },
    );
    if !policy_decision.is_accepted() {
        return Err(map_policy_decline(&policy_decision));
    }
    let Some(first) = edge_segment(brep, pair.first) else {
        return Err(BrepSewingDeclineReason::MissingEdge);
    };
    let Some(second) = edge_segment(brep, pair.second) else {
        return Err(BrepSewingDeclineReason::MissingEdge);
    };
    measurement.reversed_endpoint_gap = reversed_endpoint_gap(first, second);
    measurement.forward_endpoint_gap = forward_endpoint_gap(first, second);
    if measurement.forward_endpoint_gap < measurement.reversed_endpoint_gap {
        return Err(BrepSewingDeclineReason::EndpointGapExceeded {
            gap: measurement.reversed_endpoint_gap,
            tolerance: options.snap_tolerance,
        });
    }
    if measurement.reversed_endpoint_gap > options.snap_tolerance {
        return Err(BrepSewingDeclineReason::EndpointGapExceeded {
            gap: measurement.reversed_endpoint_gap,
            tolerance: options.snap_tolerance,
        });
    }
    if options.require_surface_agreement {
        let Some(surface_gap) = reversed_sample_gap(brep, pair.first, pair.second, options.samples)
        else {
            return Err(BrepSewingDeclineReason::MissingSurfaceSampleData);
        };
        measurement.reversed_surface_gap = Some(surface_gap);
        if surface_gap > options.snap_tolerance {
            return Err(BrepSewingDeclineReason::SurfaceGapExceeded {
                gap: surface_gap,
                tolerance: options.snap_tolerance,
            });
        }
    }
    Ok(())
}

fn topology_policy_decision(
    brep: &BRep,
    pair: &BrepSewingEdgePair,
    options: BrepSewingClassificationOptions,
) -> BrepEdgePairPolicyDecision {
    let empty_used = HashSet::new();
    assess_edge_pair_topology_policy(
        brep,
        pair,
        BrepTopologyRepairPolicy {
            allow_nonmanifold_edges: options.allow_nonmanifold_edges,
        },
        BrepEdgePairPolicyContext {
            used_edges: &empty_used,
        },
    )
}

fn map_policy_decline(decision: &BrepEdgePairPolicyDecision) -> BrepSewingDeclineReason {
    let BrepEdgePairPolicyOutcome::Declined(reason) = &decision.outcome else {
        return BrepSewingDeclineReason::AmbiguousTopology;
    };
    match reason {
        BrepEdgePairPolicyDeclineReason::SameEdge => BrepSewingDeclineReason::SameEdge,
        BrepEdgePairPolicyDeclineReason::MissingEdge => BrepSewingDeclineReason::MissingEdge,
        BrepEdgePairPolicyDeclineReason::EdgeAlreadyUsed => {
            BrepSewingDeclineReason::EdgeAlreadyUsed
        }
        BrepEdgePairPolicyDeclineReason::NotBoundaryEdge => {
            BrepSewingDeclineReason::NotBoundaryEdge
        }
        BrepEdgePairPolicyDeclineReason::AlreadyJoined => BrepSewingDeclineReason::AlreadyJoined,
        BrepEdgePairPolicyDeclineReason::NonManifoldRadialConflict => {
            BrepSewingDeclineReason::NonManifoldConflict
        }
        BrepEdgePairPolicyDeclineReason::AmbiguousTopology => {
            BrepSewingDeclineReason::AmbiguousTopology
        }
    }
}

fn classify_measured_pair(
    measurement: BrepEdgePairMeasurement,
    options: BrepSewingClassificationOptions,
) -> BrepEdgePairClassification {
    let surface_gap = measurement.reversed_surface_gap.unwrap_or(0.0);
    if measurement.reversed_endpoint_gap <= options.exact_tolerance
        && surface_gap <= options.exact_tolerance
    {
        BrepEdgePairClassification::ExactClean
    } else {
        BrepEdgePairClassification::NearClean
    }
}

fn update_report_counts(
    report: &mut BrepSewingClassificationReport,
    classified: &BrepClassifiedEdgePair,
) {
    if classified.measurement.reversed_endpoint_gap.is_finite() {
        report.max_endpoint_gap = report
            .max_endpoint_gap
            .max(classified.measurement.reversed_endpoint_gap);
    }
    if let Some(surface_gap) = classified.measurement.reversed_surface_gap {
        report.max_surface_gap = report.max_surface_gap.max(surface_gap);
    }
    match classified.classification {
        BrepEdgePairClassification::ExactClean => report.exact_clean_pairs += 1,
        BrepEdgePairClassification::NearClean => report.near_clean_pairs += 1,
        BrepEdgePairClassification::Rejected(_) => report.rejected_pairs += 1,
    }
    increment_bucket(report, classification_bucket(classified));
}

fn increment_bucket(
    report: &mut BrepSewingClassificationReport,
    bucket: BrepSewingClassificationBucket,
) {
    if let Some(existing) = report.buckets.iter_mut().find(|item| item.bucket == bucket) {
        existing.count += 1;
    } else {
        report
            .buckets
            .push(BrepSewingClassificationBucketCount { bucket, count: 1 });
    }
}

fn classification_bucket(classified: &BrepClassifiedEdgePair) -> BrepSewingClassificationBucket {
    match &classified.classification {
        BrepEdgePairClassification::ExactClean => BrepSewingClassificationBucket::ExactClean,
        BrepEdgePairClassification::NearClean => BrepSewingClassificationBucket::NearClean,
        BrepEdgePairClassification::Rejected(reason) => rejection_bucket(reason, classified),
    }
}

fn rejection_bucket(
    reason: &BrepSewingDeclineReason,
    classified: &BrepClassifiedEdgePair,
) -> BrepSewingClassificationBucket {
    match reason {
        BrepSewingDeclineReason::SameEdge => BrepSewingClassificationBucket::SameEdge,
        BrepSewingDeclineReason::MissingEdge => BrepSewingClassificationBucket::MissingEdge,
        BrepSewingDeclineReason::EdgeAlreadyUsed => BrepSewingClassificationBucket::EdgeAlreadyUsed,
        BrepSewingDeclineReason::NotBoundaryEdge => BrepSewingClassificationBucket::NotBoundaryEdge,
        BrepSewingDeclineReason::NonManifoldConflict => {
            BrepSewingClassificationBucket::NonManifoldConflict
        }
        BrepSewingDeclineReason::AlreadyJoined => BrepSewingClassificationBucket::AlreadyJoined,
        BrepSewingDeclineReason::AmbiguousTopology => {
            BrepSewingClassificationBucket::AmbiguousTopology
        }
        BrepSewingDeclineReason::EndpointGapExceeded { .. } => {
            if classified.measurement.forward_endpoint_gap
                < classified.measurement.reversed_endpoint_gap
            {
                BrepSewingClassificationBucket::ForwardOrientation
            } else {
                BrepSewingClassificationBucket::EndpointGapExceeded
            }
        }
        BrepSewingDeclineReason::MissingSurfaceSampleData => {
            BrepSewingClassificationBucket::MissingSurfaceSampleData
        }
        BrepSewingDeclineReason::SurfaceGapExceeded { .. } => {
            BrepSewingClassificationBucket::SurfaceGapExceeded
        }
    }
}

fn reversed_sample_gap(
    brep: &BRep,
    edge_a: EdgeHandle,
    edge_b: EdgeHandle,
    samples: usize,
) -> Option<f64> {
    let a = sample_edge_on_face(brep, edge_a, samples)?;
    let b = sample_edge_on_face(brep, edge_b, samples)?;
    Some(
        a.iter()
            .zip(b.iter().rev())
            .map(|(pa, pb)| pa.distance_to(pb))
            .fold(0.0_f64, f64::max),
    )
}

fn sample_edge_on_face(brep: &BRep, edge: EdgeHandle, samples: usize) -> Option<Vec<Point3>> {
    let half_edge = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let pcurve = half_edge.pcurve.as_ref()?;
    let loop_ = brep.loops.get(half_edge.loop_)?;
    let face = brep.faces.get(loop_.face)?;
    let surface = face.surface.as_ref()?;
    let samples = samples.max(2);
    let mut points = Vec::with_capacity(samples);
    for i in 0..samples {
        let t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * i as f64 / (samples - 1) as f64;
        let uv = pcurve.point_at(t);
        points.push(surface.geometry.point_at(uv.x, uv.y));
    }
    Some(points)
}

fn reversed_endpoint_gap(a: (Point3, Point3), b: (Point3, Point3)) -> f64 {
    a.0.distance_to(&b.1).max(a.1.distance_to(&b.0))
}

fn forward_endpoint_gap(a: (Point3, Point3), b: (Point3, Point3)) -> f64 {
    a.0.distance_to(&b.0).max(a.1.distance_to(&b.1))
}

#[cfg(test)]
#[path = "brep_sewing_classification_tests.rs"]
mod tests;
