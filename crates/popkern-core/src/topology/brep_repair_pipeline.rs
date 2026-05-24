//! Generic BRep boundary repair pipeline coordinator.

use std::collections::HashSet;

use crate::topology::{
    BRep, BrepCurveTrimReconstructionOptions, BrepCurveTrimReconstructionReport,
    BrepDegenerateEdgeOptions, BrepDegenerateEdgeReport, BrepEdgeCurveCompletionOptions,
    BrepEdgeCurveCompletionReport, BrepEdgePairPolicyContext, BrepEdgePairPolicyDecision,
    BrepEdgePairPolicyOutcome, BrepEdgeSurfaceCompatibilityOptions,
    BrepEdgeSurfaceCompatibilityRepairReport, BrepLoopShellDiagnosticReport,
    BrepRepairActionReport, BrepSewingClassificationOptions, BrepSewingClassificationReport,
    BrepSewingEdgePair, BrepSnapSewingOptions, BrepSnapSewingReport, BrepTopologyDiagnosticFinding,
    BrepTopologyDiagnosticKind, BrepTopologyDiagnosticSeverity, BrepTopologyRepairPolicy,
    BrepTrimRepairOptions, BrepTrimRepairReport, analyze_brep_degenerate_edges,
    analyze_brep_loop_shell_validity, assess_edge_pair_topology_policy,
    classify_boundary_edge_pairs,
};

#[path = "brep_repair_pipeline_actions.rs"]
mod pipeline_actions;
#[path = "brep_repair_pipeline_stages.rs"]
mod pipeline_stages;
use pipeline_stages::{
    run_curve_trim_reconstruction, run_edge_curve_completion, run_edge_surface_compatibility,
    run_snap_sewing, run_trim_repair,
};

/// Options for the generic BRep boundary repair pipeline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepRepairPipelineOptions {
    /// Read-only tolerant classification options.
    pub classification: BrepSewingClassificationOptions,
    /// Transactional snap-and-sew options.
    pub snap_sewing: BrepSnapSewingOptions,
    /// Transactional trim repair options.
    pub trim_repair: BrepTrimRepairOptions,
    /// Transactional missing edge-curve completion options.
    pub edge_curve_completion: BrepEdgeCurveCompletionOptions,
    /// Transactional curve-backed pcurve reconstruction options.
    pub curve_trim_reconstruction: BrepCurveTrimReconstructionOptions,
    /// Transactional direct edge-curve to adjacent-surface repair options.
    pub edge_surface_compatibility: BrepEdgeSurfaceCompatibilityOptions,
    /// Read-only degenerate edge diagnostic options.
    pub degenerate_edges: BrepDegenerateEdgeOptions,
    /// Whether any declined candidate makes the whole pipeline fail.
    pub require_all_pairs_repaired: bool,
    /// Whether missing edge-curve completion runs after trim repair.
    pub complete_edge_curves: bool,
    /// Whether 3D-curve-backed pcurve reconstruction runs before edge/surface repair.
    pub reconstruct_curve_trims: bool,
    /// Whether direct edge-curve/surface compatibility repair runs after curve completion.
    pub repair_edge_surfaces: bool,
    /// Whether trim repair runs after sewing.
    pub repair_trims: bool,
    /// Whether loop/shell diagnostics run after all repair stages.
    pub diagnose_loop_shell: bool,
    /// Topology output policy applied to loop/shell diagnostics.
    pub topology_output: BrepTopologyOutputPolicy,
}

/// Accepted topology class for repair pipeline output.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepTopologyOutputPolicy {
    /// Output must be closed-shell suitable: no diagnostic errors or warnings.
    ClosedShell,
    /// Output may be an open sheet: boundary-chain and fragmentation warnings are allowed.
    OpenSheet,
}

/// Successful BRep repair pipeline result.
#[derive(Debug, Clone)]
pub struct BrepRepairPipelineResult {
    /// Repaired clone of the input BRep.
    pub brep: BRep,
    /// Structured stage report.
    pub report: BrepRepairPipelineReport,
}

/// Structured report for a boundary repair pipeline run.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepRepairPipelineReport {
    /// Number of caller-supplied edge-pair candidates.
    pub candidate_pairs: usize,
    /// Read-only topology policy decisions before geometry checks.
    pub policy_preflight: Vec<BrepEdgePairPolicyDecision>,
    /// Read-only tolerant classification report.
    pub classification: Option<BrepSewingClassificationReport>,
    /// Transactional snap-and-sew report.
    pub snap_sewing: Option<BrepSnapSewingReport>,
    /// Transactional missing edge-curve completion report.
    pub edge_curve_completion: Option<BrepEdgeCurveCompletionReport>,
    /// Transactional curve-backed pcurve reconstruction report.
    pub curve_trim_reconstruction: Option<BrepCurveTrimReconstructionReport>,
    /// Final transactional curve-backed pcurve reconstruction report after
    /// edge/surface compatibility repair.
    pub final_curve_trim_reconstruction: Option<BrepCurveTrimReconstructionReport>,
    /// Transactional direct edge-curve to adjacent-surface compatibility report.
    pub edge_surface_compatibility: Option<BrepEdgeSurfaceCompatibilityRepairReport>,
    /// Read-only degenerate edge diagnostics.
    pub degenerate_edges: Option<BrepDegenerateEdgeReport>,
    /// Transactional trim repair report.
    pub trim_repair: Option<BrepTrimRepairReport>,
    /// Final loop/shell diagnostics.
    pub loop_shell: Option<BrepLoopShellDiagnosticReport>,
    /// Final boundary edge count when a repaired clone exists.
    pub final_boundary_edges: Option<usize>,
    /// Generic authored/declined/rolled-back action diagnostics.
    pub actions: BrepRepairActionReport,
}

impl BrepRepairPipelineReport {
    /// Build an empty report for a fixed candidate count.
    pub fn new(candidate_pairs: usize) -> Self {
        Self {
            candidate_pairs,
            policy_preflight: Vec::new(),
            classification: None,
            snap_sewing: None,
            edge_curve_completion: None,
            curve_trim_reconstruction: None,
            final_curve_trim_reconstruction: None,
            edge_surface_compatibility: None,
            degenerate_edges: None,
            trim_repair: None,
            loop_shell: None,
            final_boundary_edges: None,
            actions: BrepRepairActionReport::default(),
        }
    }
}

/// Transactional repair pipeline failure.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepRepairPipelineFailure {
    /// Stage that declined or failed.
    pub stage: BrepRepairPipelineStage,
    /// Machine-readable reason.
    pub reason: BrepRepairPipelineFailureReason,
    /// Stage reports accumulated before failure.
    pub report: BrepRepairPipelineReport,
}

/// Pipeline stage identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepRepairPipelineStage {
    /// Preflight topology policy.
    PolicyPreflight,
    /// Tolerant geometric classification.
    Classification,
    /// Transactional snap-and-sew.
    SnapSewing,
    /// Transactional missing edge-curve completion.
    EdgeCurveCompletion,
    /// Transactional curve-backed pcurve reconstruction.
    CurveTrimReconstruction,
    /// Transactional direct edge-curve to adjacent-surface compatibility repair.
    EdgeSurfaceCompatibility,
    /// Transactional trim repair.
    TrimRepair,
    /// Read-only loop/shell diagnostics.
    LoopShellDiagnostics,
}

/// Reason the repair pipeline failed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BrepRepairPipelineFailureReason {
    /// Topology policy declined one or more candidates.
    PolicyDeclined,
    /// Classification rejected one or more candidates.
    ClassificationRejected,
    /// Snap-and-sew declined one or more candidates.
    SnapSewingDeclined,
    /// Snap-and-sew failed validation.
    SnapSewingFailed(String),
    /// Missing edge-curve completion failed validation or left unresolved edges.
    EdgeCurveCompletionFailed(String),
    /// Curve-backed trim reconstruction failed validation or left unresolved trims.
    CurveTrimReconstructionFailed(String),
    /// Edge/surface compatibility repair failed validation or left unresolved edges.
    EdgeSurfaceCompatibilityFailed(String),
    /// Trim repair failed validation or left unresolved trim findings.
    TrimRepairFailed(String),
    /// Loop/shell diagnostics found error-severity topology issues.
    LoopShellErrors,
}

/// Run the generic BRep boundary repair pipeline transactionally.
pub fn repair_brep_boundary_edges_transactionally(
    brep: &BRep,
    pairs: &[BrepSewingEdgePair],
    options: BrepRepairPipelineOptions,
) -> Result<BrepRepairPipelineResult, Box<BrepRepairPipelineFailure>> {
    let mut report = BrepRepairPipelineReport::new(pairs.len());

    report.policy_preflight = preflight_policy(brep, pairs, options);
    if options.require_all_pairs_repaired
        && report
            .policy_preflight
            .iter()
            .any(|decision| !decision.is_accepted())
    {
        return Err(failure(
            BrepRepairPipelineStage::PolicyPreflight,
            BrepRepairPipelineFailureReason::PolicyDeclined,
            report,
        ));
    }

    let classification = classify_boundary_edge_pairs(brep, pairs, options.classification);
    let classified_pairs = classification.pairs.clone();
    if options.require_all_pairs_repaired && classification.rejected_pairs > 0 {
        report.classification = Some(classification);
        return Err(failure(
            BrepRepairPipelineStage::Classification,
            BrepRepairPipelineFailureReason::ClassificationRejected,
            report,
        ));
    }
    report.classification = Some(classification);

    let sewn = run_snap_sewing(brep, &classified_pairs, options, &mut report)?;
    let repaired = if options.complete_edge_curves {
        run_edge_curve_completion(&sewn, options, &mut report)?
    } else {
        sewn
    };
    let repaired = if options.reconstruct_curve_trims {
        run_curve_trim_reconstruction(&repaired, options, &mut report)?
    } else {
        repaired
    };
    let repaired = if options.repair_edge_surfaces {
        run_edge_surface_compatibility(&repaired, options, &mut report)?
    } else {
        repaired
    };
    let repaired = if options.reconstruct_curve_trims && options.repair_edge_surfaces {
        let first_reconstruction = report.curve_trim_reconstruction.take();
        let final_repaired = run_curve_trim_reconstruction(&repaired, options, &mut report)?;
        report.final_curve_trim_reconstruction = report.curve_trim_reconstruction.take();
        report.curve_trim_reconstruction = first_reconstruction;
        final_repaired
    } else {
        repaired
    };
    let repaired = if options.repair_trims {
        run_trim_repair(&repaired, options, &mut report)?
    } else {
        repaired
    };
    report.degenerate_edges = Some(analyze_brep_degenerate_edges(
        &repaired,
        options.degenerate_edges,
    ));

    if options.diagnose_loop_shell {
        let diagnostics = analyze_brep_loop_shell_validity(&repaired);
        if topology_report_blocks_output(&diagnostics, options.topology_output) {
            report.loop_shell = Some(diagnostics);
            return Err(failure(
                BrepRepairPipelineStage::LoopShellDiagnostics,
                BrepRepairPipelineFailureReason::LoopShellErrors,
                report,
            ));
        }
        report.loop_shell = Some(diagnostics);
    }

    report.final_boundary_edges = Some(repaired.boundary_edges.len());
    Ok(BrepRepairPipelineResult {
        brep: repaired,
        report,
    })
}

fn topology_report_blocks_output(
    report: &BrepLoopShellDiagnosticReport,
    policy: BrepTopologyOutputPolicy,
) -> bool {
    if policy == BrepTopologyOutputPolicy::ClosedShell && report.boundary_edges > 0 {
        return true;
    }
    report
        .findings
        .iter()
        .any(|finding| topology_finding_blocks_output(finding, policy))
}

fn topology_finding_blocks_output(
    finding: &BrepTopologyDiagnosticFinding,
    policy: BrepTopologyOutputPolicy,
) -> bool {
    if finding.severity == BrepTopologyDiagnosticSeverity::Error {
        return true;
    }
    match policy {
        BrepTopologyOutputPolicy::ClosedShell => {
            finding.severity != BrepTopologyDiagnosticSeverity::Info
        }
        BrepTopologyOutputPolicy::OpenSheet => !matches!(
            finding.kind,
            BrepTopologyDiagnosticKind::DanglingBoundaryChain
                | BrepTopologyDiagnosticKind::FragmentedShell
        ),
    }
}

fn preflight_policy(
    brep: &BRep,
    pairs: &[BrepSewingEdgePair],
    options: BrepRepairPipelineOptions,
) -> Vec<BrepEdgePairPolicyDecision> {
    let mut used_edges = HashSet::new();
    let mut decisions = Vec::with_capacity(pairs.len());
    for pair in pairs {
        let decision = assess_edge_pair_topology_policy(
            brep,
            pair,
            BrepTopologyRepairPolicy {
                allow_nonmanifold_edges: options.classification.allow_nonmanifold_edges,
            },
            BrepEdgePairPolicyContext {
                used_edges: &used_edges,
            },
        );
        if matches!(decision.outcome, BrepEdgePairPolicyOutcome::Accepted) {
            used_edges.insert(pair.first);
            used_edges.insert(pair.second);
        }
        decisions.push(decision);
    }
    decisions
}

fn failure(
    stage: BrepRepairPipelineStage,
    reason: BrepRepairPipelineFailureReason,
    report: BrepRepairPipelineReport,
) -> Box<BrepRepairPipelineFailure> {
    Box::new(BrepRepairPipelineFailure {
        stage,
        reason,
        report,
    })
}

#[cfg(test)]
#[path = "brep_repair_pipeline_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "brep_repair_pipeline_curve_trim_tests.rs"]
mod curve_trim_tests;

#[cfg(test)]
#[path = "brep_repair_pipeline_edge_curve_tests.rs"]
mod edge_curve_tests;
