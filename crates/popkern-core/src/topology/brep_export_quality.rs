//! Read-only BRep export-quality diagnostics.

use crate::topology::{
    BRep, BrepFaceUseTrimAuthorityFindingKind, BrepFaceUseTrimAuthorityOptions,
    BrepFaceUseTrimAuthorityReport, BrepLoopShellDiagnosticReport, BrepTopologyDiagnosticKind,
    BrepTopologyDiagnosticSeverity, BrepTopologyOutputPolicy, BrepTrimDiagnosticReport,
    BrepTrimRepairOptions, analyze_brep_face_use_trim_authority, analyze_brep_loop_shell_validity,
    analyze_brep_trims, edge_segment,
};
use serde::{Deserialize, Serialize};

#[path = "brep_export_quality_authority.rs"]
mod brep_export_quality_authority;
#[path = "brep_export_quality_helpers.rs"]
mod brep_export_quality_helpers;
use brep_export_quality_authority::{authoritative_degenerate_edges, export_authority_summary};
use brep_export_quality_helpers::{edge_surface_fidelity_blockers, pcurve_fidelity_blockers};

/// Options for export-quality diagnostics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepExportQualityOptions {
    /// Topology class expected by the export consumer.
    pub topology_output: BrepTopologyOutputPolicy,
    /// Trim analysis tolerance and sampling options.
    pub trim: BrepTrimRepairOptions,
    /// Whether every edge must carry a 3D curve.
    pub require_edge_curves: bool,
    /// Whether missing 3D curves on valid collapsed surface-trim edges are exportable.
    pub allow_collapsed_trim_edge_curves: bool,
    /// Optional maximum accepted edge-to-pcurve geometric deviation.
    pub max_pcurve_deviation: Option<f64>,
    /// Optional maximum accepted direct edge-curve to adjacent-surface deviation.
    pub max_edge_surface_deviation: Option<f64>,
    /// Whether to call full BRep validation.
    pub validate_brep: bool,
}

/// Export readiness classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepExportReadiness {
    /// Ready for the requested output policy.
    Ready,
    /// Export is possible as an open sheet with expected boundary diagnostics.
    ReadyOpenSheet,
    /// Export should be blocked until findings are resolved.
    Blocked,
}

/// Aggregate read-only export-quality report.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepExportQualityReport {
    /// Final readiness classification.
    pub readiness: BrepExportReadiness,
    /// Requested topology output policy.
    pub topology_output: BrepTopologyOutputPolicy,
    /// Number of edges in the BRep.
    pub edges: usize,
    /// Number of boundary edges in the BRep.
    pub boundary_edges: usize,
    /// Number of edges missing 3D curve geometry.
    pub missing_edge_curves: usize,
    /// Missing edge curves accepted as valid collapsed surface trims.
    pub collapsed_trim_missing_edge_curves: usize,
    /// Loop/shell diagnostics.
    pub loop_shell: BrepLoopShellDiagnosticReport,
    /// Trim diagnostics.
    pub trims: BrepTrimDiagnosticReport,
    /// Face-use trim authority diagnostics.
    pub face_use_trim_authority: BrepFaceUseTrimAuthorityReport,
    /// Whether full BRep validation passed.
    pub validation_passed: bool,
    /// Validation failure text, when validation was requested and failed.
    pub validation_error: Option<String>,
    /// Number of blocking topology findings under the requested output policy.
    pub topology_blockers: usize,
    /// Number of blocking trim findings.
    pub trim_blockers: usize,
    /// Number of blocking face-use trim authority findings.
    pub face_use_trim_authority_blockers: usize,
    /// Number of blocking edge-curve findings.
    pub edge_curve_blockers: usize,
    /// Number of pcurve fidelity blockers.
    pub pcurve_blockers: usize,
    /// Maximum measured edge-to-pcurve deviation.
    pub max_pcurve_deviation: f64,
    /// Number of direct edge-curve to adjacent-surface fidelity blockers.
    pub edge_surface_blockers: usize,
    /// Maximum measured direct edge-curve to adjacent-surface deviation.
    pub max_edge_surface_deviation: f64,
    /// Authority-fact summary used by export-readiness decisions.
    pub authority: BrepExportAuthoritySummary,
}

/// Serializable debug artifact for export-quality inspection.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BrepExportQualityArtifact {
    /// Compact summary for logs, sidecars, and UI panels.
    pub summary: BrepExportQualitySummary,
    /// Boundary-edge overlay geometry.
    pub boundary_edges: Vec<BrepBoundaryEdgeArtifact>,
    /// Machine-readable blocker summary.
    pub blockers: BrepExportBlockerSummary,
    /// Authority-fact summary used by export-readiness decisions.
    pub authority: BrepExportAuthoritySummary,
}

/// Serializable authority summary for export-quality inspection.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct BrepExportAuthoritySummary {
    /// Total authority facts in the all-edge export scope.
    pub facts: usize,
    /// Existing pcurves that are stale against their edge/surface dependencies.
    pub stale_pcurves: usize,
    /// Contradictory authority facts.
    pub conflicting_sources: usize,
    /// Missing non-boundary authority sources.
    pub missing_sources: usize,
    /// Explicit degenerate-topology authority facts.
    pub degenerate_topology: usize,
}

/// Serializable export-quality summary.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BrepExportQualitySummary {
    /// Readiness as stable text.
    pub readiness: String,
    /// Topology output policy as stable text.
    pub topology_output: String,
    /// Edge count.
    pub edges: usize,
    /// Boundary edge count.
    pub boundary_edges: usize,
    /// Missing edge-curve count.
    pub missing_edge_curves: usize,
    /// Missing edge curves accepted as valid collapsed surface trims.
    pub collapsed_trim_missing_edge_curves: usize,
    /// Face-use trim authority blocker count.
    pub face_use_trim_authority_blockers: usize,
    /// Maximum measured edge-to-pcurve deviation.
    pub max_pcurve_deviation: f64,
    /// Maximum measured direct edge-curve to adjacent-surface deviation.
    pub max_edge_surface_deviation: f64,
    /// Validation result.
    pub validation_passed: bool,
    /// Optional validation error.
    pub validation_error: Option<String>,
}

/// Serializable boundary edge overlay record.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BrepBoundaryEdgeArtifact {
    /// Stable debug handle string.
    pub edge: String,
    /// Stable debug half-edge handle for the surviving boundary side.
    pub half_edge: String,
    /// Stable debug loop handle that owns the surviving boundary side.
    pub loop_: Option<String>,
    /// Stable debug face handle that owns the surviving boundary side.
    pub face: Option<String>,
    /// Semantic role of the boundary edge.
    pub edge_role: String,
    /// Semantic role of the owner face, when known.
    pub face_role: Option<String>,
    /// Whether the owner face carries surface geometry.
    pub face_has_surface: bool,
    /// Start point.
    pub start: [f64; 3],
    /// End point.
    pub end: [f64; 3],
    /// Whether the edge carries 3D curve geometry.
    pub has_curve: bool,
}

/// Serializable blocker summary.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct BrepExportBlockerSummary {
    /// Topology blockers.
    pub topology: usize,
    /// Trim blockers.
    pub trims: usize,
    /// Missing edge-curve blockers.
    pub edge_curves: usize,
    /// Pcurve fidelity blockers.
    pub pcurves: usize,
    /// Face-use trim authority blockers.
    pub face_use_trim_authority: usize,
    /// Direct edge-curve to adjacent-surface fidelity blockers.
    pub edge_surfaces: usize,
    /// Validation blockers.
    pub validation: usize,
}

impl BrepExportQualityReport {
    /// True when no export blockers remain.
    pub fn is_ready(&self) -> bool {
        self.readiness != BrepExportReadiness::Blocked
    }
}

/// Analyze whether a BRep is ready for export under the requested policy.
pub fn analyze_brep_export_quality(
    brep: &BRep,
    options: BrepExportQualityOptions,
) -> BrepExportQualityReport {
    let pcurve_fidelity = pcurve_fidelity_blockers(brep, options.max_pcurve_deviation);
    analyze_brep_export_quality_with_pcurve_fidelity(brep, options, pcurve_fidelity)
}

/// Analyze export quality using a caller-supplied pcurve-fidelity summary.
pub fn analyze_brep_export_quality_with_pcurve_fidelity(
    brep: &BRep,
    options: BrepExportQualityOptions,
    pcurve_fidelity: (usize, f64),
) -> BrepExportQualityReport {
    let trims = analyze_brep_trims(brep, options.trim);
    let edge_surface_fidelity =
        edge_surface_fidelity_blockers(brep, options.max_edge_surface_deviation, options.trim);
    analyze_brep_export_quality_with_local_diagnostics(
        brep,
        options,
        trims,
        pcurve_fidelity,
        edge_surface_fidelity,
    )
}

/// Analyze export quality using caller-supplied local geometry diagnostics.
pub fn analyze_brep_export_quality_with_local_diagnostics(
    brep: &BRep,
    options: BrepExportQualityOptions,
    trims: BrepTrimDiagnosticReport,
    pcurve_fidelity: (usize, f64),
    edge_surface_fidelity: (usize, f64),
) -> BrepExportQualityReport {
    let loop_shell = analyze_brep_loop_shell_validity(brep);
    let face_use_trim_authority = analyze_brep_face_use_trim_authority(
        brep,
        BrepFaceUseTrimAuthorityOptions {
            tolerance: options.trim.endpoint_tolerance,
            samples: options.trim.samples,
        },
    );
    let missing_edge_curves = brep
        .edges
        .iter()
        .filter(|(_, edge)| edge.curve.is_none())
        .count();
    let collapsed_trim_missing_edge_curves = if options.allow_collapsed_trim_edge_curves {
        collapsed_trim_missing_edge_curve_count(brep, options.trim)
    } else {
        0
    };
    let validation_error = if options.validate_brep {
        brep.validate().err().map(|err| err.to_string())
    } else {
        None
    };
    let topology_blockers = topology_blocker_count(&loop_shell, options.topology_output);
    let trim_blockers = trims.findings.len();
    let face_use_trim_authority_blockers =
        face_use_trim_authority_blocker_count(&face_use_trim_authority);
    let edge_curve_blockers = if options.require_edge_curves {
        missing_edge_curves.saturating_sub(collapsed_trim_missing_edge_curves)
    } else {
        0
    };
    let (pcurve_blockers, max_pcurve_deviation) = pcurve_fidelity;
    let (edge_surface_blockers, max_edge_surface_deviation) = edge_surface_fidelity;
    let authority = export_authority_summary(brep, options.trim);
    let validation_passed = validation_error.is_none();
    let blocked = topology_blockers > 0
        || trim_blockers > 0
        || face_use_trim_authority_blockers > 0
        || edge_curve_blockers > 0
        || pcurve_blockers > 0
        || edge_surface_blockers > 0
        || !validation_passed;
    let readiness = if blocked {
        BrepExportReadiness::Blocked
    } else if options.topology_output == BrepTopologyOutputPolicy::OpenSheet
        && loop_shell.boundary_edges > 0
    {
        BrepExportReadiness::ReadyOpenSheet
    } else {
        BrepExportReadiness::Ready
    };
    BrepExportQualityReport {
        readiness,
        topology_output: options.topology_output,
        edges: brep.edges.len(),
        boundary_edges: loop_shell.boundary_edges,
        missing_edge_curves,
        collapsed_trim_missing_edge_curves,
        loop_shell,
        trims,
        face_use_trim_authority,
        validation_passed,
        validation_error,
        topology_blockers,
        trim_blockers,
        face_use_trim_authority_blockers,
        edge_curve_blockers,
        pcurve_blockers,
        max_pcurve_deviation,
        edge_surface_blockers,
        max_edge_surface_deviation,
        authority,
    }
}

/// Build a deterministic, serializable debug artifact from a report.
pub fn build_brep_export_quality_artifact(
    brep: &BRep,
    report: &BrepExportQualityReport,
) -> BrepExportQualityArtifact {
    let mut boundary_edges: Vec<_> = brep
        .boundary_edges
        .iter()
        .filter_map(|&edge| {
            let (start, end) = edge_segment(brep, edge)?;
            let edge_data = brep.edges.get(edge)?;
            let half_edge = brep.half_edges.get(edge_data.half_edge)?;
            let loop_ = (half_edge.loop_ != crate::topology::LoopHandle::default())
                .then(|| format!("{:?}", half_edge.loop_));
            let face = brep
                .loops
                .get(half_edge.loop_)
                .map(|loop_data| loop_data.face);
            let face_data = face.and_then(|face| brep.faces.get(face));
            Some(BrepBoundaryEdgeArtifact {
                edge: format!("{edge:?}"),
                half_edge: format!("{:?}", edge_data.half_edge),
                loop_,
                face: face.map(|face| format!("{face:?}")),
                edge_role: format!("{:?}", edge_data.role),
                face_role: face_data.map(|face| format!("{:?}", face.role)),
                face_has_surface: face_data.is_some_and(|face| face.surface.is_some()),
                start: [start.x, start.y, start.z],
                end: [end.x, end.y, end.z],
                has_curve: edge_data.curve.is_some(),
            })
        })
        .collect();
    boundary_edges.sort_by(|a, b| a.edge.cmp(&b.edge));
    BrepExportQualityArtifact {
        summary: BrepExportQualitySummary {
            readiness: readiness_label(report.readiness).to_string(),
            topology_output: topology_output_label(report.topology_output).to_string(),
            edges: report.edges,
            boundary_edges: report.boundary_edges,
            missing_edge_curves: report.missing_edge_curves,
            collapsed_trim_missing_edge_curves: report.collapsed_trim_missing_edge_curves,
            face_use_trim_authority_blockers: report.face_use_trim_authority_blockers,
            max_pcurve_deviation: report.max_pcurve_deviation,
            max_edge_surface_deviation: report.max_edge_surface_deviation,
            validation_passed: report.validation_passed,
            validation_error: report.validation_error.clone(),
        },
        boundary_edges,
        blockers: BrepExportBlockerSummary {
            topology: report.topology_blockers,
            trims: report.trim_blockers,
            face_use_trim_authority: report.face_use_trim_authority_blockers,
            edge_curves: report.edge_curve_blockers,
            pcurves: report.pcurve_blockers,
            edge_surfaces: report.edge_surface_blockers,
            validation: usize::from(!report.validation_passed),
        },
        authority: report.authority.clone(),
    }
}

fn topology_blocker_count(
    report: &BrepLoopShellDiagnosticReport,
    policy: BrepTopologyOutputPolicy,
) -> usize {
    let boundary_blocker =
        usize::from(policy == BrepTopologyOutputPolicy::ClosedShell && report.boundary_edges > 0);
    boundary_blocker
        + report
            .findings
            .iter()
            .filter(|finding| {
                topology_finding_blocks_export(finding.severity, &finding.kind, policy)
            })
            .count()
}

fn face_use_trim_authority_blocker_count(report: &BrepFaceUseTrimAuthorityReport) -> usize {
    report
        .findings
        .iter()
        .filter(|finding| {
            matches!(
                finding.kind,
                BrepFaceUseTrimAuthorityFindingKind::LoopOwnershipMismatch
                    | BrepFaceUseTrimAuthorityFindingKind::PcurveSurfaceMismatch
            )
        })
        .count()
}

fn topology_finding_blocks_export(
    severity: BrepTopologyDiagnosticSeverity,
    kind: &BrepTopologyDiagnosticKind,
    policy: BrepTopologyOutputPolicy,
) -> bool {
    if severity == BrepTopologyDiagnosticSeverity::Error {
        return true;
    }
    match policy {
        BrepTopologyOutputPolicy::ClosedShell => severity != BrepTopologyDiagnosticSeverity::Info,
        BrepTopologyOutputPolicy::OpenSheet => !matches!(
            kind,
            BrepTopologyDiagnosticKind::DanglingBoundaryChain
                | BrepTopologyDiagnosticKind::FragmentedShell
        ),
    }
}

fn readiness_label(readiness: BrepExportReadiness) -> &'static str {
    match readiness {
        BrepExportReadiness::Ready => "ready",
        BrepExportReadiness::ReadyOpenSheet => "ready_open_sheet",
        BrepExportReadiness::Blocked => "blocked",
    }
}

fn topology_output_label(policy: BrepTopologyOutputPolicy) -> &'static str {
    match policy {
        BrepTopologyOutputPolicy::ClosedShell => "closed_shell",
        BrepTopologyOutputPolicy::OpenSheet => "open_sheet",
    }
}

fn collapsed_trim_missing_edge_curve_count(brep: &BRep, trim: BrepTrimRepairOptions) -> usize {
    let authoritative_edges = authoritative_degenerate_edges(brep, trim);
    brep.edges
        .iter()
        .filter(|(edge, data)| data.curve.is_none() && authoritative_edges.contains(edge))
        .count()
}

#[cfg(test)]
#[path = "brep_export_quality_tests.rs"]
mod tests;
