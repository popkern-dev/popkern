//! Candidate reports and diagnostics for multi-patch exact-boundary caps.
#![allow(dead_code)]

use popkern_core::topology::{BrepRepairFactId, EdgeHandle, FaceHandle, HalfEdgeHandle};

use crate::apex_construction::exact_boundary_cap::sector::{
    ExactBoundaryCapSectorInput, ExactBoundaryOuterSide,
};
use crate::types::{ApexConstructionDiagnostic, ExactBoundaryCapDiagnostic};

use super::inventory::{
    ExactBoundaryMultiPatchInventoryDiagnostic, ExactBoundaryMultiPatchInventoryReport,
    ExactBoundaryOuterStripPatchInventory,
};

const MULTIPATCH_STAGE: &str = "exact_boundary_cap.multipatch";

/// Deterministic multi-patch validation order, independent of construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum ExactBoundaryMultiPatchValidationOrder {
    /// Candidate patch/interface ownership is complete and source-backed.
    Inventory = 1,
    /// First implementation is limited to one selected outer side per sector.
    OuterSideSelection = 2,
    /// Outer support strip has a valid source pcurve parameterization.
    SourceParameterization = 3,
    /// Outer boundary position agrees with the immutable source edge.
    OuterBoundaryPosition = 4,
    /// Strip samples remain inside the support-normal tolerance band.
    SupportStripBand = 5,
    /// Strip/inner interface position continuity is within tolerance.
    InterfacePosition = 6,
    /// Strip/inner interface tangent continuity is within tolerance.
    InterfaceTangent = 7,
    /// Inner pole spokes meet at the sector apex within tolerance.
    InnerPolePosition = 8,
    /// Candidate patch domains remain orientation-preserving.
    Foldover = 9,
    /// Final report is internally consistent before routing can observe it.
    ReportConsistency = 10,
}

impl ExactBoundaryMultiPatchValidationOrder {
    /// Stable numeric order for UI sorting and deterministic skip policy.
    pub(crate) fn order_index(self) -> u8 {
        self as u8
    }
}

/// Stable multi-patch decline reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryMultiPatchDeclineReason {
    /// Candidate inventory is incomplete or disagrees with source ownership.
    InventoryMismatch,
    /// More outer sides were requested than this first routed solve supports.
    TwoOuterSidesUnsupported,
    /// Source pcurve parameterization is missing or invalid.
    SourceParameterizationMissing,
    /// Candidate strip does not satisfy the support-normal strip band.
    SupportStripBandUnsatisfied,
    /// Candidate strip/inner interface exceeds the position gate.
    InterfacePositionUnsatisfied,
    /// Candidate strip/inner interface exceeds the tangent gate.
    InterfaceTangentUnsatisfied,
    /// Inner pole spoke closure exceeds the position gate.
    InnerPolePositionUnsatisfied,
    /// Candidate patch domain folds over.
    PatchFoldover,
}

impl ExactBoundaryMultiPatchDeclineReason {
    /// Stable diagnostic code for API, GUI, tests, and routing policy.
    pub(crate) fn code(self) -> &'static str {
        match self {
            Self::InventoryMismatch => "exact_boundary_cap.multipatch.inventory_mismatch",
            Self::TwoOuterSidesUnsupported => {
                "exact_boundary_cap.multipatch.two_outer_sides_unsupported"
            }
            Self::SourceParameterizationMissing => {
                "exact_boundary_cap.multipatch.source_parameterization_missing"
            }
            Self::SupportStripBandUnsatisfied => {
                "exact_boundary_cap.multipatch.support_strip_band_unsatisfied"
            }
            Self::InterfacePositionUnsatisfied => {
                "exact_boundary_cap.multipatch.interface_position_unsatisfied"
            }
            Self::InterfaceTangentUnsatisfied => {
                "exact_boundary_cap.multipatch.interface_tangent_unsatisfied"
            }
            Self::InnerPolePositionUnsatisfied => {
                "exact_boundary_cap.multipatch.inner_pole_position_unsatisfied"
            }
            Self::PatchFoldover => "exact_boundary_cap.multipatch.patch_foldover",
        }
    }

    /// Deterministic gate associated with this decline.
    pub(crate) fn validation_order(self) -> ExactBoundaryMultiPatchValidationOrder {
        match self {
            Self::InventoryMismatch => ExactBoundaryMultiPatchValidationOrder::Inventory,
            Self::TwoOuterSidesUnsupported => {
                ExactBoundaryMultiPatchValidationOrder::OuterSideSelection
            }
            Self::SourceParameterizationMissing => {
                ExactBoundaryMultiPatchValidationOrder::SourceParameterization
            }
            Self::SupportStripBandUnsatisfied => {
                ExactBoundaryMultiPatchValidationOrder::SupportStripBand
            }
            Self::InterfacePositionUnsatisfied => {
                ExactBoundaryMultiPatchValidationOrder::InterfacePosition
            }
            Self::InterfaceTangentUnsatisfied => {
                ExactBoundaryMultiPatchValidationOrder::InterfaceTangent
            }
            Self::InnerPolePositionUnsatisfied => {
                ExactBoundaryMultiPatchValidationOrder::InnerPolePosition
            }
            Self::PatchFoldover => ExactBoundaryMultiPatchValidationOrder::Foldover,
        }
    }
}

/// Candidate outcome stored without mutating live topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryMultiPatchCandidateOutcome {
    /// Report has not yet reached an accept/decline decision.
    Pending,
    /// Candidate declined with at least one structured diagnostic.
    Declined,
}

/// Deterministic counters collected while preparing a multi-patch candidate.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) struct ExactBoundaryMultiPatchPerformanceReport {
    /// Candidate sectors considered.
    pub sector_count: usize,
    /// Candidate outer strip patches considered.
    pub outer_strip_patch_count: usize,
    /// Candidate strip/inner interfaces considered.
    pub interface_count: usize,
    /// Candidate inner pole patches considered.
    pub inner_patch_count: usize,
    /// Support strip samples evaluated.
    pub support_strip_sample_count: usize,
    /// Interface samples evaluated.
    pub interface_sample_count: usize,
    /// Bounded solver iterations consumed.
    pub bounded_iteration_count: usize,
    /// Elapsed candidate preparation and validation time in microseconds.
    pub elapsed_micros: u128,
    /// First deterministic decline code suitable for routing policy.
    pub deterministic_skip_code: Option<String>,
}

/// Candidate-only report for future routing and current diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ExactBoundaryMultiPatchCandidateReport {
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable source sector name.
    pub sector_name: String,
    /// Candidate outcome.
    pub outcome: ExactBoundaryMultiPatchCandidateOutcome,
    /// Optional validated inventory counts.
    pub inventory: Option<ExactBoundaryMultiPatchInventoryReport>,
    /// Structured diagnostics emitted by deterministic validation order.
    pub diagnostics: Vec<ExactBoundaryMultiPatchDiagnostic>,
    /// Performance counters for UI/API diagnostics and future hot-path gating.
    pub performance: ExactBoundaryMultiPatchPerformanceReport,
}

impl ExactBoundaryMultiPatchCandidateReport {
    /// Build a declined report for an inventory mismatch.
    pub(crate) fn declined_for_inventory_mismatch(
        diagnostic: ExactBoundaryMultiPatchInventoryDiagnostic,
    ) -> Self {
        let sector_index = diagnostic.sector_index;
        let sector_name = diagnostic.sector_name.clone();
        let diagnostic = ExactBoundaryMultiPatchDiagnostic::inventory_mismatch(diagnostic);
        Self::declined(sector_index, sector_name, diagnostic)
    }

    /// Build a declined report for unsupported two-side routing.
    pub(crate) fn declined_for_two_outer_sides(sector: &ExactBoundaryCapSectorInput) -> Self {
        Self::declined(
            sector.sector_index,
            sector.sector_name.clone(),
            ExactBoundaryMultiPatchDiagnostic::two_outer_sides_unsupported(sector),
        )
    }

    /// Build a declined report for a failed support-normal strip band.
    pub(crate) fn declined_for_support_strip_band(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        residual: f64,
        sample_parameter: f64,
        tolerance: f64,
        source_fact_ids: Vec<BrepRepairFactId>,
    ) -> Self {
        Self::declined(
            sector.sector_index,
            sector.sector_name.clone(),
            ExactBoundaryMultiPatchDiagnostic::support_strip_band(
                sector,
                strip,
                residual,
                sample_parameter,
                tolerance,
                source_fact_ids,
            ),
        )
    }

    fn declined(
        sector_index: usize,
        sector_name: String,
        diagnostic: ExactBoundaryMultiPatchDiagnostic,
    ) -> Self {
        let deterministic_skip_code = Some(diagnostic.code.clone());
        Self {
            sector_index,
            sector_name,
            outcome: ExactBoundaryMultiPatchCandidateOutcome::Declined,
            inventory: None,
            diagnostics: vec![diagnostic],
            performance: ExactBoundaryMultiPatchPerformanceReport {
                deterministic_skip_code,
                ..ExactBoundaryMultiPatchPerformanceReport::default()
            },
        }
    }
}

/// Structured multi-patch diagnostic with source and candidate highlighting.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ExactBoundaryMultiPatchDiagnostic {
    /// Stable decline reason.
    pub reason: ExactBoundaryMultiPatchDeclineReason,
    /// Stable diagnostic code.
    pub code: String,
    /// Optional detailed sub-code from the source subsystem.
    pub detail_code: Option<String>,
    /// Deterministic validation gate that emitted this diagnostic.
    pub validation_order: ExactBoundaryMultiPatchValidationOrder,
    /// Planner-order validation gate number.
    pub validation_order_index: u8,
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable source sector name.
    pub sector_name: String,
    /// Candidate patch name for GUI/API highlighting.
    pub patch_name: Option<String>,
    /// Candidate interface name for GUI/API highlighting.
    pub interface_name: Option<String>,
    /// Source outer side involved in the decline.
    pub side: Option<ExactBoundaryOuterSide>,
    /// Source exact-boundary edge.
    pub edge: Option<EdgeHandle>,
    /// Cap-side half-edge on the source boundary.
    pub cap_half_edge: Option<HalfEdgeHandle>,
    /// Support-side half-edge on the source boundary.
    pub support_half_edge: Option<HalfEdgeHandle>,
    /// Support face across the source boundary.
    pub support_face: Option<FaceHandle>,
    /// Measured residual for the failed gate.
    pub residual: Option<f64>,
    /// Local explicit tolerance used by the failed gate.
    pub tolerance: Option<f64>,
    /// Normalized sample parameter where the failed gate was measured.
    pub sample_parameter: Option<f64>,
    /// Shared source facts that explain the diagnostic evidence.
    pub source_fact_ids: Vec<BrepRepairFactId>,
}

impl ExactBoundaryMultiPatchDiagnostic {
    fn inventory_mismatch(diagnostic: ExactBoundaryMultiPatchInventoryDiagnostic) -> Self {
        Self::new(
            ExactBoundaryMultiPatchDeclineReason::InventoryMismatch,
            diagnostic.sector_index,
            diagnostic.sector_name,
        )
        .with_detail_code(diagnostic.code)
        .with_candidate_name(diagnostic.patch_or_interface_name)
        .with_side(diagnostic.side)
    }

    fn two_outer_sides_unsupported(sector: &ExactBoundaryCapSectorInput) -> Self {
        Self::new(
            ExactBoundaryMultiPatchDeclineReason::TwoOuterSidesUnsupported,
            sector.sector_index,
            sector.sector_name.clone(),
        )
        .with_edge(sector.face_edge)
        .with_cap_half_edge(sector.cap_half_edges[1])
        .with_support_half_edge(sector.support_half_edges[1])
        .with_support_face(sector.bridge_support_face)
    }

    fn support_strip_band(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        residual: f64,
        sample_parameter: f64,
        tolerance: f64,
        source_fact_ids: Vec<BrepRepairFactId>,
    ) -> Self {
        Self::new(
            ExactBoundaryMultiPatchDeclineReason::SupportStripBandUnsatisfied,
            sector.sector_index,
            sector.sector_name.clone(),
        )
        .with_patch_name(strip.patch_name.clone())
        .with_side(Some(strip.side))
        .with_edge(strip.edge)
        .with_cap_half_edge(strip.cap_half_edge)
        .with_support_half_edge(strip.support_half_edge)
        .with_support_face(strip.support_face)
        .with_residual(residual)
        .with_tolerance(tolerance)
        .with_sample_parameter(sample_parameter)
        .with_source_fact_ids(source_fact_ids)
    }

    fn new(
        reason: ExactBoundaryMultiPatchDeclineReason,
        sector_index: usize,
        sector_name: String,
    ) -> Self {
        let validation_order = reason.validation_order();
        Self {
            reason,
            code: reason.code().to_string(),
            detail_code: None,
            validation_order,
            validation_order_index: validation_order.order_index(),
            sector_index,
            sector_name,
            patch_name: None,
            interface_name: None,
            side: None,
            edge: None,
            cap_half_edge: None,
            support_half_edge: None,
            support_face: None,
            residual: None,
            tolerance: None,
            sample_parameter: None,
            source_fact_ids: Vec::new(),
        }
    }

    /// Convert into the existing apex diagnostic surface consumed downstream.
    pub(crate) fn into_apex_diagnostic(self) -> ApexConstructionDiagnostic {
        ApexConstructionDiagnostic::ExactBoundaryCap(ExactBoundaryCapDiagnostic {
            stage: MULTIPATCH_STAGE.to_string(),
            code: self.code,
            sector_name: self.sector_name,
            edge: self.edge,
            support_half_edge: self.support_half_edge,
            support_face: self.support_face,
            residual: self.residual,
            outer_side: self.side.map(|side| format!("{side:?}")),
            sample_parameter: self.sample_parameter,
            shared_fact_ids: self.source_fact_ids,
        })
    }

    fn with_detail_code(mut self, detail_code: &'static str) -> Self {
        self.detail_code = Some(detail_code.to_string());
        self
    }

    fn with_candidate_name(mut self, name: Option<String>) -> Self {
        let Some(name) = name else {
            return self;
        };
        if name.ends_with("_interface") {
            self.interface_name = Some(name);
        } else {
            self.patch_name = Some(name);
        }
        self
    }

    fn with_patch_name(mut self, name: String) -> Self {
        self.patch_name = Some(name);
        self
    }

    fn with_side(mut self, side: Option<ExactBoundaryOuterSide>) -> Self {
        self.side = side;
        self
    }

    fn with_edge(mut self, edge: EdgeHandle) -> Self {
        self.edge = Some(edge);
        self
    }

    fn with_cap_half_edge(mut self, half_edge: HalfEdgeHandle) -> Self {
        self.cap_half_edge = Some(half_edge);
        self
    }

    fn with_support_half_edge(mut self, half_edge: HalfEdgeHandle) -> Self {
        self.support_half_edge = Some(half_edge);
        self
    }

    fn with_support_face(mut self, face: FaceHandle) -> Self {
        self.support_face = Some(face);
        self
    }

    fn with_residual(mut self, residual: f64) -> Self {
        self.residual = Some(residual);
        self
    }

    fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = Some(tolerance);
        self
    }

    fn with_sample_parameter(mut self, sample_parameter: f64) -> Self {
        self.sample_parameter = Some(sample_parameter);
        self
    }

    fn with_source_fact_ids(mut self, source_fact_ids: Vec<BrepRepairFactId>) -> Self {
        self.source_fact_ids = source_fact_ids;
        self
    }
}

#[cfg(test)]
#[path = "multipatch_report_tests.rs"]
mod tests;
