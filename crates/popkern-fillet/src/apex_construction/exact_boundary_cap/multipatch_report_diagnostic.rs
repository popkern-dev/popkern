use popkern_core::topology::{BrepRepairFactId, EdgeHandle, FaceHandle, HalfEdgeHandle};

use crate::apex_construction::exact_boundary_cap::sector::{
    ExactBoundaryCapSectorInput, ExactBoundaryOuterSide,
};
use crate::types::{ApexConstructionDiagnostic, ExactBoundaryCapDiagnostic};

use super::super::inventory::{
    ExactBoundaryMultiPatchInventoryDiagnostic, ExactBoundaryOuterStripPatchInventory,
    ExactBoundaryPatchInterfaceInventory,
};
use super::{ExactBoundaryMultiPatchDeclineReason, ExactBoundaryMultiPatchDiagnostic};

const MULTIPATCH_STAGE: &str = "exact_boundary_cap.multipatch";

impl ExactBoundaryMultiPatchDiagnostic {
    pub(crate) fn inventory_mismatch(
        diagnostic: ExactBoundaryMultiPatchInventoryDiagnostic,
    ) -> Self {
        Self::new(
            ExactBoundaryMultiPatchDeclineReason::InventoryMismatch,
            diagnostic.sector_index,
            diagnostic.sector_name,
        )
        .with_detail_code(diagnostic.code)
        .with_candidate_name(diagnostic.patch_or_interface_name)
        .with_side(diagnostic.side)
    }

    pub(super) fn two_outer_sides_unsupported(sector: &ExactBoundaryCapSectorInput) -> Self {
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

    pub(crate) fn pcurve_parameterization_invalid(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        residual: f64,
        sample_parameter: f64,
        tolerance: f64,
    ) -> Self {
        Self::strip_boundary_decline(
            ExactBoundaryMultiPatchDeclineReason::PcurveParameterizationInvalid,
            sector,
            strip,
            residual,
            sample_parameter,
            tolerance,
        )
    }

    pub(crate) fn rational_boundary_unsupported(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
    ) -> Self {
        Self::new(
            ExactBoundaryMultiPatchDeclineReason::RationalBoundaryUnsupported,
            sector.sector_index,
            sector.sector_name.clone(),
        )
        .with_patch_name(strip.patch_name.clone())
        .with_side(Some(strip.side))
        .with_edge(strip.edge)
        .with_cap_half_edge(strip.cap_half_edge)
        .with_support_half_edge(strip.support_half_edge)
        .with_support_face(strip.support_face)
    }

    pub(crate) fn outer_boundary_position(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        residual: f64,
        sample_parameter: f64,
        tolerance: f64,
    ) -> Self {
        Self::strip_boundary_decline(
            ExactBoundaryMultiPatchDeclineReason::OuterBoundaryPositionUnsatisfied,
            sector,
            strip,
            residual,
            sample_parameter,
            tolerance,
        )
    }

    pub(crate) fn patch_foldover(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        residual: f64,
        sample_parameter: f64,
        tolerance: f64,
    ) -> Self {
        Self::strip_boundary_decline(
            ExactBoundaryMultiPatchDeclineReason::PatchFoldover,
            sector,
            strip,
            residual,
            sample_parameter,
            tolerance,
        )
    }

    pub(crate) fn interface_position(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        interface: &ExactBoundaryPatchInterfaceInventory,
        residual: f64,
        sample_parameter: f64,
        tolerance: f64,
    ) -> Self {
        Self::new(
            ExactBoundaryMultiPatchDeclineReason::InterfacePositionUnsatisfied,
            sector.sector_index,
            sector.sector_name.clone(),
        )
        .with_patch_name(strip.patch_name.clone())
        .with_interface_name(interface.interface_name.clone())
        .with_side(Some(strip.side))
        .with_edge(strip.edge)
        .with_cap_half_edge(strip.cap_half_edge)
        .with_support_half_edge(strip.support_half_edge)
        .with_support_face(strip.support_face)
        .with_residual(residual)
        .with_tolerance(tolerance)
        .with_sample_parameter(sample_parameter)
    }

    pub(crate) fn interface_tangent(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        interface: &ExactBoundaryPatchInterfaceInventory,
        residual: f64,
        sample_parameter: f64,
        tolerance: f64,
    ) -> Self {
        Self::interface_decline(
            ExactBoundaryMultiPatchDeclineReason::InterfaceTangentUnsatisfied,
            sector,
            strip,
            interface,
            residual,
            sample_parameter,
            tolerance,
        )
    }

    pub(crate) fn pole_position(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        interface: &ExactBoundaryPatchInterfaceInventory,
        residual: f64,
        sample_parameter: f64,
        tolerance: f64,
    ) -> Self {
        Self::interface_decline(
            ExactBoundaryMultiPatchDeclineReason::InnerPolePositionUnsatisfied,
            sector,
            strip,
            interface,
            residual,
            sample_parameter,
            tolerance,
        )
        .with_patch_name(interface.inner_patch_name.clone())
    }

    pub(crate) fn support_strip_band(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        residual: f64,
        sample_parameter: f64,
        tolerance: f64,
        source_fact_ids: Vec<BrepRepairFactId>,
    ) -> Self {
        Self::strip_boundary_decline(
            ExactBoundaryMultiPatchDeclineReason::SupportStripBandUnsatisfied,
            sector,
            strip,
            residual,
            sample_parameter,
            tolerance,
        )
        .with_source_fact_ids(source_fact_ids)
    }

    pub(crate) fn commit_preflight_missing_surface(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        patch_name: String,
    ) -> Self {
        Self::commit_preflight_decline(
            ExactBoundaryMultiPatchDeclineReason::CommitPreflightMissingSurface,
            sector,
            strip,
        )
        .with_patch_name(patch_name)
    }

    pub(crate) fn commit_preflight_missing_pcurve(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        patch_name: String,
        interface_name: Option<String>,
    ) -> Self {
        let diagnostic = Self::commit_preflight_decline(
            ExactBoundaryMultiPatchDeclineReason::CommitPreflightMissingPcurve,
            sector,
            strip,
        )
        .with_patch_name(patch_name);
        if let Some(interface_name) = interface_name {
            diagnostic.with_interface_name(interface_name)
        } else {
            diagnostic
        }
    }

    pub(crate) fn commit_preflight_missing_interface(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        interface_name: String,
    ) -> Self {
        Self::commit_preflight_decline(
            ExactBoundaryMultiPatchDeclineReason::CommitPreflightMissingInterface,
            sector,
            strip,
        )
        .with_interface_name(interface_name)
    }

    pub(crate) fn commit_preflight_missing_generated_spoke(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        patch_name: String,
    ) -> Self {
        Self::commit_preflight_decline(
            ExactBoundaryMultiPatchDeclineReason::CommitPreflightMissingGeneratedSpoke,
            sector,
            strip,
        )
        .with_patch_name(patch_name)
    }

    pub(crate) fn commit_materialization_failed(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        detail_code: &'static str,
        residual: Option<f64>,
    ) -> Self {
        Self::commit_preflight_decline(
            ExactBoundaryMultiPatchDeclineReason::CommitMaterializationFailed,
            sector,
            strip,
        )
        .with_patch_name(strip.patch_name.clone())
        .with_detail_code(detail_code)
        .with_optional_residual(residual)
    }

    pub(crate) fn commit_validation_failed(
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        residual: Option<f64>,
    ) -> Self {
        Self::commit_preflight_decline(
            ExactBoundaryMultiPatchDeclineReason::CommitValidationFailed,
            sector,
            strip,
        )
        .with_patch_name(strip.patch_name.clone())
        .with_optional_residual(residual)
    }

    fn strip_boundary_decline(
        reason: ExactBoundaryMultiPatchDeclineReason,
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        residual: f64,
        sample_parameter: f64,
        tolerance: f64,
    ) -> Self {
        Self::new(reason, sector.sector_index, sector.sector_name.clone())
            .with_patch_name(strip.patch_name.clone())
            .with_side(Some(strip.side))
            .with_edge(strip.edge)
            .with_cap_half_edge(strip.cap_half_edge)
            .with_support_half_edge(strip.support_half_edge)
            .with_support_face(strip.support_face)
            .with_residual(residual)
            .with_tolerance(tolerance)
            .with_sample_parameter(sample_parameter)
    }

    fn commit_preflight_decline(
        reason: ExactBoundaryMultiPatchDeclineReason,
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
    ) -> Self {
        Self::new(reason, sector.sector_index, sector.sector_name.clone())
            .with_side(Some(strip.side))
            .with_edge(strip.edge)
            .with_cap_half_edge(strip.cap_half_edge)
            .with_support_half_edge(strip.support_half_edge)
            .with_support_face(strip.support_face)
    }

    fn interface_decline(
        reason: ExactBoundaryMultiPatchDeclineReason,
        sector: &ExactBoundaryCapSectorInput,
        strip: &ExactBoundaryOuterStripPatchInventory,
        interface: &ExactBoundaryPatchInterfaceInventory,
        residual: f64,
        sample_parameter: f64,
        tolerance: f64,
    ) -> Self {
        Self::new(reason, sector.sector_index, sector.sector_name.clone())
            .with_patch_name(strip.patch_name.clone())
            .with_interface_name(interface.interface_name.clone())
            .with_side(Some(strip.side))
            .with_edge(strip.edge)
            .with_cap_half_edge(strip.cap_half_edge)
            .with_support_half_edge(strip.support_half_edge)
            .with_support_face(strip.support_face)
            .with_residual(residual)
            .with_tolerance(tolerance)
            .with_sample_parameter(sample_parameter)
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
            tolerance: self.tolerance,
            shared_fact_ids: self.source_fact_ids,
        })
    }

    fn with_detail_code(mut self, detail_code: &'static str) -> Self {
        self.detail_code = Some(detail_code.to_string());
        self
    }

    fn with_optional_residual(mut self, residual: Option<f64>) -> Self {
        self.residual = residual;
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

    fn with_interface_name(mut self, name: String) -> Self {
        self.interface_name = Some(name);
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
