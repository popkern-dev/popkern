//! Candidate ownership inventory for multi-patch exact-boundary caps.
#![allow(dead_code)]

use std::collections::BTreeSet;

use popkern_core::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle};

use crate::apex_construction::exact_boundary_cap::sector::{
    ExactBoundaryCapSectorInput, ExactBoundaryOuterSeam, ExactBoundaryOuterSide,
};

/// Candidate ownership inventory for one sector split into strip/inner patches.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExactBoundaryMultiPatchInventory {
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable sector name copied from the source exact-boundary sector.
    pub sector_name: String,
    /// Expected number of outer strip patches for the first implementation.
    pub expected_outer_strip_count: usize,
    /// Expected number of strip/inner interfaces for the first implementation.
    pub expected_interface_count: usize,
    /// Expected number of inner pole patches for the first implementation.
    pub expected_inner_patch_count: usize,
    /// Candidate-owned outer strip patch inventories.
    pub outer_strips: Vec<ExactBoundaryOuterStripPatchInventory>,
    /// Candidate-owned inner pole patch inventory.
    pub inner_patch: Option<ExactBoundaryInnerPolePatchInventory>,
    /// Candidate-owned interface inventories.
    pub interfaces: Vec<ExactBoundaryPatchInterfaceInventory>,
}

impl ExactBoundaryMultiPatchInventory {
    /// Build first-slice inventory for exactly one selected failing outer side.
    pub(crate) fn for_selected_side(
        sector: &ExactBoundaryCapSectorInput,
        side: ExactBoundaryOuterSide,
    ) -> Self {
        let seam = outer_seam(sector, side);
        let strip = ExactBoundaryOuterStripPatchInventory::from_seam(sector, seam);
        let inner_patch = ExactBoundaryInnerPolePatchInventory::new(sector);
        let interface = ExactBoundaryPatchInterfaceInventory::new(
            sector,
            side,
            strip.patch_name.clone(),
            inner_patch.patch_name.clone(),
        );
        Self {
            sector_index: sector.sector_index,
            sector_name: sector.sector_name.clone(),
            expected_outer_strip_count: 1,
            expected_interface_count: 1,
            expected_inner_patch_count: 1,
            outer_strips: vec![strip],
            inner_patch: Some(inner_patch),
            interfaces: vec![interface],
        }
    }

    /// Validate candidate ownership without touching live BRep topology.
    pub(crate) fn validate_against_sector(
        &self,
        sector: &ExactBoundaryCapSectorInput,
    ) -> Result<ExactBoundaryMultiPatchInventoryReport, ExactBoundaryMultiPatchInventoryDiagnostic>
    {
        if self.outer_strips.len() != self.expected_outer_strip_count
            || self.outer_strips.is_empty()
        {
            return Err(self.diagnostic(
                ExactBoundaryMultiPatchInventoryDeclineReason::MissingOuterStrip,
                None,
                None,
            ));
        }
        if self.interfaces.len() != self.expected_interface_count || self.interfaces.is_empty() {
            return Err(self.diagnostic(
                ExactBoundaryMultiPatchInventoryDeclineReason::MissingInterface,
                None,
                None,
            ));
        }
        if self.expected_inner_patch_count == 1 && self.inner_patch.is_none() {
            return Err(self.diagnostic(
                ExactBoundaryMultiPatchInventoryDeclineReason::MissingInnerPatch,
                None,
                None,
            ));
        }
        if let Some(name) = duplicate_name(self) {
            return Err(self.diagnostic(
                ExactBoundaryMultiPatchInventoryDeclineReason::DuplicatePatchName,
                Some(name),
                None,
            ));
        }
        for strip in &self.outer_strips {
            let seam = outer_seam(sector, strip.side);
            if strip.edge != seam.edge
                || strip.cap_half_edge != seam.cap_half_edge
                || strip.support_half_edge != seam.support_half_edge
                || strip.support_face != seam.support_face
            {
                return Err(self.diagnostic(
                    ExactBoundaryMultiPatchInventoryDeclineReason::MismatchedSupportOwnership,
                    Some(strip.patch_name.clone()),
                    Some(strip.side),
                ));
            }
        }
        Ok(ExactBoundaryMultiPatchInventoryReport {
            sector_index: self.sector_index,
            sector_name: self.sector_name.clone(),
            actual_outer_strip_count: self.outer_strips.len(),
            actual_interface_count: self.interfaces.len(),
            actual_inner_patch_count: usize::from(self.inner_patch.is_some()),
        })
    }

    fn diagnostic(
        &self,
        reason: ExactBoundaryMultiPatchInventoryDeclineReason,
        patch_or_interface_name: Option<String>,
        side: Option<ExactBoundaryOuterSide>,
    ) -> ExactBoundaryMultiPatchInventoryDiagnostic {
        ExactBoundaryMultiPatchInventoryDiagnostic {
            reason,
            code: reason.code(),
            sector_index: self.sector_index,
            sector_name: self.sector_name.clone(),
            patch_or_interface_name,
            side,
        }
    }
}

/// Candidate-only outer support strip patch ownership.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExactBoundaryOuterStripPatchInventory {
    /// Stable semantic patch name.
    pub patch_name: String,
    /// Source outer side owned by this strip.
    pub side: ExactBoundaryOuterSide,
    /// Source exact-boundary edge.
    pub edge: EdgeHandle,
    /// Cap-side half-edge for the source exact-boundary edge.
    pub cap_half_edge: HalfEdgeHandle,
    /// Support-side half-edge for the source exact-boundary edge.
    pub support_half_edge: HalfEdgeHandle,
    /// Support face across the source exact-boundary edge.
    pub support_face: FaceHandle,
}

impl ExactBoundaryOuterStripPatchInventory {
    fn from_seam(sector: &ExactBoundaryCapSectorInput, seam: &ExactBoundaryOuterSeam) -> Self {
        Self {
            patch_name: outer_strip_patch_name(&sector.sector_name, seam.side),
            side: seam.side,
            edge: seam.edge,
            cap_half_edge: seam.cap_half_edge,
            support_half_edge: seam.support_half_edge,
            support_face: seam.support_face,
        }
    }
}

/// Candidate-only inner pole patch ownership.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExactBoundaryInnerPolePatchInventory {
    /// Stable semantic patch name.
    pub patch_name: String,
}

impl ExactBoundaryInnerPolePatchInventory {
    fn new(sector: &ExactBoundaryCapSectorInput) -> Self {
        Self {
            patch_name: format!("{}.inner_pole_patch", sector.sector_name),
        }
    }
}

/// Candidate-only interface ownership between one strip and the inner patch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExactBoundaryPatchInterfaceInventory {
    /// Stable semantic interface name.
    pub interface_name: String,
    /// Source side whose strip owns this interface.
    pub side: ExactBoundaryOuterSide,
    /// Adjacent outer strip patch name.
    pub strip_patch_name: String,
    /// Adjacent inner pole patch name.
    pub inner_patch_name: String,
}

impl ExactBoundaryPatchInterfaceInventory {
    fn new(
        sector: &ExactBoundaryCapSectorInput,
        side: ExactBoundaryOuterSide,
        strip_patch_name: String,
        inner_patch_name: String,
    ) -> Self {
        Self {
            interface_name: interface_name(&sector.sector_name, side),
            side,
            strip_patch_name,
            inner_patch_name,
        }
    }
}

/// Validated candidate ownership counts for reports and later routing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExactBoundaryMultiPatchInventoryReport {
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable source sector name.
    pub sector_name: String,
    /// Actual outer strip count.
    pub actual_outer_strip_count: usize,
    /// Actual interface count.
    pub actual_interface_count: usize,
    /// Actual inner patch count.
    pub actual_inner_patch_count: usize,
}

/// Stable inventory decline reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryMultiPatchInventoryDeclineReason {
    /// Expected outer strip patch is missing.
    MissingOuterStrip,
    /// Expected strip/inner interface is missing.
    MissingInterface,
    /// Expected inner pole patch is missing.
    MissingInnerPatch,
    /// Patch or interface names are not unique.
    DuplicatePatchName,
    /// Candidate strip ownership disagrees with source support seam ownership.
    MismatchedSupportOwnership,
}

impl ExactBoundaryMultiPatchInventoryDeclineReason {
    /// Stable diagnostic code for exact-boundary multi-patch reports.
    pub(crate) fn code(self) -> &'static str {
        match self {
            Self::MissingOuterStrip => {
                "exact_boundary_cap.multipatch.inventory.missing_outer_strip"
            }
            Self::MissingInterface => "exact_boundary_cap.multipatch.inventory.missing_interface",
            Self::MissingInnerPatch => {
                "exact_boundary_cap.multipatch.inventory.missing_inner_patch"
            }
            Self::DuplicatePatchName => "exact_boundary_cap.multipatch.inventory.duplicate_name",
            Self::MismatchedSupportOwnership => {
                "exact_boundary_cap.multipatch.inventory.mismatched_support_ownership"
            }
        }
    }
}

/// Structured inventory diagnostic with source sector/side identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ExactBoundaryMultiPatchInventoryDiagnostic {
    /// Stable decline reason.
    pub reason: ExactBoundaryMultiPatchInventoryDeclineReason,
    /// Stable diagnostic code.
    pub code: &'static str,
    /// Planner-order sector index.
    pub sector_index: usize,
    /// Stable source sector name.
    pub sector_name: String,
    /// Patch or interface name involved in the decline.
    pub patch_or_interface_name: Option<String>,
    /// Source side involved in the decline.
    pub side: Option<ExactBoundaryOuterSide>,
}

fn outer_seam(
    sector: &ExactBoundaryCapSectorInput,
    side: ExactBoundaryOuterSide,
) -> &ExactBoundaryOuterSeam {
    match side {
        ExactBoundaryOuterSide::FilletEnd => &sector.outer_seams[0],
        ExactBoundaryOuterSide::FaceEdge => &sector.outer_seams[1],
    }
}

fn duplicate_name(inventory: &ExactBoundaryMultiPatchInventory) -> Option<String> {
    let mut names = BTreeSet::new();
    for strip in &inventory.outer_strips {
        if !names.insert(strip.patch_name.as_str()) {
            return Some(strip.patch_name.clone());
        }
    }
    if let Some(inner) = &inventory.inner_patch
        && !names.insert(inner.patch_name.as_str())
    {
        return Some(inner.patch_name.clone());
    }
    for interface in &inventory.interfaces {
        if !names.insert(interface.interface_name.as_str()) {
            return Some(interface.interface_name.clone());
        }
    }
    None
}

fn outer_strip_patch_name(sector_name: &str, side: ExactBoundaryOuterSide) -> String {
    format!("{sector_name}.{}_strip", side_name(side))
}

fn interface_name(sector_name: &str, side: ExactBoundaryOuterSide) -> String {
    format!("{sector_name}.{}_interface", side_name(side))
}

fn side_name(side: ExactBoundaryOuterSide) -> &'static str {
    match side {
        ExactBoundaryOuterSide::FilletEnd => "fillet_end",
        ExactBoundaryOuterSide::FaceEdge => "face_edge",
    }
}

#[cfg(test)]
#[path = "multipatch_inventory_tests.rs"]
mod tests;
