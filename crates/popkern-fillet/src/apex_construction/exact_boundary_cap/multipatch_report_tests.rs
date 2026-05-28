use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BrepRepairFactId, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::{
    ExactBoundaryMultiPatchCandidateOutcome, ExactBoundaryMultiPatchCandidateReport,
    ExactBoundaryMultiPatchDeclineReason, ExactBoundaryMultiPatchValidationOrder,
};
use crate::apex_construction::exact_boundary_cap::multipatch::inventory::{
    ExactBoundaryMultiPatchInventory, ExactBoundaryMultiPatchInventoryDeclineReason,
};
use crate::apex_construction::exact_boundary_cap::sector::{
    ExactBoundaryCapSectorInput, ExactBoundaryOuterSeam, ExactBoundaryOuterSide,
};
use crate::types::ApexConstructionDiagnostic;

#[test]
fn report_promotes_inventory_mismatch_without_live_brep_mutation() {
    let sector = sector_input();
    let mut inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FaceEdge,
    );
    inventory.outer_strips.clear();
    let inventory_diagnostic = inventory
        .validate_against_sector(&sector)
        .expect_err("missing outer strip declines");

    let report = ExactBoundaryMultiPatchCandidateReport::declined_for_inventory_mismatch(
        inventory_diagnostic,
    );
    let diagnostic = &report.diagnostics[0];

    assert_eq!(
        report.outcome,
        ExactBoundaryMultiPatchCandidateOutcome::Declined
    );
    assert_eq!(
        report.performance.deterministic_skip_code.as_deref(),
        Some("exact_boundary_cap.multipatch.inventory_mismatch")
    );
    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::InventoryMismatch
    );
    assert_eq!(
        diagnostic.detail_code.as_deref(),
        Some(ExactBoundaryMultiPatchInventoryDeclineReason::MissingOuterStrip.code())
    );
    assert_eq!(
        diagnostic.validation_order,
        ExactBoundaryMultiPatchValidationOrder::Inventory
    );
    assert_eq!(diagnostic.validation_order_index, 1);
    assert_eq!(diagnostic.sector_name, "direct_apex.cap_sector.0");
}

#[test]
fn inventory_report_keeps_candidate_interface_highlight_separate() {
    let sector = sector_input();
    let mut inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FaceEdge,
    );
    let interface_name = inventory.interfaces[0].interface_name.clone();
    inventory
        .inner_patch
        .as_mut()
        .expect("inner patch")
        .patch_name = interface_name.clone();
    let inventory_diagnostic = inventory
        .validate_against_sector(&sector)
        .expect_err("duplicate interface name declines");

    let report = ExactBoundaryMultiPatchCandidateReport::declined_for_inventory_mismatch(
        inventory_diagnostic,
    );
    let diagnostic = &report.diagnostics[0];

    assert_eq!(diagnostic.patch_name, None);
    assert_eq!(
        diagnostic.interface_name.as_deref(),
        Some(interface_name.as_str())
    );
}

#[test]
fn report_describes_two_outer_side_limit_for_gui_and_api() {
    let sector = sector_input();
    let report = ExactBoundaryMultiPatchCandidateReport::declined_for_two_outer_sides(&sector);
    let diagnostic = &report.diagnostics[0];

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::TwoOuterSidesUnsupported
    );
    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.multipatch.two_outer_sides_unsupported"
    );
    assert_eq!(
        diagnostic.validation_order,
        ExactBoundaryMultiPatchValidationOrder::OuterSideSelection
    );
    assert_eq!(diagnostic.edge, Some(sector.face_edge));
    assert_eq!(diagnostic.cap_half_edge, Some(sector.cap_half_edges[1]));
    assert_eq!(
        diagnostic.support_half_edge,
        Some(sector.support_half_edges[1])
    );
    assert_eq!(diagnostic.support_face, Some(sector.bridge_support_face));
}

#[test]
fn support_strip_band_report_carries_source_and_candidate_highlights() {
    let sector = sector_input();
    let inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FaceEdge,
    );
    let strip = &inventory.outer_strips[0];
    let report = ExactBoundaryMultiPatchCandidateReport::declined_for_support_strip_band(
        &sector,
        strip,
        0.0625,
        0.75,
        0.005,
        vec![BrepRepairFactId(7), BrepRepairFactId(9)],
    );
    let diagnostic = &report.diagnostics[0];

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::SupportStripBandUnsatisfied
    );
    assert_eq!(
        diagnostic.validation_order,
        ExactBoundaryMultiPatchValidationOrder::SupportStripBand
    );
    assert_eq!(
        diagnostic.patch_name.as_deref(),
        Some(strip.patch_name.as_str())
    );
    assert_eq!(diagnostic.side, Some(ExactBoundaryOuterSide::FaceEdge));
    assert_eq!(diagnostic.edge, Some(strip.edge));
    assert_eq!(diagnostic.cap_half_edge, Some(strip.cap_half_edge));
    assert_eq!(diagnostic.support_half_edge, Some(strip.support_half_edge));
    assert_eq!(diagnostic.support_face, Some(strip.support_face));
    assert_eq!(diagnostic.residual, Some(0.0625));
    assert_eq!(diagnostic.tolerance, Some(0.005));
    assert_eq!(diagnostic.sample_parameter, Some(0.75));
    assert_eq!(
        diagnostic.source_fact_ids,
        vec![BrepRepairFactId(7), BrepRepairFactId(9)]
    );
}

#[test]
fn multipatch_diagnostic_converts_to_existing_apex_diagnostic_surface() {
    let sector = sector_input();
    let inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FaceEdge,
    );
    let report = ExactBoundaryMultiPatchCandidateReport::declined_for_support_strip_band(
        &sector,
        &inventory.outer_strips[0],
        0.04,
        0.25,
        0.005,
        vec![BrepRepairFactId(3)],
    );

    let ApexConstructionDiagnostic::ExactBoundaryCap(apex_diagnostic) =
        report.diagnostics[0].clone().into_apex_diagnostic()
    else {
        panic!("multi-patch report emits exact-boundary cap diagnostics");
    };

    assert_eq!(apex_diagnostic.stage, "exact_boundary_cap.multipatch");
    assert_eq!(
        apex_diagnostic.code,
        "exact_boundary_cap.multipatch.support_strip_band_unsatisfied"
    );
    assert_eq!(apex_diagnostic.outer_side.as_deref(), Some("FaceEdge"));
    assert_eq!(apex_diagnostic.shared_fact_ids, vec![BrepRepairFactId(3)]);
    assert_eq!(apex_diagnostic.tolerance, Some(0.005));
}

fn sector_input() -> ExactBoundaryCapSectorInput {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube");
    let edges: Vec<EdgeHandle> = brep.edges.keys().collect();
    let half_edges: Vec<HalfEdgeHandle> = brep.half_edges.keys().collect();
    let faces: Vec<FaceHandle> = brep.faces.keys().collect();
    ExactBoundaryCapSectorInput {
        sector_index: 0,
        sector_name: "direct_apex.cap_sector.0".to_string(),
        fillet_index: 0,
        cap_face: faces[0],
        fillet_end_edge: edges[0],
        face_edge: edges[1],
        fillet_support_face: faces[1],
        bridge_support_face: faces[2],
        cap_half_edges: [half_edges[0], half_edges[2]],
        support_half_edges: [half_edges[1], half_edges[3]],
        outer_seams: [
            seam(
                ExactBoundaryOuterSide::FilletEnd,
                edges[0],
                half_edges[0],
                half_edges[1],
                faces[1],
            ),
            seam(
                ExactBoundaryOuterSide::FaceEdge,
                edges[1],
                half_edges[2],
                half_edges[3],
                faces[2],
            ),
        ],
    }
}

fn seam(
    side: ExactBoundaryOuterSide,
    edge: EdgeHandle,
    cap_half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
    support_face: FaceHandle,
) -> ExactBoundaryOuterSeam {
    ExactBoundaryOuterSeam {
        side,
        edge,
        cap_half_edge,
        support_half_edge,
        support_face,
        cap_pcurve_present: true,
        support_pcurve_present: true,
    }
}
