use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::{ExactBoundaryMultiPatchInventory, ExactBoundaryMultiPatchInventoryDeclineReason};
use crate::apex_construction::exact_boundary_cap::sector::{
    ExactBoundaryCapSectorInput, ExactBoundaryOuterSeam, ExactBoundaryOuterSide,
};

#[test]
fn inventory_from_sector_preserves_source_ownership_without_brep_mutation() {
    let sector = sector_input();
    let inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FaceEdge,
    );
    let report = inventory
        .validate_against_sector(&sector)
        .expect("inventory validates");

    assert_eq!(report.actual_outer_strip_count, 1);
    assert_eq!(report.actual_interface_count, 1);
    assert_eq!(report.actual_inner_patch_count, 1);
    assert_eq!(inventory.sector_name, sector.sector_name);
    let strip = &inventory.outer_strips[0];
    assert_eq!(strip.patch_name, "direct_apex.cap_sector.0.face_edge_strip");
    assert_eq!(
        inventory.interfaces[0].interface_name,
        "direct_apex.cap_sector.0.face_edge_interface"
    );
    assert_eq!(
        inventory
            .inner_patch
            .as_ref()
            .map(|patch| patch.patch_name.as_str()),
        Some("direct_apex.cap_sector.0.inner_pole_patch")
    );
    assert_eq!(
        strip.edge,
        sector.bridge_seam.as_ref().expect("bridged fixture").edge
    );
    assert_eq!(
        strip.cap_half_edge,
        sector
            .bridge_seam
            .as_ref()
            .expect("bridged fixture")
            .cap_half_edge
    );
    assert_eq!(
        strip.support_half_edge,
        sector
            .bridge_seam
            .as_ref()
            .expect("bridged fixture")
            .support_half_edge
    );
    assert_eq!(
        strip.support_face,
        sector
            .bridge_seam
            .as_ref()
            .expect("bridged fixture")
            .support_face
    );
}

#[test]
fn inventory_detects_missing_strip_interface_and_inner_patch() {
    let sector = sector_input();
    let mut inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FilletEnd,
    );

    let mut missing_strip = inventory.clone();
    missing_strip.outer_strips.clear();
    assert_eq!(
        missing_strip
            .validate_against_sector(&sector)
            .expect_err("missing strip declines")
            .reason,
        ExactBoundaryMultiPatchInventoryDeclineReason::MissingOuterStrip
    );

    let mut missing_interface = inventory.clone();
    missing_interface.interfaces.clear();
    assert_eq!(
        missing_interface
            .validate_against_sector(&sector)
            .expect_err("missing interface declines")
            .reason,
        ExactBoundaryMultiPatchInventoryDeclineReason::MissingInterface
    );

    inventory.inner_patch = None;
    assert_eq!(
        inventory
            .validate_against_sector(&sector)
            .expect_err("missing inner patch declines")
            .reason,
        ExactBoundaryMultiPatchInventoryDeclineReason::MissingInnerPatch
    );
}

#[test]
fn inventory_detects_duplicate_names_and_support_ownership_mismatch() {
    let sector = sector_input();
    let mut inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FilletEnd,
    );
    inventory
        .inner_patch
        .as_mut()
        .expect("inner patch")
        .patch_name = inventory.outer_strips[0].patch_name.clone();
    assert_eq!(
        inventory
            .validate_against_sector(&sector)
            .expect_err("duplicate name declines")
            .reason,
        ExactBoundaryMultiPatchInventoryDeclineReason::DuplicatePatchName
    );

    let mut inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FilletEnd,
    );
    inventory.outer_strips[0].support_face = sector
        .bridge_seam
        .as_ref()
        .expect("bridged fixture")
        .support_face;
    assert_eq!(
        inventory
            .validate_against_sector(&sector)
            .expect_err("support mismatch declines")
            .reason,
        ExactBoundaryMultiPatchInventoryDeclineReason::MismatchedSupportOwnership
    );
}

#[test]
fn semantic_names_do_not_expose_planning_labels() {
    let sector = sector_input();
    let inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FilletEnd,
    );

    let mut names = vec![
        inventory.outer_strips[0].patch_name.as_str(),
        inventory.interfaces[0].interface_name.as_str(),
        inventory
            .inner_patch
            .as_ref()
            .expect("inner patch")
            .patch_name
            .as_str(),
    ];
    names.sort_unstable();
    for name in names {
        assert!(!name.contains("MP"));
        assert!(!name.contains("phase"));
    }
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
        face_edge: Some(edges[1]),
        fillet_support_face: faces[1],
        bridge_support_face: Some(faces[2]),
        fillet_seam: seam(
            ExactBoundaryOuterSide::FilletEnd,
            edges[0],
            half_edges[0],
            half_edges[1],
            faces[1],
        ),
        bridge_seam: Some(seam(
            ExactBoundaryOuterSide::FaceEdge,
            edges[1],
            half_edges[2],
            half_edges[3],
            faces[2],
        )),
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
