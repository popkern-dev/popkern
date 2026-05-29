use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::{validate_inner_pole_patch_prototype, validate_outer_strip_candidate};
use crate::apex_construction::exact_boundary_cap::multipatch::budget::ExactBoundaryMultiPatchBudget;
use crate::apex_construction::exact_boundary_cap::multipatch::inner::{
    construct_inner_pole_patch_prototype, construct_strip_inner_interface_curve,
    measure_inner_interface_g1,
};
use crate::apex_construction::exact_boundary_cap::multipatch::inventory::{
    ExactBoundaryMultiPatchInventory, ExactBoundaryMultiPatchInventoryDeclineReason,
};
use crate::apex_construction::exact_boundary_cap::multipatch::report::{
    ExactBoundaryMultiPatchDeclineReason, ExactBoundaryMultiPatchValidationOrder,
};
use crate::apex_construction::exact_boundary_cap::multipatch::strip::{
    ExactBoundaryStripBoundaryDirection, scaffold_outer_strip_boundary,
    solve_strip_support_normal_row,
};
use crate::apex_construction::exact_boundary_cap::ribbon::ExactBoundaryRibbonSample;
use crate::apex_construction::exact_boundary_cap::sector::{
    ExactBoundaryCapSectorInput, ExactBoundaryOuterSeam, ExactBoundaryOuterSide,
};

#[test]
fn strip_validation_accepts_solved_candidate_with_deterministic_order() {
    let fixture = solved_fixture(Vec3::new(0.0, 0.0, 1.0));
    let report = validate_outer_strip_candidate(
        &fixture.sector,
        &fixture.strip,
        None,
        &fixture.candidate,
        &fixture.boundary,
        &fixture.samples,
        fixture.budget,
    )
    .expect("solved strip validates");

    assert_eq!(
        report.validation_order,
        vec![
            ExactBoundaryMultiPatchValidationOrder::Inventory,
            ExactBoundaryMultiPatchValidationOrder::SourceParameterization,
            ExactBoundaryMultiPatchValidationOrder::OuterBoundaryPosition,
            ExactBoundaryMultiPatchValidationOrder::SupportStripBand,
            ExactBoundaryMultiPatchValidationOrder::Foldover,
            ExactBoundaryMultiPatchValidationOrder::ReportConsistency,
        ]
    );
    assert!(report.max_outer_boundary_residual <= fixture.budget.outer_boundary_residual);
    assert!(report.max_support_normal_deficit <= fixture.budget.support_normal_deficit);
    assert!(report.min_normal_length >= fixture.budget.min_normal_length);
}

#[test]
fn strip_validation_reports_inventory_mismatch_first() {
    let mut fixture = solved_fixture(Vec3::new(0.0, 0.0, 1.0));
    let mut inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &fixture.sector,
        ExactBoundaryOuterSide::FaceEdge,
    );
    inventory.outer_strips.clear();
    let inventory_diagnostic = inventory
        .validate_against_sector(&fixture.sector)
        .expect_err("inventory mismatch");
    fixture.candidate.surface.control_points[0].z = 0.5;

    let diagnostic = validate_outer_strip_candidate(
        &fixture.sector,
        &fixture.strip,
        Some(inventory_diagnostic),
        &fixture.candidate,
        &fixture.boundary,
        &fixture.samples,
        fixture.budget,
    )
    .expect_err("inventory mismatch wins");

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
}

#[test]
fn strip_validation_reports_pcurve_invalid_before_outer_residual() {
    let mut fixture = solved_fixture(Vec3::new(0.0, 0.0, 1.0));
    fixture.candidate.outer_pcurve =
        popkern_core::geometry::pcurve::Pcurve::linear((0.0, 0.0), (0.0, 0.0), 0.0, 1.0);
    fixture.candidate.surface.control_points[0].z = 0.5;

    let diagnostic = validate_outer_strip_candidate(
        &fixture.sector,
        &fixture.strip,
        None,
        &fixture.candidate,
        &fixture.boundary,
        &fixture.samples,
        fixture.budget,
    )
    .expect_err("invalid pcurve wins");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::PcurveParameterizationInvalid
    );
    assert_eq!(
        diagnostic.validation_order,
        ExactBoundaryMultiPatchValidationOrder::SourceParameterization
    );
}

#[test]
fn strip_validation_reports_outer_residual_before_support_band() {
    let mut fixture = solved_fixture(Vec3::new(1.0, 0.0, 0.0));
    fixture.candidate.surface.control_points[0].z = 0.5;

    let diagnostic = validate_outer_strip_candidate(
        &fixture.sector,
        &fixture.strip,
        None,
        &fixture.candidate,
        &fixture.boundary,
        &fixture.samples,
        fixture.budget,
    )
    .expect_err("outer residual wins before support band");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::OuterBoundaryPositionUnsatisfied
    );
    assert_eq!(
        diagnostic.validation_order,
        ExactBoundaryMultiPatchValidationOrder::OuterBoundaryPosition
    );
    assert!(diagnostic.residual.expect("residual") > fixture.budget.outer_boundary_residual);
}

#[test]
fn strip_validation_reports_support_strip_band_with_source_handles() {
    let fixture = solved_fixture(Vec3::new(1.0, 0.0, 0.0));

    let diagnostic = validate_outer_strip_candidate(
        &fixture.sector,
        &fixture.strip,
        None,
        &fixture.candidate,
        &fixture.boundary,
        &fixture.samples,
        fixture.budget,
    )
    .expect_err("support band fails");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::SupportStripBandUnsatisfied
    );
    assert_eq!(
        diagnostic.support_half_edge,
        Some(fixture.strip.support_half_edge)
    );
    assert_eq!(diagnostic.support_face, Some(fixture.strip.support_face));
    assert_eq!(
        diagnostic.tolerance,
        Some(fixture.budget.support_normal_deficit)
    );
}

#[test]
fn strip_validation_reports_foldover_after_support_band_passes() {
    let mut fixture = solved_fixture(Vec3::new(0.0, 0.0, 1.0));
    fixture.samples.clear();
    for index in 0..fixture.candidate.surface.num_u {
        let row0 =
            fixture.candidate.surface.control_points[index * fixture.candidate.surface.num_v];
        for row in 1..fixture.candidate.surface.num_v {
            fixture.candidate.surface.control_points
                [index * fixture.candidate.surface.num_v + row] = row0;
        }
    }

    let diagnostic = validate_outer_strip_candidate(
        &fixture.sector,
        &fixture.strip,
        None,
        &fixture.candidate,
        &fixture.boundary,
        &fixture.samples,
        fixture.budget,
    )
    .expect_err("foldover is detected");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::PatchFoldover
    );
    assert_eq!(
        diagnostic.validation_order,
        ExactBoundaryMultiPatchValidationOrder::Foldover
    );
}

#[test]
fn inner_pole_validation_declines_g0_only_candidate_before_commit() {
    let fixture = solved_fixture(Vec3::new(0.0, 0.0, 1.0));
    let inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &fixture.sector,
        ExactBoundaryOuterSide::FaceEdge,
    );
    let interface = construct_strip_inner_interface_curve(
        &fixture.sector,
        &fixture.strip,
        &inventory.interfaces[0],
        &fixture.candidate,
        fixture.budget.interface_g0,
    )
    .expect("interface constructs");
    let prototype = construct_inner_pole_patch_prototype(
        &fixture.sector,
        &fixture.strip,
        &inventory.interfaces[0],
        &fixture.candidate,
        &interface,
        Point3::new(0.5, 0.7, 0.0),
        fixture.budget.pole_g0,
    )
    .expect("inner G0 prototype constructs");

    let diagnostic = validate_inner_pole_patch_prototype(
        &fixture.sector,
        &fixture.strip,
        &inventory.interfaces[0],
        &interface,
        &prototype,
        fixture.budget,
    )
    .expect_err("G0-only prototype declines before commit");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::InterfaceTangentUnsatisfied
    );
    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.multipatch.interface_g1_exceeded"
    );
    assert_eq!(
        diagnostic.interface_name.as_deref(),
        Some(inventory.interfaces[0].interface_name.as_str())
    );
}

#[test]
fn inner_pole_validation_accepts_candidate_after_interface_g1_measurement() {
    let fixture = solved_fixture(Vec3::new(0.0, 0.0, 1.0));
    let inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &fixture.sector,
        ExactBoundaryOuterSide::FaceEdge,
    );
    let interface = construct_strip_inner_interface_curve(
        &fixture.sector,
        &fixture.strip,
        &inventory.interfaces[0],
        &fixture.candidate,
        fixture.budget.interface_g0,
    )
    .expect("interface constructs");
    let prototype = construct_inner_pole_patch_prototype(
        &fixture.sector,
        &fixture.strip,
        &inventory.interfaces[0],
        &fixture.candidate,
        &interface,
        Point3::new(0.5, 0.7, 0.0),
        fixture.budget.pole_g0,
    )
    .expect("inner G0 prototype constructs");
    let measured = measure_inner_interface_g1(
        &fixture.sector,
        &fixture.strip,
        &inventory.interfaces[0],
        &fixture.candidate,
        &interface,
        &prototype,
        fixture.budget.interface_g1_deficit,
    )
    .expect("interface G1 passes");

    let report = validate_inner_pole_patch_prototype(
        &fixture.sector,
        &fixture.strip,
        &inventory.interfaces[0],
        &interface,
        &measured,
        fixture.budget,
    )
    .expect("measured inner prototype validates");

    assert!(report.interface_g1_validated);
    assert!(report.max_interface_g1_deficit <= fixture.budget.interface_g1_deficit);
    assert!(report.max_interface_g0 <= fixture.budget.interface_g0);
    assert!(report.max_pole_g0 <= fixture.budget.pole_g0);
}

struct StripFixture {
    sector: ExactBoundaryCapSectorInput,
    strip: crate::apex_construction::exact_boundary_cap::multipatch::inventory::ExactBoundaryOuterStripPatchInventory,
    boundary: NurbsCurve,
    samples: Vec<ExactBoundaryRibbonSample>,
    candidate: crate::apex_construction::exact_boundary_cap::multipatch::strip::ExactBoundaryOuterStripBoundaryScaffold,
    budget: ExactBoundaryMultiPatchBudget,
}

fn solved_fixture(support_normal: Vec3) -> StripFixture {
    let sector = sector_input();
    let inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FaceEdge,
    );
    let strip = inventory.outer_strips[0].clone();
    let boundary = cubic_boundary();
    let samples = support_samples(support_normal);
    let scaffold = scaffold_outer_strip_boundary(
        &sector,
        &strip,
        &boundary,
        ExactBoundaryStripBoundaryDirection::Forward,
        1.0e-10,
    )
    .expect("scaffold");
    let mut budget = ExactBoundaryMultiPatchBudget::default();
    budget.outer_boundary_residual = 1.0e-10;
    budget.support_normal_deficit = 1.0e-8;
    budget.min_normal_length = 1.0e-12;
    let candidate =
        solve_strip_support_normal_row(&sector, &strip, &scaffold, &boundary, &samples, budget)
            .unwrap_or(scaffold);
    StripFixture {
        sector,
        strip,
        boundary,
        samples,
        candidate,
        budget,
    }
}

fn cubic_boundary() -> NurbsCurve {
    NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.33, 0.1, 0.0),
            Point3::new(0.66, 0.1, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0; 4],
    }
}

fn support_samples(support_normal: Vec3) -> Vec<ExactBoundaryRibbonSample> {
    [0.0, 0.5, 1.0]
        .into_iter()
        .map(|parameter| ExactBoundaryRibbonSample {
            parameter,
            support_uv: (parameter, 0.0),
            boundary_point: Point3::new(parameter, 0.0, 0.0),
            support_point: Point3::new(parameter, 0.0, 0.0),
            boundary_tangent: Vec3::new(1.0, 0.0, 0.0),
            support_normal,
            cross_derivative: Vec3::new(0.0, 0.2, 0.0),
        })
        .collect()
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
