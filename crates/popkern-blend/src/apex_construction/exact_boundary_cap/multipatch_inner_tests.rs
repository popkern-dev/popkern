use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::{
    construct_inner_pole_patch_prototype, construct_strip_inner_interface_curve,
    decline_g0_only_inner_prototype_before_commit, measure_inner_interface_g1,
};
use crate::apex_construction::exact_boundary_cap::multipatch::budget::ExactBoundaryMultiPatchBudget;
use crate::apex_construction::exact_boundary_cap::multipatch::inventory::ExactBoundaryMultiPatchInventory;
use crate::apex_construction::exact_boundary_cap::multipatch::report::{
    ExactBoundaryMultiPatchDeclineReason, ExactBoundaryMultiPatchValidationOrder,
};
use crate::apex_construction::exact_boundary_cap::multipatch::strip::{
    ExactBoundaryOuterStripBoundaryScaffold, ExactBoundaryStripBoundaryDirection,
    scaffold_outer_strip_boundary, solve_strip_support_normal_row,
};
use crate::apex_construction::exact_boundary_cap::ribbon::ExactBoundaryRibbonSample;
use crate::apex_construction::exact_boundary_cap::sector::{
    ExactBoundaryCapSectorInput, ExactBoundaryOuterSeam, ExactBoundaryOuterSide,
};

#[test]
fn interface_curve_matches_strip_inner_boundary_with_owned_names() {
    let fixture = strip_fixture(ExactBoundaryStripBoundaryDirection::Forward);
    let interface = construct_strip_inner_interface_curve(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &fixture.inventory.interfaces[0],
        &fixture.strip_candidate,
        1.0e-6,
    )
    .expect("interface curve constructs");

    assert_eq!(
        interface.interface_name,
        "direct_apex.cap_sector.0.face_edge_interface"
    );
    assert_eq!(
        interface.strip_patch_name,
        fixture.inventory.outer_strips[0].patch_name
    );
    assert_eq!(
        interface.inner_patch_name,
        fixture
            .inventory
            .inner_patch
            .as_ref()
            .expect("inner patch")
            .patch_name
    );
    assert_eq!(
        interface.curve.degree,
        fixture.strip_candidate.surface.degree_u
    );
    assert_eq!(
        interface.curve.knots,
        fixture.strip_candidate.surface.knots_u
    );
    assert!(interface.max_strip_boundary_residual <= 1.0e-6);
    assert_eq!(interface.sample_count, 17);
}

#[test]
fn interface_pcurve_is_finite_and_monotone_for_forward_and_reversed_strips() {
    let forward =
        interface_for_fixture(&strip_fixture(ExactBoundaryStripBoundaryDirection::Forward));
    let reversed = interface_for_fixture(&strip_fixture(
        ExactBoundaryStripBoundaryDirection::Reversed,
    ));

    let f0 = forward.strip_pcurve.point_at(0.0).x;
    let f1 = forward.strip_pcurve.point_at(1.0).x;
    let r0 = reversed.strip_pcurve.point_at(0.0).x;
    let r1 = reversed.strip_pcurve.point_at(1.0).x;

    assert!(f0.is_finite() && f1.is_finite() && f0 < f1);
    assert!(r0.is_finite() && r1.is_finite() && r0 > r1);
    assert!(!forward.strip_parameter_reversed);
    assert!(reversed.strip_parameter_reversed);
}

#[test]
fn interface_curve_evaluates_to_strip_inner_boundary_samples() {
    let fixture = strip_fixture(ExactBoundaryStripBoundaryDirection::Forward);
    let interface = construct_strip_inner_interface_curve(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &fixture.inventory.interfaces[0],
        &fixture.strip_candidate,
        1.0e-6,
    )
    .expect("interface curve constructs");

    for parameter in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let uv = interface.strip_pcurve.point_at(parameter);
        let strip_point = fixture.strip_candidate.surface.point_at(uv.x, uv.y);
        let curve_point = interface.curve.point_at(uv.x);
        assert!(strip_point.distance_to(&curve_point) <= 1.0e-6);
    }
}

#[test]
fn interface_ownership_mismatch_declines_with_interface_diagnostic() {
    let fixture = strip_fixture(ExactBoundaryStripBoundaryDirection::Forward);
    let mut mismatched = fixture.inventory.interfaces[0].clone();
    mismatched.strip_patch_name = "direct_apex.cap_sector.0.other_strip".to_string();

    let diagnostic = construct_strip_inner_interface_curve(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &mismatched,
        &fixture.strip_candidate,
        1.0e-6,
    )
    .expect_err("mismatched interface ownership declines");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::InterfacePositionUnsatisfied
    );
    assert_eq!(
        diagnostic.validation_order,
        ExactBoundaryMultiPatchValidationOrder::InterfacePosition
    );
    assert_eq!(
        diagnostic.interface_name.as_deref(),
        Some(mismatched.interface_name.as_str())
    );
    assert_eq!(
        diagnostic.patch_name.as_deref(),
        Some(fixture.inventory.outer_strips[0].patch_name.as_str())
    );
}

#[test]
fn inner_pole_patch_prototype_connects_interface_curve_and_pole_g0() {
    let fixture = strip_fixture(ExactBoundaryStripBoundaryDirection::Forward);
    let interface = interface_for_fixture(&fixture);
    let pole = Point3::new(0.5, 0.7, 0.0);

    let prototype = construct_inner_pole_patch_prototype(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &fixture.inventory.interfaces[0],
        &fixture.strip_candidate,
        &interface,
        pole,
        1.0e-6,
    )
    .expect("inner G0 prototype constructs");

    assert_eq!(prototype.patch_name, interface.inner_patch_name);
    assert_eq!(prototype.interface_name, interface.interface_name);
    assert_eq!(prototype.surface.degree_u, interface.curve.degree);
    assert_eq!(prototype.surface.knots_u, interface.curve.knots);
    assert_eq!(prototype.surface.degree_v, 2);
    assert!(prototype.max_interface_g0 <= 1.0e-6);
    assert!(prototype.max_pole_g0 <= 1.0e-6);
    assert!(!prototype.interface_g1_validated);
    for parameter in [0.0, 0.5, 1.0] {
        let uv = prototype.interface_pcurve.point_at(parameter);
        let patch_point = prototype.surface.point_at(uv.x, uv.y);
        let interface_point = interface.curve.point_at(uv.x);
        assert!(patch_point.distance_to(&interface_point) <= 1.0e-6);
        assert!(
            prototype
                .surface
                .point_at(parameter, 1.0)
                .distance_to(&pole)
                <= 1.0e-6
        );
    }
}

#[test]
fn g0_only_inner_prototype_declines_before_commit_with_g1_diagnostic() {
    let fixture = strip_fixture(ExactBoundaryStripBoundaryDirection::Forward);
    let interface = interface_for_fixture(&fixture);
    let prototype = construct_inner_pole_patch_prototype(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &fixture.inventory.interfaces[0],
        &fixture.strip_candidate,
        &interface,
        Point3::new(0.5, 0.7, 0.0),
        1.0e-6,
    )
    .expect("inner G0 prototype constructs");

    let diagnostic = decline_g0_only_inner_prototype_before_commit(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &fixture.inventory.interfaces[0],
        &interface,
        &prototype,
        5.0e-2,
    );

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::InterfaceTangentUnsatisfied
    );
    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.multipatch.interface_g1_exceeded"
    );
    assert_eq!(
        diagnostic.validation_order,
        ExactBoundaryMultiPatchValidationOrder::InterfaceTangent
    );
    assert_eq!(
        diagnostic.interface_name.as_deref(),
        Some(fixture.inventory.interfaces[0].interface_name.as_str())
    );
    assert_eq!(
        diagnostic.sample_parameter,
        Some(interface.sample_parameter)
    );
    assert_eq!(diagnostic.tolerance, Some(5.0e-2));
}

#[test]
fn inner_pole_patch_reports_interface_g0_and_pole_g0_diagnostics() {
    let fixture = strip_fixture(ExactBoundaryStripBoundaryDirection::Forward);
    let mut interface = interface_for_fixture(&fixture);
    interface.curve.control_points[0].z = f64::NAN;
    let interface_diag = construct_inner_pole_patch_prototype(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &fixture.inventory.interfaces[0],
        &fixture.strip_candidate,
        &interface,
        Point3::new(0.5, 0.7, 0.0),
        1.0e-6,
    )
    .expect_err("interface G0 mismatch declines");

    assert_eq!(
        interface_diag.reason,
        ExactBoundaryMultiPatchDeclineReason::InterfacePositionUnsatisfied
    );
    assert_eq!(
        interface_diag.code,
        "exact_boundary_cap.multipatch.interface_g0_exceeded"
    );

    let interface = interface_for_fixture(&fixture);
    let pole_diag = construct_inner_pole_patch_prototype(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &fixture.inventory.interfaces[0],
        &fixture.strip_candidate,
        &interface,
        Point3::new(f64::NAN, 0.7, 0.0),
        1.0e-6,
    )
    .expect_err("invalid pole declines");

    assert_eq!(
        pole_diag.reason,
        ExactBoundaryMultiPatchDeclineReason::InnerPolePositionUnsatisfied
    );
    assert_eq!(
        pole_diag.code,
        "exact_boundary_cap.multipatch.pole_g0_exceeded"
    );
}

#[test]
fn interface_g1_measurement_marks_compatible_candidate_validated() {
    let fixture = strip_fixture(ExactBoundaryStripBoundaryDirection::Forward);
    let interface = interface_for_fixture(&fixture);
    let prototype = construct_inner_pole_patch_prototype(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &fixture.inventory.interfaces[0],
        &fixture.strip_candidate,
        &interface,
        Point3::new(0.5, 0.7, 0.0),
        1.0e-6,
    )
    .expect("inner G0 prototype constructs");

    let measured = measure_inner_interface_g1(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &fixture.inventory.interfaces[0],
        &fixture.strip_candidate,
        &interface,
        &prototype,
        5.0e-2,
    )
    .expect("coplanar prototype passes G1");
    let report = measured.interface_g1_report.expect("G1 report");

    assert!(measured.interface_g1_validated);
    assert!(report.max_interface_g1_deficit <= 5.0e-2);
    assert_eq!(report.sample_count, 17);
    assert_eq!(report.tolerance, 5.0e-2);
}

#[test]
fn interface_g1_measurement_is_decoupled_from_pole_row_height() {
    let fixture = strip_fixture(ExactBoundaryStripBoundaryDirection::Forward);
    let interface = interface_for_fixture(&fixture);
    let prototype = construct_inner_pole_patch_prototype(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &fixture.inventory.interfaces[0],
        &fixture.strip_candidate,
        &interface,
        Point3::new(0.5, 0.7, 0.8),
        1.0e-6,
    )
    .expect("tilted inner G0 prototype constructs");

    let measured = measure_inner_interface_g1(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &fixture.inventory.interfaces[0],
        &fixture.strip_candidate,
        &interface,
        &prototype,
        5.0e-2,
    )
    .expect("interface G1 row is solved independently from pole height");

    assert!(measured.interface_g1_validated);
    assert!(
        measured
            .interface_g1_report
            .is_some_and(|report| report.max_interface_g1_deficit <= 5.0e-2)
    );
}

#[test]
fn interface_g1_measurement_leaves_endpoints_to_pole_closure() {
    let fixture = strip_fixture(ExactBoundaryStripBoundaryDirection::Forward);
    let interface = interface_for_fixture(&fixture);
    let prototype = construct_inner_pole_patch_prototype(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &fixture.inventory.interfaces[0],
        &fixture.strip_candidate,
        &interface,
        Point3::new(0.5, 0.7, 1.5),
        1.0e-6,
    )
    .expect("tilted pole prototype constructs");

    assert!(
        measure_inner_interface_g1(
            &fixture.sector,
            &fixture.inventory.outer_strips[0],
            &fixture.inventory.interfaces[0],
            &fixture.strip_candidate,
            &interface,
            &prototype,
            5.0e-2,
        )
        .expect("endpoint-only pole closure should not block open interface G1")
        .interface_g1_report
        .is_some_and(|report| report.max_interface_g1_deficit <= 5.0e-2)
    );
}

fn interface_for_fixture(fixture: &StripFixture) -> super::ExactBoundaryStripInnerInterfaceCurve {
    construct_strip_inner_interface_curve(
        &fixture.sector,
        &fixture.inventory.outer_strips[0],
        &fixture.inventory.interfaces[0],
        &fixture.strip_candidate,
        1.0e-6,
    )
    .expect("interface curve constructs")
}

struct StripFixture {
    sector: ExactBoundaryCapSectorInput,
    inventory: ExactBoundaryMultiPatchInventory,
    strip_candidate: ExactBoundaryOuterStripBoundaryScaffold,
}

fn strip_fixture(direction: ExactBoundaryStripBoundaryDirection) -> StripFixture {
    let sector = sector_input();
    let inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FaceEdge,
    );
    let boundary = cubic_boundary();
    let scaffold = scaffold_outer_strip_boundary(
        &sector,
        &inventory.outer_strips[0],
        &boundary,
        direction,
        1.0e-10,
    )
    .expect("scaffold");
    let budget = ExactBoundaryMultiPatchBudget {
        support_normal_deficit: 1.0e-8,
        ..Default::default()
    };
    let strip_candidate = solve_strip_support_normal_row(
        &sector,
        &inventory.outer_strips[0],
        &scaffold,
        &boundary,
        &support_samples(),
        budget,
    )
    .expect("strip support-normal row");
    StripFixture {
        sector,
        inventory,
        strip_candidate,
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

fn support_samples() -> Vec<ExactBoundaryRibbonSample> {
    [0.0, 0.5, 1.0]
        .into_iter()
        .map(|parameter| ExactBoundaryRibbonSample {
            parameter,
            support_uv: (parameter, 0.0),
            boundary_point: Point3::new(parameter, 0.0, 0.0),
            support_point: Point3::new(parameter, 0.0, 0.0),
            boundary_tangent: Vec3::new(1.0, 0.0, 0.0),
            support_normal: Vec3::new(0.0, 0.0, 1.0),
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
