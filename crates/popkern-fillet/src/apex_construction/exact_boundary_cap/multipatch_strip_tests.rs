use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::{ExactBoundaryStripBoundaryDirection, scaffold_outer_strip_boundary};
use crate::apex_construction::exact_boundary_cap::multipatch::inventory::ExactBoundaryMultiPatchInventory;
use crate::apex_construction::exact_boundary_cap::multipatch::report::{
    ExactBoundaryMultiPatchDeclineReason, ExactBoundaryMultiPatchValidationOrder,
};
use crate::apex_construction::exact_boundary_cap::sector::{
    ExactBoundaryCapSectorInput, ExactBoundaryOuterSeam, ExactBoundaryOuterSide,
};

#[test]
fn strip_scaffold_reuses_boundary_layout_and_rational_weights() {
    let sector = sector_input();
    let inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FaceEdge,
    );
    let strip = &inventory.outer_strips[0];
    let boundary = rational_quadratic_boundary();

    let scaffold = scaffold_outer_strip_boundary(
        &sector,
        strip,
        &boundary,
        ExactBoundaryStripBoundaryDirection::Forward,
        1.0e-12,
    )
    .expect("rational boundary should preserve weights");

    assert_eq!(scaffold.patch_name, strip.patch_name);
    assert_eq!(scaffold.surface.degree_u, boundary.degree);
    assert_eq!(scaffold.surface.knots_u, boundary.knots);
    assert_eq!(scaffold.surface.degree_v, 1);
    assert_eq!(scaffold.surface.num_u, boundary.control_points.len());
    assert_eq!(scaffold.surface.num_v, 2);
    for index in 0..boundary.control_points.len() {
        assert_eq!(
            scaffold.surface.control_points[index * 2],
            boundary.control_points[index]
        );
        assert_eq!(scaffold.surface.weights[index * 2], boundary.weights[index]);
    }
    assert!(scaffold.rational_boundary);
    assert!(scaffold.support_normal_row_pending);
    assert!(scaffold.max_outer_boundary_residual <= 1.0e-12);
}

#[test]
fn reversed_boundary_orientation_authors_reversed_outer_pcurve() {
    let sector = sector_input();
    let inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FilletEnd,
    );
    let strip = &inventory.outer_strips[0];
    let boundary = cubic_boundary();
    let (t_start, t_end) = boundary.parameter_range();

    let scaffold = scaffold_outer_strip_boundary(
        &sector,
        strip,
        &boundary,
        ExactBoundaryStripBoundaryDirection::Reversed,
        1.0e-12,
    )
    .expect("reversed parameterization should validate");

    assert_eq!(
        scaffold.boundary_direction,
        ExactBoundaryStripBoundaryDirection::Reversed
    );
    let start_uv = scaffold.outer_pcurve.point_at(t_start);
    let end_uv = scaffold.outer_pcurve.point_at(t_end);
    assert_eq!(start_uv.x, t_end);
    assert_eq!(end_uv.x, t_start);
    assert!(scaffold.max_outer_boundary_residual <= 1.0e-12);
}

#[test]
fn invalid_boundary_layout_declines_with_stable_diagnostic() {
    let sector = sector_input();
    let inventory = ExactBoundaryMultiPatchInventory::for_selected_side(
        &sector,
        ExactBoundaryOuterSide::FaceEdge,
    );
    let strip = &inventory.outer_strips[0];
    let mut boundary = cubic_boundary();
    boundary.weights.pop();

    let diagnostic = scaffold_outer_strip_boundary(
        &sector,
        strip,
        &boundary,
        ExactBoundaryStripBoundaryDirection::Forward,
        1.0e-12,
    )
    .expect_err("invalid boundary layout must decline");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryMultiPatchDeclineReason::PcurveParameterizationInvalid
    );
    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.multipatch.pcurve_parameterization_invalid"
    );
    assert_eq!(
        diagnostic.validation_order,
        ExactBoundaryMultiPatchValidationOrder::SourceParameterization
    );
    assert_eq!(
        diagnostic.patch_name.as_deref(),
        Some(strip.patch_name.as_str())
    );
    assert_eq!(diagnostic.edge, Some(strip.edge));
    assert_eq!(diagnostic.tolerance, Some(1.0e-12));
}

fn rational_quadratic_boundary() -> NurbsCurve {
    NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 0.25, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 0.7071067811865476, 1.0],
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
