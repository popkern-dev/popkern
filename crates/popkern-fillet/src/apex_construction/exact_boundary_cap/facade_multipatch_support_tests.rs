use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::*;
use crate::apex_construction::exact_boundary_cap::boundary::PreparedExactBoundaryCurve;
use crate::apex_construction::exact_boundary_cap::multipatch::report::ExactBoundaryMultiPatchDeclineReason;
use crate::apex_construction::exact_boundary_cap::pcurves::PreparedExactBoundaryPcurvePair;
use crate::apex_construction::exact_boundary_cap::sector::ExactBoundaryOuterSide;

#[test]
fn support_authority_boundary_reauthors_from_support_pcurve_image() {
    let boundary = prepared_boundary(curve_with_midpoint_lift(0.25, 1.0));
    let support_pair = support_pair(Point3::new(1.0, 0.0, 0.0));

    let candidate = support_authority_boundary_curve(&boundary, &support_pair)
        .expect("support pcurve image should reauthor coherent boundary");

    assert!(candidate.point_at(0.5).z.abs() <= 1.0e-9);
    assert!(tangent_support_normal_dot(&candidate, &support_pair) <= 1.0e-9);
}

#[test]
fn support_authority_boundary_declines_when_endpoints_change_topology() {
    let boundary = prepared_boundary(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.1, 0.0), Point3::new(1.0, 0.1, 0.0)],
        weights: vec![1.0, 1.0],
    });
    let support_pair = support_pair(Point3::new(1.0, 0.0, 0.0));

    let attempt = support_authority_boundary_attempt(&boundary, &support_pair);

    assert!(attempt.candidate.is_none());
    let decline = attempt.decline.expect("endpoint mismatch decline");
    assert_eq!(
        decline.reason,
        ExactBoundaryMultiPatchDeclineReason::SupportBoundaryReauthorEndpointMismatch
    );
    assert!(decline.residual > decline.tolerance);
}

fn prepared_boundary(curve: NurbsCurve) -> PreparedExactBoundaryCurve {
    PreparedExactBoundaryCurve {
        sector_index: 0,
        sector_name: "direct_apex.cap_sector.0".to_string(),
        side: ExactBoundaryOuterSide::FaceEdge,
        edge: EdgeHandle::default(),
        curve,
    }
}

fn support_pair(end: Point3) -> PreparedExactBoundaryPcurvePair {
    PreparedExactBoundaryPcurvePair {
        sector_index: 0,
        sector_name: "direct_apex.cap_sector.0".to_string(),
        side: ExactBoundaryOuterSide::FaceEdge,
        edge: EdgeHandle::default(),
        support_half_edge: HalfEdgeHandle::default(),
        support_face: FaceHandle::default(),
        support_surface: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        cap_pcurve: Pcurve::linear((0.0, 0.0), (end.x, end.y), 0.0, 1.0),
        support_pcurve: Pcurve::linear((0.0, 0.0), (end.x, end.y), 0.0, 1.0),
        max_cap_residual: 0.0,
        max_support_residual: 0.0,
        max_pair_residual: 0.0,
    }
}

fn curve_with_midpoint_lift(lift: f64, end_x: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 0.0, lift),
            Point3::new(end_x, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    }
}
