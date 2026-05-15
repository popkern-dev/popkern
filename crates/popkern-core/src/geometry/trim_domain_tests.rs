use super::*;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::{Point3, SurfaceGeometry};

fn unit_square_surface() -> SurfaceGeometry {
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    })
}

#[test]
fn classifies_boundary_subspan_without_requiring_full_edge() {
    let surface = unit_square_surface();
    let pcurve = Pcurve::linear((0.0, 0.2), (0.0, 0.8), 0.0, 1.0);
    let anchors = (Point3::new(0.0, 0.2, 0.0), Point3::new(0.0, 0.8, 0.0));

    let candidate = classify_trim_domain_candidate(&surface, &pcurve, Some(anchors), 1.0e-9, 9)
        .expect("boundary subspan should classify");

    assert_eq!(
        candidate.kind,
        TrimDomainCandidateKind::BoundarySubspan {
            side: BoundarySide::UMin
        }
    );
    assert_eq!(candidate.endpoint_gap, 0.0);
}

#[test]
fn classifies_valid_interior_trim_separately_from_boundary_ownership() {
    let surface = unit_square_surface();
    let pcurve = Pcurve::linear((0.25, 0.2), (0.75, 0.8), 0.0, 1.0);
    let anchors = (Point3::new(0.25, 0.2, 0.0), Point3::new(0.75, 0.8, 0.0));

    let candidate = classify_trim_domain_candidate(&surface, &pcurve, Some(anchors), 1.0e-9, 9)
        .expect("interior trim should classify");

    assert_eq!(candidate.kind, TrimDomainCandidateKind::Interior);
    assert!(candidate.best_boundary_deviation > 0.0);
}

#[test]
fn declines_inside_domain_pcurve_with_wrong_topology_endpoints() {
    let surface = unit_square_surface();
    let pcurve = Pcurve::linear((0.25, 0.2), (0.75, 0.8), 0.0, 1.0);
    let anchors = (Point3::new(0.25, 0.2, 0.0), Point3::new(0.75, 0.9, 0.0));

    let decline = classify_trim_domain_candidate(&surface, &pcurve, Some(anchors), 1.0e-9, 9)
        .expect_err("endpoint mismatch must decline");

    assert_eq!(decline.reason, TrimDomainDeclineReason::EndpointMismatch);
    assert!(decline.endpoint_gap > 0.09);
}
