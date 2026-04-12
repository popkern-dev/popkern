use super::*;
use crate::geometry::Point3;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::healing::repair::hole_fill::DetectedHole;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn fill_n_sided_hole_rejects_fewer_than_3_sides() {
    let brep = BRep::new(Tolerance::default());
    let hole = DetectedHole {
        boundary_edges: vec![],
        num_sides: 2,
    };
    assert!(fill_n_sided_hole(&brep, &hole, 1e-6).is_none());
}

#[test]
fn centroid_of_triangle() {
    let pts = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(3.0, 0.0, 0.0),
        Point3::new(0.0, 3.0, 0.0),
    ];
    let c = centroid(&pts);
    assert!(c.is_some());
    let c = c.as_ref();
    assert!((c.map(|c| c.x).unwrap_or(0.0) - 1.0).abs() < 1e-10);
    assert!((c.map(|c| c.y).unwrap_or(0.0) - 1.0).abs() < 1e-10);
    assert!((c.map(|c| c.z).unwrap_or(0.0)).abs() < 1e-10);
}

#[test]
fn centroid_of_empty_returns_none() {
    assert!(centroid(&[]).is_none());
}

#[test]
fn line_nurbs_endpoints_are_correct() {
    let start = Point3::new(1.0, 2.0, 3.0);
    let end = Point3::new(4.0, 5.0, 6.0);
    let curve = line_nurbs(start, end);

    use crate::geometry::curve::CurveEval;
    let p0 = curve.point_at(0.0);
    let p1 = curve.point_at(1.0);

    assert!(p0.distance_to(&start) < 1e-10, "start mismatch");
    assert!(p1.distance_to(&end) < 1e-10, "end mismatch");
}

#[test]
fn line_nurbs_is_degree_1() {
    let curve = line_nurbs(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
    assert_eq!(curve.degree, 1);
    assert_eq!(curve.control_points.len(), 2);
    assert_eq!(curve.weights.len(), 2);
    assert_eq!(curve.knots.len(), 4);
}
