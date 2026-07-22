use super::nearest_branch_delta;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn nearest_branch_delta_normalizes_into_half_open_period() {
    let period = std::f64::consts::TAU;
    assert!((nearest_branch_delta(0.0, period)).abs() < 1e-12);
    assert!((nearest_branch_delta(period, period)).abs() < 1e-12);
    assert!((nearest_branch_delta(-period, period)).abs() < 1e-12);
    let d = nearest_branch_delta(period / 2.0 + 0.1, period);
    assert!((d - (0.1 - period / 2.0)).abs() < 1e-12);
    let d = nearest_branch_delta(-1.5 * period, period);
    assert!((d - period / 2.0).abs() < 1e-12);
}

#[test]
fn polyline_pcurve_collapses_straight_samples_to_two_points() {
    let samples: Vec<Point3> = (0..17)
        .map(|i| Point3::new(1.0 + 0.1 * f64::from(i), 2.0, 0.0))
        .collect();
    let pcurve = super::polyline_pcurve(&samples);
    let nurbs = pcurve.as_nurbs().expect("polyline pcurve is a NURBS");
    assert_eq!(nurbs.control_points.len(), 2);
    let (t0, t1) = pcurve.parameter_range();
    assert!((pcurve.point_at(t0).x - 1.0).abs() < 1e-12);
    assert!((pcurve.point_at(t1).x - 2.6).abs() < 1e-12);
}

#[test]
fn polyline_pcurve_keeps_curved_samples_as_polyline() {
    let samples: Vec<Point3> = (0..17)
        .map(|i| {
            let t = f64::from(i) / 16.0;
            Point3::new(t, t * t, 0.0)
        })
        .collect();
    let pcurve = super::polyline_pcurve(&samples);
    let nurbs = pcurve.as_nurbs().expect("polyline pcurve is a NURBS");
    assert_eq!(nurbs.control_points.len(), 17);
}

#[test]
fn reconcile_is_a_no_op_on_non_periodic_faces() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let faces: Vec<_> = cube.faces.keys().collect();
    for face in faces {
        assert!(!cube.reconcile_periodic_pcurve_branches(face));
    }
}
