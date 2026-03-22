// This file is only compiled under #[cfg(test)] (see mod.rs: `#[cfg(test)] mod tests;`).
#[cfg(test)]
use super::*;
use crate::geometry::{CurveEval, Point3};

#[test]
fn degree_1_nurbs_is_polyline() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let p = curve.point_at(0.5);
    assert!((p.x - 0.5).abs() < 1e-10);
    assert!((p.y).abs() < 1e-10);
}

#[test]
fn degree_1_endpoints() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(2.0, 3.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let p0 = curve.point_at(0.0);
    assert!((p0.x).abs() < 1e-10);
    let p1 = curve.point_at(1.0 - 1e-12);
    assert!((p1.x - 2.0).abs() < 1e-6);
}

#[test]
fn degree_2_quarter_circle() {
    let w = 1.0_f64 / 2.0_f64.sqrt();
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        weights: vec![1.0, w, 1.0],
    };
    // At t=0.5, point should be on the unit circle
    let p = curve.point_at(0.5);
    let radius = (p.x * p.x + p.y * p.y).sqrt();
    assert!(
        (radius - 1.0).abs() < 1e-10,
        "Expected radius 1.0, got {}",
        radius
    );
}

#[test]
fn degree_2_bezier_midpoint() {
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let p = curve.point_at(0.5);
    assert!((p.x - 0.5).abs() < 1e-10);
    assert!((p.y - 0.5).abs() < 1e-10);
}

#[test]
fn derivative_degree_1_is_constant() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(3.0, 4.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let d = curve.derivative_at(0.5, 1);
    assert!((d.x - 3.0).abs() < 1e-8);
    assert!((d.y - 4.0).abs() < 1e-8);
}

#[test]
fn length_degree_1() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(3.0, 4.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let len = curve.length_between(0.0, 1.0);
    assert!((len - 5.0).abs() < 1e-6, "Expected length 5.0, got {}", len);
}

#[test]
fn parameter_range() {
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let (t0, t1) = curve.parameter_range();
    assert_eq!(t0, 0.0);
    assert_eq!(t1, 1.0);
}

#[test]
fn tangent_at_degree_1() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(3.0, 4.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let t = curve.tangent_at(0.5);
    assert!((t.length() - 1.0).abs() < 1e-6, "Tangent should be unit");
}

#[test]
fn derivative_order_0() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let d0 = curve.derivative_at(0.5, 0);
    // order 0 returns position as vec
    assert!((d0.x - 0.5).abs() < 1e-6);
}

#[test]
fn derivative_order_2() {
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    // Second derivative should be non-zero for a quadratic
    let d2 = curve.derivative_at(0.5, 2);
    assert!(d2.length() > 1e-3);
}

#[test]
fn find_span_clamped_at_end() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    // Evaluate at exactly 1.0 (clamped)
    let p = curve.point_at(1.0);
    assert!((p.x - 1.0).abs() < 1e-10);
}

#[test]
fn find_span_clamped_before_start() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    // Evaluate before start — clamped to span p, so result is at t=0 origin
    let p = curve.point_at(-0.5);
    // The find_span clamps to p=1, so de Boor uses span 1 with t=-0.5
    // which extrapolates. Just check it doesn't panic.
    assert!(p.x.is_finite());
}

#[test]
fn multi_span_nurbs_curve() {
    // Degree 2 with 4 control points = 2 spans
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.25, 1.0, 0.0),
            Point3::new(0.75, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0, 1.0],
    };
    let p = curve.point_at(0.5);
    // Should be roughly at (0.5, something, 0)
    assert!((p.x - 0.5).abs() < 0.1);
    // Length should be computable
    let len = curve.length_between(0.0, 1.0);
    assert!(len > 1.0);
}

#[test]
fn curve_insert_knot_preserves_shape() {
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let new_curve = curve.insert_knot(0.5).unwrap();
    for &t in &[0.0, 0.25, 0.5, 0.75, 0.99] {
        let p1 = curve.point_at(t);
        let p2 = new_curve.point_at(t);
        assert!(
            (p1.x - p2.x).abs() < 1e-10,
            "x mismatch at t={t}: {p1:?} vs {p2:?}"
        );
        assert!(
            (p1.y - p2.y).abs() < 1e-10,
            "y mismatch at t={t}: {p1:?} vs {p2:?}"
        );
    }
}

#[test]
fn curve_insert_knot_increases_counts() {
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let new_curve = curve.insert_knot(0.3).unwrap();
    assert_eq!(new_curve.knots.len(), curve.knots.len() + 1);
    assert_eq!(
        new_curve.control_points.len(),
        curve.control_points.len() + 1
    );
    assert_eq!(new_curve.weights.len(), curve.weights.len() + 1);
}

#[test]
fn curve_insert_knot_out_of_range() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    assert!(curve.insert_knot(2.0).is_err());
    assert!(curve.insert_knot(-0.5).is_err());
}

#[test]
fn curve_insert_knots_preserves_shape() {
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let new_curve = curve.insert_knots(&[0.25, 0.5, 0.75]).unwrap();
    for &t in &[0.0, 0.25, 0.5, 0.75, 0.99] {
        let p1 = curve.point_at(t);
        let p2 = new_curve.point_at(t);
        assert!(
            (p1.x - p2.x).abs() < 1e-10,
            "x mismatch at t={t}: {p1:?} vs {p2:?}"
        );
        assert!(
            (p1.y - p2.y).abs() < 1e-10,
            "y mismatch at t={t}: {p1:?} vs {p2:?}"
        );
    }
}

#[test]
fn curve_insert_knots_out_of_range() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    assert!(curve.insert_knots(&[0.5, 2.0]).is_err());
}

#[test]
fn curve_remove_knot_roundtrip() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let inserted = curve.insert_knot(0.5).unwrap();
    assert_eq!(inserted.knots.len(), 5);
    let removed = inserted.remove_knot(2, 0.0);
    assert!(removed.is_some());
    assert_eq!(removed.unwrap().knots.len(), 4);
}

#[test]
fn curve_remove_knot_times() {
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let inserted = curve.insert_knots(&[0.5, 0.5]).unwrap();
    // After inserting 0.5 twice: knots = [0,0,0, 0.5,0.5, 1,1,1], index 3 holds first 0.5
    let (result, count) = inserted.remove_knot_times(3, 3, 0.0);
    assert_eq!(count, 2);
    // Two knots removed from the inserted curve
    assert_eq!(result.knots.len(), inserted.knots.len() - 2);
}

#[test]
fn curve_remove_knot_none_on_boundary() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    // Boundary knots cannot be removed
    assert!(curve.remove_knot(0, 0.0).is_none());
    assert!(curve.remove_knot(3, 0.0).is_none());
}

#[test]
fn curve_elevate_preserves_shape() {
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let elevated = curve.elevate_degree();
    assert_eq!(elevated.degree, 3);
    for t in [0.0, 0.25, 0.5, 0.75, 0.99] {
        let p1 = curve.point_at(t);
        let p2 = elevated.point_at(t);
        assert!(
            (p1.x - p2.x).abs() < 1e-10,
            "x mismatch at t={t}: {p1:?} vs {p2:?}"
        );
        assert!(
            (p1.y - p2.y).abs() < 1e-10,
            "y mismatch at t={t}: {p1:?} vs {p2:?}"
        );
    }
}

#[test]
fn curve_elevate_to_target() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let elevated = curve.elevate_degree_to(4).unwrap();
    assert_eq!(elevated.degree, 4);
}

#[test]
fn curve_elevate_to_same_fails() {
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    assert!(curve.elevate_degree_to(2).is_err());
    assert!(curve.elevate_degree_to(1).is_err());
}

#[test]
fn curve_split_preserves_shape() {
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let (left, right) = curve.split_at(0.5).unwrap();

    for t in [0.0, 0.1, 0.2, 0.3, 0.4, 0.49] {
        let po = curve.point_at(t);
        let pl = left.point_at(t);
        assert!(
            (po.x - pl.x).abs() < 1e-8,
            "left x mismatch at t={t}: {po:?} vs {pl:?}"
        );
        assert!(
            (po.y - pl.y).abs() < 1e-8,
            "left y mismatch at t={t}: {po:?} vs {pl:?}"
        );
    }
    for t in [0.51, 0.6, 0.7, 0.8, 0.9, 0.99] {
        let po = curve.point_at(t);
        let pr = right.point_at(t);
        assert!(
            (po.x - pr.x).abs() < 1e-8,
            "right x mismatch at t={t}: {po:?} vs {pr:?}"
        );
        assert!(
            (po.y - pr.y).abs() < 1e-8,
            "right y mismatch at t={t}: {po:?} vs {pr:?}"
        );
    }
}

#[test]
fn curve_split_at_endpoint_fails() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    assert!(curve.split_at(0.0).is_err());
    assert!(curve.split_at(1.0).is_err());
}

#[test]
fn derivative_at_boundary_zero_dt_guard() {
    // Very short parameter range so that t+h = t-h after clamping => dt ≈ 0
    // This shouldn't happen with normal curves, but the guard exists
    // We can trigger it indirectly with order 2 at the boundary
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1e-40, 1e-40], // extremely tiny range
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    // derivative_at with order 2 at t=0 should hit the dt guard for order 2
    let d = curve.derivative_at(0.0, 2);
    assert!(d.x.is_finite());
}

#[test]
fn derivative_order_1_at_start_boundary() {
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    // At t=0, t_bwd is clamped to t_min=0, so forward difference only
    let d = curve.derivative_at(0.0, 1);
    assert!(d.length() > 0.1); // should be a valid non-zero tangent
}
