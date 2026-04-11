//! Newton-Raphson refinement for intersection parameter convergence.

use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::tolerance::Tolerance;

/// Refine an approximate intersection using Newton-Raphson.
/// Returns Some((t1, t2, point)) if converged, None otherwise.
#[allow(dead_code)]
pub(crate) fn newton_refine(
    curve1: &CurveGeometry,
    t1_approx: f64,
    curve2: &CurveGeometry,
    t2_approx: f64,
    tolerance: &Tolerance,
) -> Option<(f64, f64, Point3)> {
    let max_iter = 20;
    let mut t1 = t1_approx;
    let mut t2 = t2_approx;

    for _ in 0..max_iter {
        let p1 = curve1.point_at(t1);
        let p2 = curve2.point_at(t2);
        let f = p1 - p2; // Vec3 residual

        if f.length() < tolerance.linear() {
            return Some((t1, t2, p1));
        }

        let d1 = curve1.derivative_at(t1, 1);
        let d2 = curve2.derivative_at(t2, 1);

        // Solve projected 2x2 system
        let a = d1.dot(&d1);
        let b = -d1.dot(&d2);
        let c_val = -d2.dot(&d1);
        let d = d2.dot(&d2);
        let det = a * d - b * c_val;

        if det.abs() < 1e-30 {
            return None;
        }

        let rhs1 = d1.dot(&f);
        let rhs2 = -d2.dot(&f);

        let dt1 = (d * rhs1 - b * rhs2) / det;
        let dt2 = (a * rhs2 - c_val * rhs1) / det;

        t1 -= dt1;
        t2 -= dt2;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::CurveGeometry;
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;

    #[test]
    fn refines_line_line() {
        let c1 = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        let c2 = CurveGeometry::Line {
            origin: Point3::new(0.5, -1.0, 0.0),
            direction: Vec3::new(0.0, 1.0, 0.0),
        };
        let r = newton_refine(&c1, 0.4, &c2, 0.9, &Tolerance::default());
        assert!(r.is_some());
        let (t1, t2, _) = r.unwrap();
        assert!((t1 - 0.5).abs() < 1e-12);
        assert!((t2 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn fails_for_parallel_lines() {
        let c1 = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        let c2 = CurveGeometry::Line {
            origin: Point3::new(0.0, 10.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        let r = newton_refine(&c1, 0.5, &c2, 0.5, &Tolerance::default());
        assert!(r.is_none());
    }

    #[test]
    fn refines_with_nurbs_curves() {
        use crate::geometry::nurbs_curve::NurbsCurve;
        // Quadratic Bezier: x(t)=t, y(t)=2t(1-t).
        // At y=0.375 the curve crosses transversally at t≈0.25 and t≈0.75.
        let c1 = CurveGeometry::Nurbs(NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        });
        let c2 = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.375, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        // Start near the first crossing (t1≈0.25, t2≈0.25).
        let r = newton_refine(&c1, 0.20, &c2, 0.20, &Tolerance::default());
        assert!(r.is_some());
        let (t1, t2, point) = r.unwrap();
        let _ = (t1, t2);
        // The intersection point should be at y ≈ 0.375
        assert!((point.y - 0.375).abs() < 1e-8);
    }
}
