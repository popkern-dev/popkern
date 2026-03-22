//! NURBS splitting via knot insertion to full multiplicity.

use super::knot_insertion::insert_knot_1d;
use crate::error::CoreError;
use crate::geometry::Point3;

/// Split a 1D NURBS at parameter `t`.
/// Returns (left_knots, left_pts, left_wts, right_knots, right_pts, right_wts).
///
/// The algorithm inserts the knot `t` until its multiplicity equals `degree + 1`,
/// which causes the curve to split naturally at that parameter. The resulting
/// control points and knots are then partitioned into two independent curves
/// that together reproduce the original curve exactly.
#[allow(dead_code, clippy::type_complexity)]
pub(crate) fn split_1d(
    degree: usize,
    knots: &[f64],
    control_points: &[Point3],
    weights: &[f64],
    t: f64,
) -> Result<
    (
        Vec<f64>,
        Vec<Point3>,
        Vec<f64>,
        Vec<f64>,
        Vec<Point3>,
        Vec<f64>,
    ),
    CoreError,
> {
    let t_min = knots[degree];
    let t_max = knots[knots.len() - degree - 1];
    if t <= t_min || t >= t_max {
        return Err(CoreError::PreconditionFailed(format!(
            "split parameter {} must be strictly inside ({}, {})",
            t, t_min, t_max
        )));
    }

    // Count current multiplicity of t
    let current_mult = knots.iter().filter(|&&k| (k - t).abs() < 1e-15).count();
    let insertions_needed = degree + 1 - current_mult;

    // Insert t until multiplicity = degree + 1
    let mut k = knots.to_vec();
    let mut pts = control_points.to_vec();
    let mut wts = weights.to_vec();
    for _ in 0..insertions_needed {
        let (nk, np, nw) = insert_knot_1d(degree, &k, &pts, &wts, t);
        k = nk;
        pts = np;
        wts = nw;
    }

    // Find the first index where t appears in the knot vector.
    let first_t_idx = k.iter().position(|&x| (x - t).abs() < 1e-15).unwrap_or(0);

    // After inserting t to multiplicity (degree+1), the knot vector has a block
    // of (degree+1) copies of t starting at first_t_idx. This block acts as the
    // boundary between left and right curves:
    //
    //   Left knots:  k[0 .. first_t_idx + degree + 1]  (ends with degree+1 copies of t)
    //   Right knots: k[first_t_idx ..]                  (starts with degree+1 copies of t)
    //   Left CPs:    pts[0 .. first_t_idx]
    //   Right CPs:   pts[first_t_idx ..]
    //
    // Both halves are clamped at t: the last CP of the left and the first CP of
    // the right both evaluate to C(t).

    let left_pts = pts[..first_t_idx].to_vec();
    let left_wts = wts[..first_t_idx].to_vec();
    let left_knots: Vec<f64> = k[..first_t_idx + degree + 1].to_vec();

    let right_pts = pts[first_t_idx..].to_vec();
    let right_wts = wts[first_t_idx..].to_vec();
    let right_knots: Vec<f64> = k[first_t_idx..].to_vec();

    // Defensive checks
    debug_assert_eq!(
        left_knots.len(),
        left_pts.len() + degree + 1,
        "Left knot vector length mismatch"
    );
    debug_assert_eq!(
        right_knots.len(),
        right_pts.len() + degree + 1,
        "Right knot vector length mismatch"
    );

    Ok((
        left_knots,
        left_pts,
        left_wts,
        right_knots,
        right_pts,
        right_wts,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;
    use crate::geometry::curve::CurveEval;
    use crate::geometry::nurbs_curve::NurbsCurve;

    #[test]
    fn split_degree_1_line() {
        let pts = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        let wts = vec![1.0, 1.0];
        let knots = vec![0.0, 0.0, 1.0, 1.0];
        let result = split_1d(1, &knots, &pts, &wts, 0.5);
        assert!(result.is_ok());
        let (lk, lp, lw, rk, rp, rw) = result.unwrap();

        // Verify knot vector lengths
        assert_eq!(lk.len(), lp.len() + 2, "Left knots wrong length");
        assert_eq!(rk.len(), rp.len() + 2, "Right knots wrong length");

        // Verify shape
        let original = NurbsCurve {
            degree: 1,
            knots,
            control_points: pts,
            weights: wts,
        };
        let left = NurbsCurve {
            degree: 1,
            knots: lk,
            control_points: lp,
            weights: lw,
        };
        let right = NurbsCurve {
            degree: 1,
            knots: rk,
            control_points: rp,
            weights: rw,
        };

        // Left half: [0, 0.5]
        for t in [0.0, 0.1, 0.2, 0.3, 0.4] {
            let po = original.point_at(t);
            let pl = left.point_at(t);
            assert!((po.x - pl.x).abs() < 1e-10, "Left differs at t={t}");
        }

        // Right half: [0.5, 1.0]
        for t in [0.6, 0.7, 0.8, 0.9, 0.99] {
            let po = original.point_at(t);
            let pr = right.point_at(t);
            assert!((po.x - pr.x).abs() < 1e-10, "Right differs at t={t}");
        }
    }

    #[test]
    fn split_degree_2_bezier() {
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0];
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let (lk, lp, lw, rk, rp, rw) = split_1d(2, &knots, &pts, &wts, 0.5).unwrap();

        assert_eq!(lk.len(), lp.len() + 3);
        assert_eq!(rk.len(), rp.len() + 3);

        let original = NurbsCurve {
            degree: 2,
            knots,
            control_points: pts,
            weights: wts,
        };
        let left = NurbsCurve {
            degree: 2,
            knots: lk,
            control_points: lp,
            weights: lw,
        };
        let right = NurbsCurve {
            degree: 2,
            knots: rk,
            control_points: rp,
            weights: rw,
        };

        for t in [0.0, 0.1, 0.2, 0.3, 0.4, 0.49] {
            let po = original.point_at(t);
            let pl = left.point_at(t);
            assert!((po.x - pl.x).abs() < 1e-8, "Left x at t={t}");
            assert!((po.y - pl.y).abs() < 1e-8, "Left y at t={t}");
        }
        for t in [0.51, 0.6, 0.7, 0.8, 0.9, 0.99] {
            let po = original.point_at(t);
            let pr = right.point_at(t);
            assert!((po.x - pr.x).abs() < 1e-8, "Right x at t={t}");
            assert!((po.y - pr.y).abs() < 1e-8, "Right y at t={t}");
        }
    }

    #[test]
    fn split_at_endpoint_fails() {
        let pts = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        let wts = vec![1.0, 1.0];
        let knots = vec![0.0, 0.0, 1.0, 1.0];
        assert!(split_1d(1, &knots, &pts, &wts, 0.0).is_err());
        assert!(split_1d(1, &knots, &pts, &wts, 1.0).is_err());
    }
}
