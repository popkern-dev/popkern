//! Knot insertion via Boehm's algorithm.

use super::{Homogeneous, find_span};
use crate::geometry::Point3;

/// Insert a single knot at parameter `t` into a 1D NURBS.
/// Returns (new_knots, new_control_points, new_weights).
/// The curve shape is mathematically unchanged.
#[allow(dead_code)]
pub(crate) fn insert_knot_1d(
    degree: usize,
    knots: &[f64],
    control_points: &[Point3],
    weights: &[f64],
    t: f64,
) -> (Vec<f64>, Vec<Point3>, Vec<f64>) {
    let n_controls = control_points.len();
    let p = degree;
    let k = find_span(p, knots, n_controls, t);

    // Convert to homogeneous coordinates
    let hom: Vec<Homogeneous> = control_points
        .iter()
        .zip(weights.iter())
        .map(|(pt, &w)| Homogeneous::from_point_weight(pt, w))
        .collect();

    // Build new knot vector: insert t at position k+1
    let mut new_knots = Vec::with_capacity(knots.len() + 1);
    new_knots.extend_from_slice(&knots[..=k]);
    new_knots.push(t);
    new_knots.extend_from_slice(&knots[k + 1..]);

    // Compute new control points in homogeneous space
    let mut new_hom = Vec::with_capacity(n_controls + 1);

    for i in 0..=n_controls {
        if i <= k.saturating_sub(p) {
            // Unchanged: copy directly
            new_hom.push(hom[i]);
        } else if i > k {
            // Shifted: take from previous index
            new_hom.push(hom[i - 1]);
        } else {
            // Blended: k-p+1 <= i <= k
            let denom = knots[i + p] - knots[i];
            let alpha = if denom.abs() < 1e-30 {
                0.0
            } else {
                (t - knots[i]) / denom
            };
            new_hom.push(Homogeneous::lerp(&hom[i - 1], &hom[i], alpha));
        }
    }

    // Convert back from homogeneous
    let mut new_pts = Vec::with_capacity(new_hom.len());
    let mut new_wts = Vec::with_capacity(new_hom.len());
    for h in &new_hom {
        let (pt, w) = h.to_point3_and_weight();
        new_pts.push(pt);
        new_wts.push(w);
    }

    (new_knots, new_pts, new_wts)
}

/// Insert multiple sorted knots (repeated single insertion).
#[allow(dead_code)]
pub(crate) fn insert_knots_1d(
    degree: usize,
    knots: &[f64],
    control_points: &[Point3],
    weights: &[f64],
    new_knots_to_insert: &[f64],
) -> (Vec<f64>, Vec<Point3>, Vec<f64>) {
    let mut k = knots.to_vec();
    let mut pts = control_points.to_vec();
    let mut wts = weights.to_vec();
    for &t in new_knots_to_insert {
        let (nk, np, nw) = insert_knot_1d(degree, &k, &pts, &wts, t);
        k = nk;
        pts = np;
        wts = nw;
    }
    (k, pts, wts)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{CurveEval, NurbsCurve};

    const TOL: f64 = 1e-10;

    fn points_close(a: &Point3, b: &Point3) -> bool {
        let dx = a.x - b.x;
        let dy = a.y - b.y;
        let dz = a.z - b.z;
        (dx * dx + dy * dy + dz * dz).sqrt() < TOL
    }

    /// Helper: evaluate a NurbsCurve constructed from raw data at param t.
    fn eval_raw(degree: usize, knots: &[f64], pts: &[Point3], wts: &[f64], t: f64) -> Point3 {
        let curve = NurbsCurve {
            degree,
            knots: knots.to_vec(),
            control_points: pts.to_vec(),
            weights: wts.to_vec(),
        };
        curve.point_at(t)
    }

    #[test]
    fn degree1_line_preserves_shape() {
        // Degree-1 line from (0,0,0) to (1,1,0)
        let degree = 1;
        let knots = vec![0.0, 0.0, 1.0, 1.0];
        let pts = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)];
        let wts = vec![1.0, 1.0];

        let (nk, np, nw) = insert_knot_1d(degree, &knots, &pts, &wts, 0.5);

        // Check sizes
        assert_eq!(nk.len(), knots.len() + 1);
        assert_eq!(np.len(), pts.len() + 1);
        assert_eq!(nw.len(), wts.len() + 1);

        // Evaluate at several parameters and compare
        for &t in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let before = eval_raw(degree, &knots, &pts, &wts, t);
            let after = eval_raw(degree, &nk, &np, &nw, t);
            assert!(
                points_close(&before, &after),
                "Mismatch at t={t}: before={before:?}, after={after:?}"
            );
        }
    }

    #[test]
    fn degree2_bezier_preserves_shape() {
        // Degree-2 Bezier: 3 control points
        let degree = 2;
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0];

        let (nk, np, nw) = insert_knot_1d(degree, &knots, &pts, &wts, 0.5);

        assert_eq!(nk.len(), knots.len() + 1);
        assert_eq!(np.len(), pts.len() + 1);
        assert_eq!(nw.len(), wts.len() + 1);

        for &t in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let before = eval_raw(degree, &knots, &pts, &wts, t);
            let after = eval_raw(degree, &nk, &np, &nw, t);
            assert!(
                points_close(&before, &after),
                "Mismatch at t={t}: before={before:?}, after={after:?}"
            );
        }
    }

    #[test]
    fn knot_and_control_point_counts() {
        let degree = 2;
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0];

        let (nk, np, nw) = insert_knot_1d(degree, &knots, &pts, &wts, 0.3);

        assert_eq!(nk.len(), knots.len() + 1);
        assert_eq!(np.len(), pts.len() + 1);
        assert_eq!(nw.len(), wts.len() + 1);

        // Knot vector should still satisfy n_knots = n_controls + degree + 1
        assert_eq!(nk.len(), np.len() + degree + 1);
    }

    #[test]
    fn multiple_knot_insertion() {
        let degree = 2;
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let pts = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];
        let wts = vec![1.0, 1.0, 1.0];

        let to_insert = vec![0.25, 0.5, 0.75];
        let (nk, np, nw) = insert_knots_1d(degree, &knots, &pts, &wts, &to_insert);

        // 3 insertions = 3 extra knots and control points
        assert_eq!(nk.len(), knots.len() + 3);
        assert_eq!(np.len(), pts.len() + 3);
        assert_eq!(nw.len(), wts.len() + 3);

        // Shape preserved
        for &t in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
            let before = eval_raw(degree, &knots, &pts, &wts, t);
            let after = eval_raw(degree, &nk, &np, &nw, t);
            assert!(
                points_close(&before, &after),
                "Mismatch at t={t}: before={before:?}, after={after:?}"
            );
        }
    }

    #[test]
    fn rational_quarter_circle() {
        // Quarter circle NURBS (degree 2, rational)
        // Control points: (1,0,0), (1,1,0), (0,1,0) with weights 1, 1/sqrt(2), 1
        let degree = 2;
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let pts = vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        let w_mid = 1.0 / std::f64::consts::SQRT_2;
        let wts = vec![1.0, w_mid, 1.0];

        let (nk, np, nw) = insert_knot_1d(degree, &knots, &pts, &wts, 0.5);

        assert_eq!(nk.len(), knots.len() + 1);
        assert_eq!(np.len(), pts.len() + 1);
        assert_eq!(nw.len(), wts.len() + 1);

        // Shape preserved — check that points lie on unit circle
        for &t in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let before = eval_raw(degree, &knots, &pts, &wts, t);
            let after = eval_raw(degree, &nk, &np, &nw, t);
            assert!(
                points_close(&before, &after),
                "Mismatch at t={t}: before={before:?}, after={after:?}"
            );
            // Verify it's on the unit circle (x^2 + y^2 ≈ 1)
            let r2 = after.x * after.x + after.y * after.y;
            assert!(
                (r2 - 1.0).abs() < 1e-8,
                "Point not on unit circle at t={t}: r^2={r2}"
            );
        }
    }
}
