//! Least-squares B-spline curve fitting.

pub(crate) mod numerics;
use numerics::{basis_function, solve_linear_system};

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::nurbs_curve::NurbsCurve;

/// **Interpolate** (not approximate) a non-rational B-spline curve
/// through every data point exactly.
///
/// Solves the square N×N basis system `N · CP = data` directly via
/// Gaussian elimination — no normal equations, so no condition-number
/// squaring. With well-spaced data and a moderate degree, mid-curve
/// interpolation error scales as O(h^(d+1)) and is typically far
/// below ULP in 3-D for d=5 with h ~ 1/64.
///
/// Knot vector: clamped, with interior knots placed by the de Boor
/// averaging formula (Piegl & Tiller §9.2.1).
///
/// # Arguments
/// * `points` - Data points to interpolate (at least `degree + 1`).
/// * `params` - Strictly increasing parameter values, one per point.
/// * `degree` - B-spline degree (typically 3 or 5).
///
/// # Errors
/// Returns `CoreError::PreconditionFailed` for invalid inputs or a
/// singular interpolation matrix.
pub fn interpolate_curve(
    points: &[Point3],
    params: &[f64],
    degree: usize,
) -> Result<NurbsCurve, CoreError> {
    let m = points.len();
    if m < degree + 1 {
        return Err(CoreError::PreconditionFailed(format!(
            "interpolate_curve: need at least {} points for degree {}",
            degree + 1,
            degree
        )));
    }
    if params.len() != m {
        return Err(CoreError::PreconditionFailed(
            "params length must match points length".to_string(),
        ));
    }
    let n = m; // num control points = num data points (interpolation)

    // Knot vector: clamped at endpoints, interior by de Boor averaging.
    let num_knots = n + degree + 1;
    let mut knots = vec![0.0; num_knots];
    for k in &mut knots[..=degree] {
        *k = params[0];
    }
    for k in &mut knots[n..] {
        *k = params[m - 1];
    }
    if n > degree + 1 {
        let num_interior = n - degree - 1;
        for j in 1..=num_interior {
            let sum: f64 = params[j..(j + degree)].iter().sum();
            knots[j + degree] = sum / degree as f64;
        }
    }

    // Build the n×n basis matrix N where N[i][j] = B_j(params[i]).
    let mut a: Vec<Vec<f64>> = vec![vec![0.0; n]; n];
    for (i, &t) in params.iter().enumerate() {
        for (j, slot) in a[i].iter_mut().enumerate() {
            *slot = basis_function(j, degree, t, &knots);
        }
    }

    // Solve A · CP = points, separately per (x, y, z) coordinate.
    let bx: Vec<f64> = points.iter().map(|p| p.x).collect();
    let by: Vec<f64> = points.iter().map(|p| p.y).collect();
    let bz: Vec<f64> = points.iter().map(|p| p.z).collect();
    let cx = solve_linear_system(&a, &bx)?;
    let cy = solve_linear_system(&a, &by)?;
    let cz = solve_linear_system(&a, &bz)?;
    let control_points: Vec<Point3> = (0..n).map(|i| Point3::new(cx[i], cy[i], cz[i])).collect();
    let weights = vec![1.0; n];

    Ok(NurbsCurve {
        degree,
        knots,
        control_points,
        weights,
    })
}

/// Fit a non-rational B-spline curve through a sequence of points.
///
/// Endpoints are interpolated exactly. Interior points are approximated
/// in the least-squares sense (Piegl & Tiller, Chapter 9).
///
/// # Arguments
/// * `points` - Data points to fit (at least 2).
/// * `params` - Parameter value for each data point (e.g., chord-length parameterization).
/// * `degree` - Desired polynomial degree (typically 3).
/// * `num_control_points` - Number of control points in the fitted curve.
///
/// # Errors
/// Returns `CoreError::PreconditionFailed` if inputs are invalid.
pub fn fit_curve(
    points: &[Point3],
    params: &[f64],
    degree: usize,
    num_control_points: usize,
) -> Result<NurbsCurve, CoreError> {
    let m = points.len(); // number of data points

    // Validation
    if m < 2 {
        return Err(CoreError::PreconditionFailed(
            "need at least 2 data points for curve fitting".to_string(),
        ));
    }
    if params.len() != m {
        return Err(CoreError::PreconditionFailed(
            "params length must match points length".to_string(),
        ));
    }
    if num_control_points < degree + 1 {
        return Err(CoreError::PreconditionFailed(format!(
            "need at least {} control points for degree {}",
            degree + 1,
            degree
        )));
    }
    if m < num_control_points {
        return Err(CoreError::PreconditionFailed(format!(
            "need at least {} data points for {} control points",
            num_control_points, num_control_points
        )));
    }

    let n = num_control_points; // number of control points

    // Special case: n == 2 and degree >= 1 => just interpolate endpoints
    // (no interior control points to solve for)
    // Also handles: n == degree + 1 (minimum control points)

    // Step 1: Compute knot vector by averaging parameter values (P&T Eq. 9.69)
    let num_knots = n + degree + 1;
    let mut knots = vec![0.0; num_knots];

    // Clamped: first (degree+1) knots = params[0], last (degree+1) knots = params[m-1]
    for k in &mut knots[..=degree] {
        *k = params[0];
    }
    for k in &mut knots[n..] {
        *k = params[m - 1];
    }

    // Interior knots. Two regimes with different textbook rules:
    //  • m == n (interpolation): knot averaging (P&T Eq. 9.8) —
    //    knots[j+degree] = mean(params[j..j+degree]) — satisfies
    //    Schoenberg–Whitney by construction.
    //  • m > n (least-squares approximation): P&T Eq. 9.68/9.69 —
    //    knots distributed over the WHOLE parameter range in
    //    proportion to data density. Using the interpolation rule
    //    here crams every interior knot into the first n params and
    //    the fit runs wild over the rest of the curve (measured: a
    //    513-sample cylinder×cylinder pierce loop fitted with 102
    //    control points deviated 1.4 from the surfaces). Conversely
    //    the approximation rule at m == n drifts off the
    //    Schoenberg–Whitney positions and the interpolation system
    //    degrades (measured: every edge_triple analytic fit failed
    //    its post-construction verify).
    if n > degree + 1 {
        let num_interior = n - degree - 1;
        if m == n {
            for j in 1..=num_interior {
                let sum: f64 = params[j..(j + degree)].iter().sum();
                knots[j + degree] = sum / degree as f64;
            }
        } else {
            let d = m as f64 / (n - degree) as f64;
            for j in 1..=num_interior {
                let pos = j as f64 * d;
                let i = (pos.floor() as usize).clamp(1, m - 1);
                let alpha = pos - i as f64;
                knots[j + degree] = (1.0 - alpha) * params[i - 1] + alpha * params[i];
            }
        }
    }

    // Step 2: Set endpoint control points
    let mut control_points = vec![Point3::new(0.0, 0.0, 0.0); n];
    control_points[0] = points[0];
    control_points[n - 1] = points[m - 1];

    // If there are no interior control points (n == 2), we're done
    if n <= 2 {
        let weights = vec![1.0; n];
        return Ok(NurbsCurve {
            degree,
            knots,
            control_points,
            weights,
        });
    }

    // Step 3: Build basis function matrix N for interior points and interior CPs
    // N is (m-2) x (n-2), stored BANDED: a degree-p basis row has at most
    // p+1 consecutive nonzero values (the active span), so evaluating and
    // accumulating the dense rows×cols matrix wastes rows×cols work on
    // structural zeros. Measured on the wavy-plate offset ladder (PK-329):
    // 240 ms of a 250 ms plane×NURBS SSI solve sat in this fit at ~1000
    // samples × 200 control points; banding removes that term. Skipped
    // entries are EXACTLY 0.0 (Cox–de Boor returns structural zeros
    // outside the support window, padded ±1 against its 1e-14 endpoint
    // capture), and every retained (i,j,k) product accumulates in the
    // same increasing-k order as the dense triple loop, so fitted curves
    // are bit-identical to the dense construction.
    let rows = m - 2; // interior data points
    let cols = n - 2; // interior control points

    // Per-row band: first active interior column and the basis values
    // for the active columns. Also N_0[i] and N_{n-1}[i] for the RHS
    // adjustment (nonzero only when the window reaches an endpoint).
    let mut band_start: Vec<usize> = Vec::with_capacity(rows);
    let mut band_vals: Vec<Vec<f64>> = Vec::with_capacity(rows);
    let mut n0_vec: Vec<f64> = Vec::with_capacity(rows);
    let mut nn_vec: Vec<f64> = Vec::with_capacity(rows);

    for i in 0..rows {
        let t = params[i + 1];
        // Largest knot index with knots[span] <= t, clamped to the
        // non-degenerate span range [degree, n-1].
        let span = match knots.partition_point(|&k| k <= t) {
            0 => degree,
            p => (p - 1).clamp(degree, n - 1),
        };
        // Basis indices with support containing t are [span-degree, span];
        // pad one index each side for the recursion's endpoint capture.
        let lo = span.saturating_sub(degree + 1);
        let hi = (span + 1).min(n - 1);
        // Interior column j holds basis index j+1.
        let col_lo = lo.saturating_sub(1).min(cols.saturating_sub(1));
        let col_hi = hi.saturating_sub(1).min(cols.saturating_sub(1));
        let vals: Vec<f64> = (col_lo..=col_hi)
            .map(|j| basis_function(j + 1, degree, t, &knots))
            .collect();
        band_start.push(col_lo);
        band_vals.push(vals);
        n0_vec.push(if lo == 0 {
            basis_function(0, degree, t, &knots)
        } else {
            0.0
        });
        nn_vec.push(if hi >= n - 1 {
            basis_function(n - 1, degree, t, &knots)
        } else {
            0.0
        });
    }

    // Step 4: Compute RHS: R[i] = points[i+1] - N_0[i]*P[0] - N_{n-1}[i]*P[n-1]
    let mut rhs_x = vec![0.0; rows];
    let mut rhs_y = vec![0.0; rows];
    let mut rhs_z = vec![0.0; rows];

    for i in 0..rows {
        rhs_x[i] =
            points[i + 1].x - n0_vec[i] * control_points[0].x - nn_vec[i] * control_points[n - 1].x;
        rhs_y[i] =
            points[i + 1].y - n0_vec[i] * control_points[0].y - nn_vec[i] * control_points[n - 1].y;
        rhs_z[i] =
            points[i + 1].z - n0_vec[i] * control_points[0].z - nn_vec[i] * control_points[n - 1].z;
    }

    // Step 5: Form normal equations N^T * N and N^T * R over the banded
    // rows. Row-outer accumulation adds each (i,j) product in increasing
    // row order — the same accumulation chain as the dense triple loop,
    // minus terms that are exactly ±0.0 (bit-neutral on these
    // non-negative chains).
    let mut ntn: Vec<Vec<f64>> = vec![vec![0.0; cols]; cols];
    let mut ntr_x = vec![0.0; cols];
    let mut ntr_y = vec![0.0; cols];
    let mut ntr_z = vec![0.0; cols];

    for k in 0..rows {
        let s = band_start[k];
        let vals = &band_vals[k];
        for (di, &vi) in vals.iter().enumerate() {
            let i = s + di;
            for (dj, &vj) in vals.iter().enumerate() {
                ntn[i][s + dj] += vi * vj;
            }
            ntr_x[i] += vi * rhs_x[k];
            ntr_y[i] += vi * rhs_y[k];
            ntr_z[i] += vi * rhs_z[k];
        }
    }

    // Step 6: Solve for each coordinate separately
    let sol_x = solve_linear_system(&ntn, &ntr_x)?;
    let sol_y = solve_linear_system(&ntn, &ntr_y)?;
    let sol_z = solve_linear_system(&ntn, &ntr_z)?;

    // Step 7: Assemble full control point array
    for j in 0..cols {
        control_points[j + 1] = Point3::new(sol_x[j], sol_y[j], sol_z[j]);
    }

    let weights = vec![1.0; n];

    Ok(NurbsCurve {
        degree,
        knots,
        control_points,
        weights,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::CurveEval;

    #[test]
    fn fit_two_points_degree_1() {
        let points = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)];
        let params = vec![0.0, 1.0];
        let curve = fit_curve(&points, &params, 1, 2).unwrap();
        assert_eq!(curve.degree, 1);
        let mid = curve.point_at(0.5);
        assert!((mid.x - 0.5).abs() < 1e-10);
        assert!((mid.y - 0.5).abs() < 1e-10);
    }

    #[test]
    fn fit_line_points_exact() {
        let points: Vec<Point3> = (0..10)
            .map(|i| Point3::new(i as f64 * 0.1, 0.0, 0.0))
            .collect();
        let params: Vec<f64> = (0..10).map(|i| i as f64 / 9.0).collect();
        let curve = fit_curve(&points, &params, 3, 6).unwrap();
        for (&p, &t) in points.iter().zip(params.iter()) {
            let fitted = curve.point_at(t);
            assert!((p.x - fitted.x).abs() < 1e-6, "Deviation at t={t}");
        }
    }

    #[test]
    fn fit_preserves_endpoints() {
        let points: Vec<Point3> = (0..10)
            .map(|i| {
                let t = i as f64 / 9.0;
                Point3::new(t, t * t, 0.0)
            })
            .collect();
        let params: Vec<f64> = (0..10).map(|i| i as f64 / 9.0).collect();
        let curve = fit_curve(&points, &params, 3, 6).unwrap();
        assert!(points[0].distance_to(&curve.point_at(0.0)) < 1e-10);
        assert!(points[9].distance_to(&curve.point_at(1.0 - 1e-12)) < 1e-4);
    }

    #[test]
    fn fit_too_few_points_fails() {
        let points = vec![Point3::new(0.0, 0.0, 0.0)];
        let params = vec![0.0];
        assert!(fit_curve(&points, &params, 3, 4).is_err());
    }

    #[test]
    fn fit_quadratic_curve() {
        // y = x^2 sampled
        let n = 20;
        let points: Vec<Point3> = (0..n)
            .map(|i| {
                let x = i as f64 / (n - 1) as f64;
                Point3::new(x, x * x, 0.0)
            })
            .collect();
        let params: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        let curve = fit_curve(&points, &params, 3, 8).unwrap();
        // Check deviation at sample points
        for (&p, &t) in points.iter().zip(params.iter()) {
            let fitted = curve.point_at(t);
            assert!(p.distance_to(&fitted) < 0.05, "Large deviation at t={t}");
        }
    }

    #[test]
    fn fit_mismatched_params_length_fails() {
        let points = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
        let params = vec![0.0]; // wrong length
        assert!(fit_curve(&points, &params, 1, 2).is_err());
    }

    #[test]
    fn fit_too_few_control_points_fails() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];
        let params = vec![0.0, 0.5, 1.0];
        // degree 3 needs at least 4 control points
        assert!(fit_curve(&points, &params, 3, 3).is_err());
    }

    #[test]
    fn fit_more_control_points_than_data_fails() {
        let points = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ];
        let params = vec![0.0, 0.5, 1.0];
        // 5 control points but only 3 data points
        assert!(fit_curve(&points, &params, 1, 5).is_err());
    }

    #[test]
    fn fit_empty_system_solve() {
        // Test solve_linear_system with empty input
        let result = solve_linear_system(&[], &[]);
        assert!(result.is_ok());
        assert!(result.unwrap().is_empty());
    }

    #[test]
    fn fit_basis_function_last_span() {
        // Test the special case where t == knots[i+1] at the last span
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        // At t=1.0, last span, degree 0
        let val = basis_function(2, 0, 1.0, &knots);
        assert!(
            (val - 1.0).abs() < 1e-10,
            "basis at last knot should be 1.0"
        );
    }

    #[test]
    fn fit_basis_function_middle_span() {
        // Test that t in a middle span returns non-zero
        let knots = vec![0.0, 0.0, 0.5, 1.0, 1.0];
        let val = basis_function(0, 1, 0.25, &knots);
        assert!(val > 0.0);
    }
}
