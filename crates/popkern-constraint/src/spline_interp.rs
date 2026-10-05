//! B-spline interpolation for through-point spline creation.
//!
//! Given a set of through-points, computes control points for a cubic
//! B-spline that interpolates all of them.  Uses chord-length
//! parameterization and the knot averaging method.

use faer::Mat;

/// Result of B-spline interpolation: (control_points, knot_vector, parameter_values).
type InterpResult = Result<(Vec<(f64, f64)>, Vec<f64>, Vec<f64>), String>;

/// Compute control points for a cubic B-spline interpolating the given points.
///
/// Returns `(control_point_positions, knot_vector, parameter_values)`.
///
/// # Algorithm
/// 1. Chord-length parameterization
/// 2. Knot averaging for interior knots
/// 3. Build interpolation matrix (basis function values at each parameter)
/// 4. LU solve for control point coordinates
pub(crate) fn interpolate_cubic(through_points: &[(f64, f64)]) -> InterpResult {
    let n = through_points.len();
    if n < 2 {
        return Err("need at least 2 through-points".into());
    }
    if n == 2 {
        return interpolate_two_points(through_points);
    }

    let degree = 3;

    // If fewer points than degree+1, reduce to what we can handle.
    if n < degree + 1 {
        return interpolate_low_order(through_points);
    }

    // Step 1: chord-length parameterization.
    let t_values = chord_length_params(through_points);

    // Step 2: clamped knot vector via averaging method.
    let knots = averaging_knots(&t_values, n, degree);

    // Step 3: build interpolation matrix.
    let basis_mat = build_interpolation_matrix(&t_values, &knots, n, degree);

    // Step 4: solve for control points.
    let control_points = solve_interpolation(&basis_mat, through_points, n)?;

    Ok((control_points, knots, t_values))
}

/// Compute control points for a periodic cubic B-spline interpolating
/// the given points.
///
/// Returns `(control_point_positions, knot_vector, parameter_values)`.
///
/// The spline wraps around: the last point connects smoothly back to the
/// first.  Uses uniform parameterization and a periodic knot vector.
pub(crate) fn interpolate_cubic_periodic(through_points: &[(f64, f64)]) -> InterpResult {
    let n = through_points.len();
    if n < 3 {
        return Err("need at least 3 through-points for periodic spline".into());
    }

    let degree = 3;

    // Uniform parameterization: t_i = i (period = n).
    let t_values: Vec<f64> = (0..n).map(|i| i as f64).collect();

    // Uniform periodic knot vector: [0, 1, 2, ..., n].
    let knots: Vec<f64> = (0..=n).map(|i| i as f64).collect();

    // Build the periodic interpolation matrix using the extended knot vector.
    let ext_knots = build_periodic_ext_knots(&knots, n, degree);
    let n_ext = n + degree;

    let mut mat = Mat::<f64>::zeros(n, n);
    for (row, &ti) in t_values.iter().enumerate() {
        let (span, tc) = crate::equations::spline_eval::find_span(&ext_knots, degree, n_ext, ti);
        let basis = crate::equations::spline_eval::basis_functions(&ext_knots, span, degree, tc);
        for (k, &nk) in basis.iter().enumerate().take(degree + 1) {
            let raw_idx = span - degree + k;
            let col = raw_idx % n;
            mat.write(row, col, mat.read(row, col) + nk);
        }
    }

    // Solve M * Px = Dx, M * Py = Dy.

    let mut rhs_x = Mat::<f64>::zeros(n, 1);
    let mut rhs_y = Mat::<f64>::zeros(n, 1);
    for (i, &(px, py)) in through_points.iter().enumerate() {
        rhs_x.write(i, 0, px);
        rhs_y.write(i, 0, py);
    }

    let sol_x = crate::dense_solve::lu(&mat, &rhs_x);
    let sol_y = crate::dense_solve::lu(&mat, &rhs_y);

    let mut control_points = Vec::with_capacity(n);
    for i in 0..n {
        control_points.push((sol_x.read(i, 0), sol_y.read(i, 0)));
    }

    Ok((control_points, knots, t_values))
}

// ── Internal helpers ─────────────────────────────────────────────────────

/// Chord-length parameterization: t_0=0, t_N=1.
fn chord_length_params(pts: &[(f64, f64)]) -> Vec<f64> {
    let n = pts.len();
    let mut dists = vec![0.0; n];
    for i in 1..n {
        let dx = pts[i].0 - pts[i - 1].0;
        let dy = pts[i].1 - pts[i - 1].1;
        dists[i] = dists[i - 1] + (dx * dx + dy * dy).sqrt();
    }
    let total = dists[n - 1];
    if total < f64::EPSILON {
        // Degenerate: all points coincident, use uniform.
        return (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
    }
    dists.iter().map(|d| d / total).collect()
}

/// Knot averaging method for interior knots.
fn averaging_knots(t_values: &[f64], n: usize, degree: usize) -> Vec<f64> {
    let m = n + degree + 1;
    let mut knots = Vec::with_capacity(m);

    // First degree+1 knots = 0.
    knots.extend(std::iter::repeat_n(0.0, degree + 1));

    // Interior knots: j = 1..n-degree.
    let n_interior = n.saturating_sub(degree);
    for j in 1..n_interior {
        let sum: f64 = t_values[j..(j + degree)].iter().sum();
        knots.push(sum / degree as f64);
    }

    // Last degree+1 knots = 1.
    knots.extend(std::iter::repeat_n(1.0, degree + 1));

    knots
}

/// Build the N x N interpolation matrix: M[i][j] = N_j(t_i).
fn build_interpolation_matrix(
    t_values: &[f64],
    knots: &[f64],
    n: usize,
    degree: usize,
) -> Mat<f64> {
    let mut mat = Mat::<f64>::zeros(n, n);
    for (row, &ti) in t_values.iter().enumerate() {
        let (span, tc) = crate::equations::spline_eval::find_span(knots, degree, n, ti);
        let basis = crate::equations::spline_eval::basis_functions(knots, span, degree, tc);
        for (k, &nk) in basis.iter().enumerate().take(degree + 1) {
            let col = span - degree + k;
            if col < n {
                mat.write(row, col, mat.read(row, col) + nk);
            }
        }
    }
    mat
}

/// Solve the interpolation system M * P = D.
fn solve_interpolation(
    mat: &Mat<f64>,
    through_points: &[(f64, f64)],
    n: usize,
) -> Result<Vec<(f64, f64)>, String> {
    let mut rhs_x = Mat::<f64>::zeros(n, 1);
    let mut rhs_y = Mat::<f64>::zeros(n, 1);
    for (i, &(px, py)) in through_points.iter().enumerate() {
        rhs_x.write(i, 0, px);
        rhs_y.write(i, 0, py);
    }

    let sol_x = crate::dense_solve::lu(mat, &rhs_x);
    let sol_y = crate::dense_solve::lu(mat, &rhs_y);

    let mut control_points = Vec::with_capacity(n);
    for i in 0..n {
        control_points.push((sol_x.read(i, 0), sol_y.read(i, 0)));
    }
    Ok(control_points)
}

/// Handle the degenerate case of exactly 2 through-points.
///
/// For a cubic with 2 points we need 4 control points. Place them
/// along the straight line to produce a linear segment.
fn interpolate_two_points(pts: &[(f64, f64)]) -> InterpResult {
    let (x0, y0) = pts[0];
    let (x1, y1) = pts[1];
    let control_points = vec![
        (x0, y0),
        (x0 + (x1 - x0) / 3.0, y0 + (y1 - y0) / 3.0),
        (x0 + 2.0 * (x1 - x0) / 3.0, y0 + 2.0 * (y1 - y0) / 3.0),
        (x1, y1),
    ];
    let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
    let t_values = vec![0.0, 1.0];
    Ok((control_points, knots, t_values))
}

/// Handle 3 through-points with a cubic (degree 3, 3 CPs is not
/// enough for a cubic, so we use degree=2 quadratic and upsample).
///
/// Actually with 3 through-points and degree 3, we need 3 CPs but
/// a cubic needs at least 4 CPs (degree+1). So we interpolate with
/// degree = min(degree, n-1) = 2, then return as-is.
fn interpolate_low_order(through_points: &[(f64, f64)]) -> InterpResult {
    let n = through_points.len();
    let degree = n - 1; // n=3 -> degree=2

    let t_values = chord_length_params(through_points);
    let knots = averaging_knots(&t_values, n, degree);
    let basis_mat = build_interpolation_matrix(&t_values, &knots, n, degree);
    let control_points = solve_interpolation(&basis_mat, through_points, n)?;

    Ok((control_points, knots, t_values))
}

/// Build extended periodic knot vector (same as in spline_eval).
fn build_periodic_ext_knots(knots: &[f64], n: usize, degree: usize) -> Vec<f64> {
    crate::equations::spline_eval::periodic_knot_vector(knots, n, degree)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interp_collinear_3_points() {
        let pts = vec![(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)];
        let (cps, _knots, t_values) = interpolate_cubic(&pts).unwrap();

        // With collinear points, the resulting spline should be a straight line.
        // Verify by evaluating at the through-point parameters.
        let mut sk = crate::Sketch::new();
        let mut handles = Vec::new();
        for &(x, y) in &cps {
            handles.push(sk.add_point(x, y));
        }
        let degree = if cps.len() >= 4 { 3 } else { cps.len() - 1 };
        let sh = sk.add_spline(&handles, degree);
        // Override knots with the interpolation knots.
        sk.splines[sh].knots = _knots;

        let sp = &sk.splines[sh];
        for (i, &ti) in t_values.iter().enumerate() {
            let r = crate::equations::spline_eval::eval_point(sp, ti, &sk);
            let tol = 1e-10;
            assert!(
                (r.x - pts[i].0).abs() < tol && (r.y - pts[i].1).abs() < tol,
                "point {} at t={}: got ({}, {}), expected ({}, {})",
                i,
                ti,
                r.x,
                r.y,
                pts[i].0,
                pts[i].1
            );
        }

        // Also check midpoint is on the line y=x.
        let t_mid = t_values[1];
        let r = crate::equations::spline_eval::eval_point(sp, t_mid, &sk);
        assert!(
            (r.x - r.y).abs() < 1e-10,
            "midpoint should be on y=x: ({}, {})",
            r.x,
            r.y
        );
    }

    #[test]
    fn interp_4_points_passes_through() {
        let pts = vec![(0.0, 0.0), (1.0, 2.0), (3.0, 1.0), (4.0, 0.0)];
        let (cps, knots, t_values) = interpolate_cubic(&pts).unwrap();

        assert_eq!(cps.len(), 4);

        let mut sk = crate::Sketch::new();
        let mut handles = Vec::new();
        for &(x, y) in &cps {
            handles.push(sk.add_point(x, y));
        }
        let sh = sk.add_spline(&handles, 3);
        sk.splines[sh].knots = knots;

        let sp = &sk.splines[sh];
        for (i, &ti) in t_values.iter().enumerate() {
            let r = crate::equations::spline_eval::eval_point(sp, ti, &sk);
            let tol = 1e-10;
            assert!(
                (r.x - pts[i].0).abs() < tol && (r.y - pts[i].1).abs() < tol,
                "point {} at t={}: got ({}, {}), expected ({}, {})",
                i,
                ti,
                r.x,
                r.y,
                pts[i].0,
                pts[i].1
            );
        }
    }

    #[test]
    fn interp_6_points_passes_through() {
        let pts = vec![
            (0.0, 0.0),
            (1.0, 3.0),
            (2.0, 1.0),
            (3.0, 4.0),
            (4.0, 2.0),
            (5.0, 0.0),
        ];
        let (cps, knots, t_values) = interpolate_cubic(&pts).unwrap();

        assert_eq!(cps.len(), 6);
        assert_eq!(knots.len(), 6 + 3 + 1); // n + degree + 1

        let mut sk = crate::Sketch::new();
        let mut handles = Vec::new();
        for &(x, y) in &cps {
            handles.push(sk.add_point(x, y));
        }
        let sh = sk.add_spline(&handles, 3);
        sk.splines[sh].knots = knots;

        let sp = &sk.splines[sh];
        for (i, &ti) in t_values.iter().enumerate() {
            let r = crate::equations::spline_eval::eval_point(sp, ti, &sk);
            let tol = 1e-9;
            assert!(
                (r.x - pts[i].0).abs() < tol && (r.y - pts[i].1).abs() < tol,
                "point {} at t={}: got ({}, {}), expected ({}, {})",
                i,
                ti,
                r.x,
                r.y,
                pts[i].0,
                pts[i].1
            );
        }
    }

    #[test]
    fn interp_2_points() {
        let pts = vec![(0.0, 0.0), (4.0, 3.0)];
        let (cps, knots, t_values) = interpolate_cubic(&pts).unwrap();
        assert_eq!(cps.len(), 4);
        assert_eq!(t_values.len(), 2);
        assert!((t_values[0] - 0.0).abs() < 1e-12);
        assert!((t_values[1] - 1.0).abs() < 1e-12);

        // All control points should be on the line from (0,0) to (4,3).
        for &(cx, cy) in &cps {
            let cross = cx * 3.0 - cy * 4.0;
            assert!(cross.abs() < 1e-10, "cp ({cx}, {cy}) not on line");
        }

        // Verify endpoints.
        let mut sk = crate::Sketch::new();
        let mut handles = Vec::new();
        for &(x, y) in &cps {
            handles.push(sk.add_point(x, y));
        }
        let sh = sk.add_spline(&handles, 3);
        sk.splines[sh].knots = knots;

        let sp = &sk.splines[sh];
        let r0 = crate::equations::spline_eval::eval_point(sp, 0.0, &sk);
        let r1 = crate::equations::spline_eval::eval_point(sp, 1.0, &sk);
        assert!((r0.x).abs() < 1e-10 && (r0.y).abs() < 1e-10);
        assert!((r1.x - 4.0).abs() < 1e-10 && (r1.y - 3.0).abs() < 1e-10);
    }

    #[test]
    fn interp_periodic_triangle() {
        let pts = vec![(1.0, 0.0), (-0.5, 0.866), (-0.5, -0.866)];
        let (cps, knots, t_values) = interpolate_cubic_periodic(&pts).unwrap();
        assert_eq!(cps.len(), 3);
        assert_eq!(knots.len(), 4); // n + 1

        // Verify the spline passes through the points.
        let mut sk = crate::Sketch::new();
        let mut handles = Vec::new();
        for &(x, y) in &cps {
            handles.push(sk.add_point(x, y));
        }
        let sh = sk.add_spline(&handles, 3);
        sk.splines[sh].knots = knots;
        sk.splines[sh].periodic = true;

        let sp = &sk.splines[sh];
        for (i, &ti) in t_values.iter().enumerate() {
            let r = crate::equations::spline_eval::eval_point(sp, ti, &sk);
            let tol = 1e-9;
            assert!(
                (r.x - pts[i].0).abs() < tol && (r.y - pts[i].1).abs() < tol,
                "periodic point {} at t={}: got ({}, {}), expected ({}, {})",
                i,
                ti,
                r.x,
                r.y,
                pts[i].0,
                pts[i].1
            );
        }
    }

    #[test]
    fn interp_error_on_too_few() {
        assert!(interpolate_cubic(&[]).is_err());
        assert!(interpolate_cubic(&[(1.0, 2.0)]).is_err());
        assert!(interpolate_cubic_periodic(&[]).is_err());
        assert!(interpolate_cubic_periodic(&[(0.0, 0.0), (1.0, 1.0)]).is_err());
    }
}
