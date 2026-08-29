//! Least-squares degree-3 B-spline fitting for sampled offset surfaces.
//!
//! Replaces the former degree-1 "control points = samples" construction:
//! exact offset samples on a uniform parameter grid are approximated by
//! a degree-3 tensor-product B-spline via the standard two-pass global
//! least-squares fit (rows in u, then the intermediate control points in
//! v). The fit itself claims nothing about accuracy — the caller
//! measures deviation against independent dense samples and gates on it.

use popkern_core::geometry::{NurbsSurface, Point3, Vec3};

use super::error::OffsetError;

/// Fitted B-spline degree in both directions.
pub(crate) const FIT_DEGREE: usize = 3;

/// Build a uniform clamped knot vector for `n_pts` control points at the
/// given `degree`, spanning `[start, end]`.
pub(crate) fn clamped_knots(n_pts: usize, degree: usize, start: f64, end: f64) -> Vec<f64> {
    let n_knots = n_pts + degree + 1;
    let n_interior = n_pts - degree; // >= 1
    let mut knots = Vec::with_capacity(n_knots);
    for _ in 0..=degree {
        knots.push(start);
    }
    for i in 1..n_interior {
        let t = start + (end - start) * (i as f64) / (n_interior as f64);
        knots.push(t);
    }
    for _ in 0..=degree {
        knots.push(end);
    }
    knots
}

/// Find the knot span index for parameter `t` (Piegl & Tiller A2.1).
fn find_span(knots: &[f64], degree: usize, n_pts: usize, t: f64) -> usize {
    if t >= knots[n_pts] {
        return n_pts - 1;
    }
    if t <= knots[degree] {
        return degree;
    }
    let mut low = degree;
    let mut high = n_pts;
    let mut mid = (low + high) / 2;
    while t < knots[mid] || t >= knots[mid + 1] {
        if t < knots[mid] {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    mid
}

/// Evaluate all non-vanishing B-spline basis functions at `t` and place
/// them into a dense length-`n_pts` row (Cox–de Boor, P&T A2.2).
fn basis_row(knots: &[f64], degree: usize, n_pts: usize, t: f64) -> Vec<f64> {
    let span = find_span(knots, degree, n_pts, t);
    let mut n = vec![0.0_f64; degree + 1];
    let mut left = vec![0.0_f64; degree + 1];
    let mut right = vec![0.0_f64; degree + 1];
    n[0] = 1.0;
    for j in 1..=degree {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom != 0.0 { n[r] / denom } else { 0.0 };
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    let mut row = vec![0.0_f64; n_pts];
    row[(span - degree)..=span].copy_from_slice(&n);
    row
}

/// Cholesky factorization of a symmetric positive-definite matrix
/// (row-major, in place: lower triangle receives L).
fn cholesky_factor(a: &mut [f64], n: usize) -> Result<(), OffsetError> {
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i * n + j];
            for k in 0..j {
                sum -= a[i * n + k] * a[j * n + k];
            }
            if i == j {
                if sum <= 0.0 {
                    return Err(OffsetError::FitFailed {
                        reason: format!("normal equations not positive definite at row {i}"),
                    });
                }
                a[i * n + i] = sum.sqrt();
            } else {
                a[i * n + j] = sum / a[j * n + j];
            }
        }
    }
    Ok(())
}

/// Solve `L Lᵀ x = b` for a point-valued right-hand side, in place.
fn cholesky_solve(l: &[f64], n: usize, b: &mut [Vec3]) {
    for i in 0..n {
        let mut x = b[i];
        for k in 0..i {
            x = x - b[k] * l[i * n + k];
        }
        b[i] = x * (1.0 / l[i * n + i]);
    }
    for i in (0..n).rev() {
        let mut x = b[i];
        for k in (i + 1)..n {
            x = x - b[k] * l[k * n + i];
        }
        b[i] = x * (1.0 / l[i * n + i]);
    }
}

/// One-directional least-squares setup: basis rows for every sample
/// parameter plus the Cholesky-factored normal matrix `BᵀB`.
struct DirectionFit {
    rows: Vec<Vec<f64>>,
    factored: Vec<f64>,
    n_cp: usize,
}

impl DirectionFit {
    fn new(params: &[f64], knots: &[f64], degree: usize, n_cp: usize) -> Result<Self, OffsetError> {
        let rows: Vec<Vec<f64>> = params
            .iter()
            .map(|&t| basis_row(knots, degree, n_cp, t))
            .collect();
        let mut ata = vec![0.0_f64; n_cp * n_cp];
        for row in &rows {
            for i in 0..n_cp {
                if row[i] == 0.0 {
                    continue;
                }
                for j in 0..n_cp {
                    ata[i * n_cp + j] += row[i] * row[j];
                }
            }
        }
        cholesky_factor(&mut ata, n_cp)?;
        Ok(Self {
            rows,
            factored: ata,
            n_cp,
        })
    }

    /// Least-squares fit of one sample sequence (length = params).
    fn solve(&self, samples: &[Vec3]) -> Vec<Vec3> {
        let mut rhs = vec![Vec3::zero(); self.n_cp];
        for (row, s) in self.rows.iter().zip(samples) {
            for i in 0..self.n_cp {
                if row[i] != 0.0 {
                    rhs[i] = rhs[i] + *s * row[i];
                }
            }
        }
        cholesky_solve(&self.factored, self.n_cp, &mut rhs);
        rhs
    }
}

/// Fit a degree-3 tensor-product B-spline through exact offset samples.
///
/// `samples` is row-major `[i_u][i_v]` over the uniform parameter grids
/// `u_params × v_params`. Returns a non-rational (all weights 1) NURBS
/// surface with `n_u × n_v` control points on uniform clamped knots.
///
/// # Errors
///
/// [`OffsetError::FitFailed`] when the normal equations degenerate
/// (requires `n_u, n_v > FIT_DEGREE` and at least as many samples as
/// control points per direction).
pub(crate) fn fit_degree3_surface(
    u_params: &[f64],
    v_params: &[f64],
    samples: &[Point3],
    n_u: usize,
    n_v: usize,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> Result<NurbsSurface, OffsetError> {
    let m_u = u_params.len();
    let m_v = v_params.len();
    if n_u <= FIT_DEGREE || n_v <= FIT_DEGREE || m_u < n_u || m_v < n_v {
        return Err(OffsetError::FitFailed {
            reason: format!("underdetermined fit: {m_u}×{m_v} samples for {n_u}×{n_v} controls"),
        });
    }
    let knots_u = clamped_knots(n_u, FIT_DEGREE, u_range.0, u_range.1);
    let knots_v = clamped_knots(n_v, FIT_DEGREE, v_range.0, v_range.1);
    let fit_u = DirectionFit::new(u_params, &knots_u, FIT_DEGREE, n_u)?;
    let fit_v = DirectionFit::new(v_params, &knots_v, FIT_DEGREE, n_v)?;

    // Pass 1: fit each v-column of samples along u → n_u × m_v
    // intermediate points.
    let origin = samples[0];
    let as_vec = |p: &Point3| *p - origin;
    let mut intermediate = vec![Vec3::zero(); n_u * m_v];
    let mut column = vec![Vec3::zero(); m_u];
    for j in 0..m_v {
        for i in 0..m_u {
            column[i] = as_vec(&samples[i * m_v + j]);
        }
        for (k, cp) in fit_u.solve(&column).into_iter().enumerate() {
            intermediate[k * m_v + j] = cp;
        }
    }

    // Pass 2: fit each intermediate row along v → n_u × n_v controls.
    let mut control_points = vec![Point3::new(0.0, 0.0, 0.0); n_u * n_v];
    let mut row = vec![Vec3::zero(); m_v];
    for k in 0..n_u {
        row.copy_from_slice(&intermediate[k * m_v..(k + 1) * m_v]);
        for (l, cp) in fit_v.solve(&row).into_iter().enumerate() {
            control_points[k * n_v + l] = origin + cp;
        }
    }

    Ok(NurbsSurface {
        degree_u: FIT_DEGREE,
        degree_v: FIT_DEGREE,
        knots_u,
        knots_v,
        control_points,
        num_u: n_u,
        num_v: n_v,
        weights: vec![1.0; n_u * n_v],
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn uniform(range: (f64, f64), m: usize) -> Vec<f64> {
        (0..m)
            .map(|i| range.0 + (range.1 - range.0) * (i as f64) / ((m - 1) as f64))
            .collect()
    }

    #[test]
    fn basis_rows_partition_unity() {
        let n = 8;
        let knots = clamped_knots(n, FIT_DEGREE, 0.0, 2.0);
        for &t in &[0.0, 0.3, 1.0, 1.999, 2.0] {
            let row = basis_row(&knots, FIT_DEGREE, n, t);
            let sum: f64 = row.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12, "partition of unity at t={t}");
        }
    }

    #[test]
    fn fit_reproduces_bilinear_data_exactly() {
        // Degree-3 spline space contains all bilinear functions, so the
        // least-squares fit of samples from z = x + 2y is exact.
        let u = uniform((0.0, 1.0), 16);
        let v = uniform((0.0, 1.0), 16);
        let samples: Vec<Point3> = u
            .iter()
            .flat_map(|&x| v.iter().map(move |&y| Point3::new(x, y, x + 2.0 * y)))
            .collect();
        let surf = fit_degree3_surface(&u, &v, &samples, 8, 8, (0.0, 1.0), (0.0, 1.0)).unwrap();
        use popkern_core::geometry::surface::SurfaceEval;
        for &x in &[0.05, 0.4, 0.93] {
            for &y in &[0.11, 0.5, 0.87] {
                let p = surf.point_at(x, y);
                assert!(
                    (p.z - (x + 2.0 * y)).abs() < 1e-10,
                    "bilinear reproduction at ({x},{y}): {}",
                    p.z
                );
            }
        }
    }

    #[test]
    fn underdetermined_fit_declines() {
        let u = uniform((0.0, 1.0), 4);
        let v = uniform((0.0, 1.0), 4);
        let samples = vec![Point3::new(0.0, 0.0, 0.0); 16];
        let result = fit_degree3_surface(&u, &v, &samples, 8, 8, (0.0, 1.0), (0.0, 1.0));
        assert!(matches!(result, Err(OffsetError::FitFailed { .. })));
    }
}
