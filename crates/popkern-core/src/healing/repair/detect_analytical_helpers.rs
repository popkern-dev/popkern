//! Helper functions for analytical surface detection.
//!
//! Contains axis estimation, linear algebra, circle fitting, and
//! center estimation used by the detection routines.

use crate::geometry::{Point3, Vec3};

/// Solve a 3x3 linear system `A * x = b` via Cramer's rule.
///
/// Returns `None` if the determinant is near-zero.
fn solve_3x3(a: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() < 1e-30 {
        return None;
    }
    let inv = 1.0 / det;
    let x = (b[0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (b[1] * a[2][2] - a[1][2] * b[2])
        + a[0][2] * (b[1] * a[2][1] - a[1][1] * b[2]))
        * inv;
    let y = (a[0][0] * (b[1] * a[2][2] - a[1][2] * b[2])
        - b[0] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * b[2] - b[1] * a[2][0]))
        * inv;
    let z = (a[0][0] * (a[1][1] * b[2] - b[1] * a[2][1])
        - a[0][1] * (a[1][0] * b[2] - b[1] * a[2][0])
        + b[0] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
        * inv;
    Some([x, y, z])
}

/// Find the eigenvector corresponding to the smallest eigenvalue of a 3x3
/// symmetric matrix, using inverse power iteration.
fn smallest_eigenvector_3x3(cov: [[f64; 3]; 3]) -> Option<Vec3> {
    let mut v = [1.0, 0.5, 0.25];
    for _ in 0..20 {
        let w = solve_3x3(cov, v)?;
        let len = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
        if len < 1e-30 {
            return None;
        }
        v = [w[0] / len, w[1] / len, w[2] / len];
    }
    Some(Vec3::new(v[0], v[1], v[2]))
}

/// Estimate the symmetry axis from point positions via covariance eigenvectors.
///
/// For rotationally symmetric surfaces (cylinder, cone, torus), the axis
/// is the direction of least variance in the point cloud — the smallest
/// eigenvector of the position covariance matrix.
pub(super) fn estimate_axis_from_positions(points: &[Point3]) -> Option<Vec3> {
    if points.len() < 3 {
        return None;
    }
    let n = points.len() as f64;
    let mx = points.iter().map(|p| p.x).sum::<f64>() / n;
    let my = points.iter().map(|p| p.y).sum::<f64>() / n;
    let mz = points.iter().map(|p| p.z).sum::<f64>() / n;
    let mut cov = [[0.0_f64; 3]; 3];
    for p in points {
        let d = [p.x - mx, p.y - my, p.z - mz];
        for i in 0..3 {
            for j in 0..3 {
                cov[i][j] += d[i] * d[j];
            }
        }
    }
    for row in &mut cov {
        for val in row.iter_mut() {
            *val /= n;
        }
    }
    smallest_eigenvector_3x3(cov)
}

/// Estimate axis direction from cross products of normal pairs.
pub(super) fn estimate_axis(normals: &[Vec3]) -> Option<Vec3> {
    let unit_normals: Vec<Vec3> = normals.iter().filter_map(|n| n.normalized()).collect();
    if unit_normals.len() < 2 {
        return None;
    }
    let mut accum = Vec3::zero();
    let mut count = 0usize;
    for i in 0..unit_normals.len() {
        for j in (i + 1)..unit_normals.len() {
            let c = unit_normals[i].cross(&unit_normals[j]);
            if c.length_squared() > 1e-12 {
                let to_add = if count > 0 && c.dot(&accum) < 0.0 {
                    c * -1.0
                } else {
                    c
                };
                accum = accum + to_add;
                count += 1;
            }
        }
    }
    accum.normalized()
}

/// Estimate cone axis from normal differences (for analytical surfaces).
pub(super) fn estimate_cone_axis_from_normals(normals: &[Vec3]) -> Option<Vec3> {
    let mut un: Vec<Vec3> = Vec::new();
    for n in normals {
        if let Some(unit) = n.normalized()
            && !un
                .iter()
                .any(|prev| (unit - *prev).length_squared() < 1e-10)
        {
            un.push(unit);
        }
    }
    if un.len() < 3 {
        return None;
    }
    let mut accum = Vec3::zero();
    let mut count = 0usize;
    for i in 1..un.len() {
        let d1 = un[i] - un[0];
        if d1.length_squared() < 1e-12 {
            continue;
        }
        for j in (i + 1)..un.len() {
            let d2 = un[j] - un[0];
            let c = d1.cross(&d2);
            if c.length_squared() < 1e-12 {
                continue;
            }
            let to_add = if count > 0 && c.dot(&accum) < 0.0 {
                c * -1.0
            } else {
                c
            };
            accum = accum + to_add;
            count += 1;
        }
    }
    accum.normalized()
}

/// Fit a circle to 2D points using the algebraic Kasa method.
///
/// Returns `Some((cx, cy, radius))` or `None` if the system is degenerate.
pub(super) fn fit_circle_2d(pts: &[(f64, f64)]) -> Option<(f64, f64, f64)> {
    let n = pts.len() as f64;
    if n < 3.0 {
        return None;
    }
    let (mut sx, mut sy) = (0.0, 0.0);
    let (mut sxx, mut syy, mut sxy) = (0.0, 0.0, 0.0);
    let (mut sxxx, mut syyy) = (0.0, 0.0);
    let (mut sxxy, mut sxyy) = (0.0, 0.0);
    for &(x, y) in pts {
        sx += x;
        sy += y;
        sxx += x * x;
        syy += y * y;
        sxy += x * y;
        sxxx += x * x * x;
        syyy += y * y * y;
        sxxy += x * x * y;
        sxyy += x * y * y;
    }
    let rhs0 = sxxx + sxyy;
    let rhs1 = sxxy + syyy;
    let rhs2 = sxx + syy;
    let det = sxx * (syy * n - sy * sy) - sxy * (sxy * n - sy * sx) + sx * (sxy * sy - syy * sx);
    if det.abs() < 1e-30 {
        return None;
    }
    let a = (rhs0 * (syy * n - sy * sy) - sxy * (rhs1 * n - sy * rhs2)
        + sx * (rhs1 * sy - syy * rhs2))
        / det;
    let b = (sxx * (rhs1 * n - sy * rhs2) - rhs0 * (sxy * n - sy * sx)
        + sx * (sxy * rhs2 - rhs1 * sx))
        / det;
    let c = (sxx * (syy * rhs2 - rhs1 * sy) - sxy * (sxy * rhs2 - rhs1 * sx)
        + rhs0 * (sxy * sy - syy * sx))
        / det;
    let cx = a / 2.0;
    let cy = b / 2.0;
    let r_sq = c + cx * cx + cy * cy;
    if r_sq <= 0.0 {
        return None;
    }
    Some((cx, cy, r_sq.sqrt()))
}

/// Estimate center from closest approach of point+normal ray pairs.
pub(super) fn estimate_center_from_rays(
    points: &[Point3],
    normals: &[Vec3],
) -> Option<(Point3, Vec<f64>)> {
    let pairs: Vec<(Point3, Vec3)> = points
        .iter()
        .zip(normals.iter())
        .filter_map(|(p, n)| n.normalized().map(|un| (*p, un)))
        .collect();
    if pairs.len() < 3 {
        return None;
    }
    let mut accum = Vec3::zero();
    let mut cnt = 0usize;
    'outer: for i in 0..pairs.len() {
        for j in (i + 1)..pairs.len() {
            if cnt >= 50 {
                break 'outer;
            }
            let (pi, ni) = &pairs[i];
            let (pj, nj) = &pairs[j];
            let d = *pi - *pj;
            let a = ni.dot(ni);
            let b = ni.dot(nj);
            let c = nj.dot(nj);
            let denom = a * c - b * b;
            if denom.abs() < 1e-12 {
                continue;
            }
            let e = ni.dot(&d);
            let f = nj.dot(&d);
            let ci = *pi + *ni * ((b * f - c * e) / denom);
            let cj = *pj + *nj * ((a * f - b * e) / denom);
            accum = accum
                + Vec3::new(
                    (ci.x + cj.x) * 0.5,
                    (ci.y + cj.y) * 0.5,
                    (ci.z + cj.z) * 0.5,
                );
            cnt += 1;
        }
    }
    if cnt == 0 {
        return None;
    }
    let inv = 1.0 / cnt as f64;
    let center = Point3::new(accum.x * inv, accum.y * inv, accum.z * inv);
    let distances: Vec<f64> = points.iter().map(|p| (*p - center).length()).collect();
    Some((center, distances))
}
