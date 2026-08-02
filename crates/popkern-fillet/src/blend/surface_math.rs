//! Math utilities for blend surface construction.
//!
//! Point cloud centroid, least-squares plane fitting, circumcenter
//! computation, and eigenvector extraction.

use popkern_core::geometry::{Point3, Vec3};

/// Solve a 3x3 linear system via Gaussian elimination with partial
/// pivoting.
///
/// Takes 3 rows of [a, b, c, rhs] and returns (x, y, z) or None if
/// singular (pivot below 1e-15 — callers should scale their rows to
/// O(1) so this acts as a relative rank test).
#[allow(clippy::needless_range_loop)]
pub(super) fn solve_linear_3x3(rows: &[[f64; 4]]) -> Option<(f64, f64, f64)> {
    if rows.len() < 3 {
        return None;
    }
    let mut m = [[0.0f64; 4]; 3];
    m[0] = rows[0];
    m[1] = rows[1];
    m[2] = rows[2];

    for col in 0..3 {
        let mut max_row = col;
        let mut max_val = m[col][col].abs();
        for row in (col + 1)..3 {
            if m[row][col].abs() > max_val {
                max_val = m[row][col].abs();
                max_row = row;
            }
        }
        if max_val < 1e-15 {
            return None;
        }
        m.swap(col, max_row);

        for row in (col + 1)..3 {
            let factor = m[row][col] / m[col][col];
            for c in col..4 {
                m[row][c] -= factor * m[col][c];
            }
        }
    }

    let cz = m[2][3] / m[2][2];
    let cy = (m[1][3] - m[1][2] * cz) / m[1][1];
    let cx = (m[0][3] - m[0][1] * cy - m[0][2] * cz) / m[0][0];

    if cx.is_finite() && cy.is_finite() && cz.is_finite() {
        Some((cx, cy, cz))
    } else {
        None
    }
}

/// Least-squares sphere fit over ALL samples (PK-233).
///
/// Solves the equidistance system `|p_i - c|² = r²` linearly: with
/// centroid-shifted coordinates `q_i = p_i - p̄` the per-point
/// equation `2·q_i·c' + (r² - |c'|²) = |q_i|²` has its affine
/// unknown eliminated by subtracting the mean equation (Σq_i = 0),
/// leaving the 3×3 normal equations
/// `(Σ q_i q_iᵀ)·c' = ½·Σ q_i·|q_i|²`.
///
/// The scatter matrix is singular exactly when the samples are
/// coplanar — i.e. when they genuinely cannot determine a unique
/// sphere (co-circular data lies on infinitely many spheres). This
/// replaces the first-3-differences circumcenter rows, which were
/// adjacent samples of ONE boundary arc and therefore co-circular —
/// singular even when the full sample set determined the sphere
/// exactly (PK-229 probe: regular pyramid apexes with concurrent
/// cylinder axes never promoted).
///
/// Coordinates are scale-normalized before the solve so the pivot
/// threshold is meaningful at any corner radius. Returns
/// `(center, radius)` — the caller certifies with its own deviation
/// gate; this function imposes no acceptance window.
pub(super) fn fit_sphere_least_squares(points: &[Point3]) -> Option<(Point3, f64)> {
    if points.len() < 4 {
        return None;
    }

    let centroid = compute_centroid(points);
    let scale = points
        .iter()
        .map(|p| p.distance_to(&centroid))
        .fold(0.0_f64, f64::max);
    if scale <= 0.0 || !scale.is_finite() {
        return None;
    }

    // Normal equations over unit-scale shifted coordinates,
    // averaged so entries are O(1) and the solver's pivot floor acts
    // as a relative rank test.
    let inv_n = 1.0 / points.len() as f64;
    let mut m = [[0.0f64; 3]; 3];
    let mut rhs = [0.0f64; 3];
    for p in points {
        let q = [
            (p.x - centroid.x) / scale,
            (p.y - centroid.y) / scale,
            (p.z - centroid.z) / scale,
        ];
        let q_sq = q[0] * q[0] + q[1] * q[1] + q[2] * q[2];
        for r in 0..3 {
            for c in 0..3 {
                m[r][c] += q[r] * q[c] * inv_n;
            }
            rhs[r] += 0.5 * q[r] * q_sq * inv_n;
        }
    }

    let rows = [
        [m[0][0], m[0][1], m[0][2], rhs[0]],
        [m[1][0], m[1][1], m[1][2], rhs[1]],
        [m[2][0], m[2][1], m[2][2], rhs[2]],
    ];
    let (cx, cy, cz) = solve_linear_3x3(&rows)?;

    let center = Point3::new(
        centroid.x + cx * scale,
        centroid.y + cy * scale,
        centroid.z + cz * scale,
    );
    let radius = points.iter().map(|p| p.distance_to(&center)).sum::<f64>() / points.len() as f64;
    if !radius.is_finite() {
        return None;
    }
    Some((center, radius))
}

/// Compute the centroid of a point cloud.
pub(super) fn compute_centroid(points: &[Point3]) -> Point3 {
    let n = points.len() as f64;
    let x = points.iter().map(|p| p.x).sum::<f64>() / n;
    let y = points.iter().map(|p| p.y).sum::<f64>() / n;
    let z = points.iter().map(|p| p.z).sum::<f64>() / n;
    Point3::new(x, y, z)
}

/// Fit a least-squares plane through a set of points.
///
/// Returns `(centroid, normal)`. The normal is the eigenvector corresponding
/// to the smallest eigenvalue of the covariance matrix.
pub(super) fn fit_plane_least_squares(points: &[Point3]) -> Result<(Point3, Vec3), String> {
    if points.len() < 3 {
        return Err("need at least 3 points to fit a plane".into());
    }

    let centroid = compute_centroid(points);

    let mut cov = [[0.0f64; 3]; 3];
    for p in points {
        let dx = p.x - centroid.x;
        let dy = p.y - centroid.y;
        let dz = p.z - centroid.z;
        let d = [dx, dy, dz];
        for r in 0..3 {
            for c in 0..3 {
                cov[r][c] += d[r] * d[c];
            }
        }
    }

    let normal = smallest_eigenvector_3x3(&cov)?;
    Ok((centroid, normal))
}

/// Find the eigenvector corresponding to the smallest eigenvalue of a
/// symmetric 3x3 matrix via power iteration + deflation.
#[allow(clippy::needless_range_loop)]
fn smallest_eigenvector_3x3(cov: &[[f64; 3]; 3]) -> Result<Vec3, String> {
    let mut v = [1.0, 0.0, 0.0];
    for _ in 0..50 {
        let mut w = [0.0f64; 3];
        for r in 0..3 {
            for c in 0..3 {
                w[r] += cov[r][c] * v[c];
            }
        }
        let len = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
        if len < 1e-30 {
            break;
        }
        v = [w[0] / len, w[1] / len, w[2] / len];
    }
    let e1 = Vec3::new(v[0], v[1], v[2]);

    let lambda1: f64 = {
        let mut w = [0.0f64; 3];
        for r in 0..3 {
            for c in 0..3 {
                w[r] += cov[r][c] * v[c];
            }
        }
        (w[0] * v[0] + w[1] * v[1] + w[2] * v[2]).max(0.0)
    };

    let mut cov2 = *cov;
    for r in 0..3 {
        for c in 0..3 {
            cov2[r][c] -= lambda1 * v[r] * v[c];
        }
    }

    let mut v2 = [0.0, 1.0, 0.0];
    for _ in 0..50 {
        let mut w = [0.0f64; 3];
        for r in 0..3 {
            for c in 0..3 {
                w[r] += cov2[r][c] * v2[c];
            }
        }
        let len = (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt();
        if len < 1e-30 {
            break;
        }
        v2 = [w[0] / len, w[1] / len, w[2] / len];
    }
    let e2 = Vec3::new(v2[0], v2[1], v2[2]);

    let normal = e1.cross(&e2);
    normal
        .normalized()
        .ok_or_else(|| "degenerate covariance — cannot fit plane".to_string())
}

/// Pick a vector perpendicular to the given vector.
pub(super) fn pick_perpendicular(v: &Vec3) -> Vec3 {
    let candidate = if v.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let perp = v.cross(&candidate);
    perp.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Latitude-circle samples of the sphere `center + radius·ŝ(φ, θ)`.
    /// All samples of one call are co-circular (constant polar angle φ).
    fn latitude_arc(
        center: Point3,
        radius: f64,
        phi: f64,
        theta_range: (f64, f64),
        n: usize,
    ) -> Vec<Point3> {
        (0..n)
            .map(|i| {
                let t = i as f64 / (n - 1) as f64;
                let theta = theta_range.0 + t * (theta_range.1 - theta_range.0);
                Point3::new(
                    center.x + radius * phi.sin() * theta.cos(),
                    center.y + radius * phi.sin() * theta.sin(),
                    center.z + radius * phi.cos(),
                )
            })
            .collect()
    }

    /// PK-233 acceptance: the first boundary arc's samples are
    /// co-circular, which made the predecessor's 3×3 (built from the
    /// first 3 point differences) singular. The least-squares fit
    /// over all samples must recover the sphere exactly.
    #[test]
    fn co_circular_first_arc_samples_recover_known_sphere() {
        let center = Point3::new(1.25, -2.5, 3.75);
        let radius = 0.35;
        let mut points = latitude_arc(center, radius, 0.9, (0.0, 1.5), 8);
        points.extend(latitude_arc(center, radius, 1.4, (0.3, 1.8), 8));
        points.extend(latitude_arc(center, radius, 0.4, (0.6, 2.1), 8));

        let (c, r) =
            fit_sphere_least_squares(&points).expect("three distinct arcs determine the sphere");
        assert!(
            c.distance_to(&center) < 1e-9,
            "center off by {}",
            c.distance_to(&center)
        );
        assert!(
            (r - radius).abs() < 1e-9,
            "radius off by {}",
            (r - radius).abs()
        );
    }

    /// A single circle lies on infinitely many spheres: the scatter
    /// matrix is rank-2 and the fit must decline, not guess.
    #[test]
    fn single_circle_samples_decline_underdetermined() {
        let points = latitude_arc(Point3::new(0.0, 0.0, 0.0), 1.0, 1.1, (0.0, 6.0), 16);
        assert!(fit_sphere_least_squares(&points).is_none());
    }

    /// Coplanar non-concyclic samples lie on no sphere and leave the
    /// center's normal component undetermined — must decline.
    #[test]
    fn coplanar_grid_samples_decline() {
        let mut points = Vec::new();
        for i in 0..4 {
            for j in 0..4 {
                points.push(Point3::new(i as f64 * 0.3, j as f64 * 0.2, 0.7));
            }
        }
        assert!(fit_sphere_least_squares(&points).is_none());
    }

    /// Fit noise tolerance: samples with ~1e-6 displacement (NURBS
    /// boundary fit drift scale) recover the center to the same
    /// order — the caller's deviation gate does the acceptance.
    #[test]
    fn noisy_arc_samples_recover_center_to_noise_order() {
        let center = Point3::new(-0.5, 0.25, 1.0);
        let radius = 0.08;
        let mut points = latitude_arc(center, radius, 0.8, (0.0, 1.2), 8);
        points.extend(latitude_arc(center, radius, 1.5, (0.4, 1.6), 8));
        points.extend(latitude_arc(center, radius, 0.35, (0.8, 2.0), 8));
        for (i, p) in points.iter_mut().enumerate() {
            let s = if i % 2 == 0 { 1.0 } else { -1.0 };
            p.x += s * 1.0e-6;
            p.z -= s * 0.7e-6;
        }
        let (c, r) =
            fit_sphere_least_squares(&points).expect("noisy arcs still determine the sphere");
        assert!(
            c.distance_to(&center) < 1e-4,
            "center off by {}",
            c.distance_to(&center)
        );
        assert!((r - radius).abs() < 1e-4);
    }
}
