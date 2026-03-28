//! Math utilities for blend surface construction.
//!
//! Point cloud centroid, least-squares plane fitting, circumcenter
//! computation, and eigenvector extraction.

use popkern_core::geometry::{Point3, Vec3};

/// Solve 3x3 linear system for circumcenter via Gaussian elimination.
///
/// Takes 3 rows of [a, b, c, rhs] and returns (cx, cy, cz) or None if singular.
#[allow(clippy::needless_range_loop)]
pub(super) fn solve_circumcenter_3x3(rows: &[[f64; 4]]) -> Option<(f64, f64, f64)> {
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
