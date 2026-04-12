//! Analytical surface detection from sampled points and normals.
//!
//! Detects planar, cylindrical, spherical, conical, and toroidal surfaces
//! by fitting geometric primitives to sampled surface data.

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};

/// Check if all sampled normals are parallel and all points lie on a plane.
pub(super) fn try_detect_plane(
    points: &[Point3],
    normals: &[Vec3],
    tolerance: f64,
) -> Option<(Vec3, Point3, f64)> {
    if points.is_empty() || normals.is_empty() {
        return None;
    }
    let ref_normal = normals.iter().find_map(|n| n.normalized())?;
    for n in normals {
        if let Some(unit_n) = n.normalized() {
            let angle = ref_normal.angle_to(&unit_n);
            if angle.min(std::f64::consts::PI - angle) > 1e-6 {
                return None;
            }
        }
    }
    let count = points.len() as f64;
    let centroid = Point3::new(
        points.iter().map(|p| p.x).sum::<f64>() / count,
        points.iter().map(|p| p.y).sum::<f64>() / count,
        points.iter().map(|p| p.z).sum::<f64>() / count,
    );
    let mut max_dev = 0.0_f64;
    for p in points {
        max_dev = max_dev.max((*p - centroid).dot(&ref_normal).abs());
    }
    if max_dev > tolerance {
        return None;
    }
    Some((ref_normal, centroid, max_dev))
}

/// Build two orthonormal axes perpendicular to a given normal.
pub(super) fn build_plane_axes(normal: &Vec3) -> (Vec3, Vec3) {
    let reference = if normal.x.abs() <= normal.y.abs() && normal.x.abs() <= normal.z.abs() {
        Vec3::new(1.0, 0.0, 0.0)
    } else if normal.y.abs() <= normal.z.abs() {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let proj = reference - *normal * reference.dot(normal);
    let u_axis = proj.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v_axis = normal.cross(&u_axis);
    (u_axis, v_axis)
}

/// Estimate axis direction from cross products of normal pairs.
fn estimate_axis(normals: &[Vec3]) -> Option<Vec3> {
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

/// Attempt to detect a cylindrical surface from sampled points and normals.
pub(super) fn try_detect_cylinder(
    points: &[Point3],
    normals: &[Vec3],
    tolerance: f64,
) -> Option<(Vec3, Point3, f64, f64)> {
    if points.len() < 3 || normals.len() < 3 {
        return None;
    }
    let axis = estimate_axis(normals)?;
    let (perp1, perp2) = orthonormal_basis(axis);
    let proj: Vec<(f64, f64)> = points
        .iter()
        .map(|p| {
            let v = Vec3::new(p.x, p.y, p.z);
            (v.dot(&perp1), v.dot(&perp2))
        })
        .collect();
    let (cx, cy, radius) = fit_circle_2d(&proj)?;
    if radius < tolerance {
        return None;
    }
    let origin = Point3::new(0.0, 0.0, 0.0) + perp1 * cx + perp2 * cy;
    let mut max_dev = 0.0_f64;
    for p in points {
        let v = *p - origin;
        let along = v.dot(&axis);
        max_dev = max_dev.max(((v - axis * along).length() - radius).abs());
    }
    if max_dev > tolerance {
        return None;
    }
    Some((axis, origin, radius, max_dev))
}

/// Estimate center from closest approach of point+normal ray pairs.
fn estimate_center_from_rays(points: &[Point3], normals: &[Vec3]) -> Option<(Point3, Vec<f64>)> {
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

/// Attempt to detect a spherical surface from sampled points and normals.
pub(super) fn try_detect_sphere(
    points: &[Point3],
    normals: &[Vec3],
    tolerance: f64,
) -> Option<(Point3, f64, f64)> {
    if points.len() < 3 || normals.len() < 3 {
        return None;
    }
    let (center, distances) = estimate_center_from_rays(points, normals)?;
    let radius: f64 = distances.iter().sum::<f64>() / distances.len() as f64;
    if radius < tolerance {
        return None;
    }
    let mut max_dev = 0.0_f64;
    for &dist in &distances {
        max_dev = max_dev.max((dist - radius).abs());
    }
    if max_dev > tolerance {
        return None;
    }
    Some((center, radius, max_dev))
}

/// Try to detect a cone from sampled surface points and normals.
///
/// Estimates axis via cross-products of normal differences, derives
/// half-angle from normal-axis dot products, and locates the apex
/// via linear regression of radial distance vs height.
pub(super) fn try_detect_cone(
    surface: &dyn SurfaceEval,
    points: &[Point3],
    normals: &[Vec3],
    tolerance: f64,
) -> Option<(Point3, Vec3, f64, f64)> {
    let _ = surface;
    if points.len() < 4 || normals.len() < 4 {
        return None;
    }
    // Try both axis estimation methods and keep the one that yields a valid cone.
    let axis_candidates = [
        estimate_cone_axis_from_normals(normals),
        estimate_axis(normals),
    ];
    for axis_opt in &axis_candidates {
        let Some(axis) = axis_opt else { continue };
        if let Some(result) = try_fit_cone_with_axis(points, normals, *axis, tolerance) {
            return Some(result);
        }
    }
    None
}

/// Estimate cone axis from normal differences (for analytical surfaces).
fn estimate_cone_axis_from_normals(normals: &[Vec3]) -> Option<Vec3> {
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

/// Try fitting a cone with a given axis direction.
fn try_fit_cone_with_axis(
    points: &[Point3],
    normals: &[Vec3],
    axis: Vec3,
    tolerance: f64,
) -> Option<(Point3, Vec3, f64, f64)> {
    // Half-angle from normal-axis dot products.
    let un: Vec<Vec3> = normals.iter().filter_map(|n| n.normalized()).collect();
    let mean_dot = un.iter().map(|n| n.dot(&axis).abs()).sum::<f64>() / un.len() as f64;
    let half_angle = mean_dot.clamp(0.0, 1.0).asin();
    if !(1e-6..=std::f64::consts::FRAC_PI_2 - 1e-6).contains(&half_angle) {
        return None;
    }
    let tan_ha = half_angle.tan();
    let n = points.len() as f64;
    let (perp1, perp2) = orthonormal_basis(axis);
    // Linear regression: R = slope * h + intercept (per-point h, R from axis).
    let mut sh = 0.0_f64;
    let mut sr = 0.0_f64;
    let mut shh = 0.0_f64;
    let mut shr = 0.0_f64;
    for p in points {
        let pv = Vec3::new(p.x, p.y, p.z);
        let h = pv.dot(&axis);
        let r = (pv - axis * h).length();
        sh += h;
        sr += r;
        shh += h * h;
        shr += h * r;
    }
    let det = n * shh - sh * sh;
    if det.abs() < 1e-30 {
        return None;
    }
    let slope = (n * shr - sh * sr) / det;
    let intercept = (sr * shh - sh * shr) / det;
    if slope.abs() < 1e-6 {
        return None;
    }
    let h_apex = -intercept / slope;
    // Axis 2D position: extrapolate each point's radial direction to R=0.
    let (mut ax, mut ay) = (0.0_f64, 0.0_f64);
    let mut cnt = 0usize;
    for p in points {
        let pv = Vec3::new(p.x, p.y, p.z);
        let h = pv.dot(&axis);
        let px = pv.dot(&perp1);
        let py = pv.dot(&perp2);
        let r = (px * px + py * py).sqrt();
        if r > 1e-15 {
            let r_exp = slope * h + intercept;
            ax += px - px / r * r_exp;
            ay += py - py / r * r_exp;
            cnt += 1;
        }
    }
    if cnt == 0 {
        return None;
    }
    let apex = Point3::new(0.0, 0.0, 0.0)
        + perp1 * (ax / cnt as f64)
        + perp2 * (ay / cnt as f64)
        + axis * h_apex;
    // Verify all points lie on the cone surface.
    let mut max_dev = 0.0_f64;
    for p in points {
        let v = *p - apex;
        let along = v.dot(&axis);
        max_dev = max_dev.max(((v - axis * along).length() - along.abs() * tan_ha).abs());
    }
    if max_dev > tolerance {
        return None;
    }
    Some((apex, axis, half_angle, max_dev))
}

/// Compute distance from a point offset to the torus major circle.
fn dist_to_major_circle(v: Vec3, axis: Vec3, major_radius: f64) -> f64 {
    let along = v.dot(&axis);
    let in_plane = v - axis * along;
    let d = in_plane.length();
    if d < 1e-15 {
        return major_radius;
    }
    (v - in_plane * (major_radius / d)).length()
}

/// Try to detect a torus from sampled surface points and normals.
///
/// Estimates axis via normal cross-products, fits the major circle,
/// then derives minor_radius from distances to that circle.
pub(super) fn try_detect_torus(
    surface: &dyn SurfaceEval,
    points: &[Point3],
    normals: &[Vec3],
    tolerance: f64,
) -> Option<(Point3, Vec3, f64, f64, f64)> {
    let _ = surface;
    if points.len() < 6 || normals.len() < 6 {
        return None;
    }
    let axis = estimate_cone_axis_from_normals(normals).or_else(|| estimate_axis(normals))?;
    let n = points.len() as f64;
    let ctr = Point3::new(
        points.iter().map(|p| p.x).sum::<f64>() / n,
        points.iter().map(|p| p.y).sum::<f64>() / n,
        points.iter().map(|p| p.z).sum::<f64>() / n,
    );
    let (perp1, perp2) = orthonormal_basis(axis);
    let proj: Vec<(f64, f64)> = points
        .iter()
        .map(|p| {
            let v = *p - ctr;
            (v.dot(&perp1), v.dot(&perp2))
        })
        .collect();
    let (cx, cy, major_radius) = fit_circle_2d(&proj)?;
    if major_radius < tolerance {
        return None;
    }
    let center = ctr + perp1 * cx + perp2 * cy;
    // Compute minor_radius as average distance from each point to major circle.
    let dists: Vec<f64> = points
        .iter()
        .map(|p| dist_to_major_circle(*p - center, axis, major_radius))
        .collect();
    let minor_radius = dists.iter().sum::<f64>() / dists.len() as f64;
    if minor_radius < tolerance || minor_radius >= major_radius {
        return None;
    }
    let mut max_dev = 0.0_f64;
    for &d in &dists {
        max_dev = max_dev.max((d - minor_radius).abs());
    }
    if max_dev > tolerance {
        return None;
    }
    Some((center, axis, major_radius, minor_radius, max_dev))
}

#[cfg(test)]
#[path = "detect_analytical_tests.rs"]
mod tests;
