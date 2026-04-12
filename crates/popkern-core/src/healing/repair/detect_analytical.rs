//! Analytical surface detection from sampled points and normals.
//!
//! Detects planar, cylindrical, spherical, conical, and toroidal surfaces
//! by fitting geometric primitives to sampled surface data.

#[path = "detect_analytical_helpers.rs"]
mod helpers;

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};
use helpers::{
    estimate_axis, estimate_axis_from_positions, estimate_center_from_rays,
    estimate_cone_axis_from_normals, fit_circle_2d,
};

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
    // Try multiple axis estimation methods and keep the first that yields a valid cone.
    let axis_candidates = [
        estimate_cone_axis_from_normals(normals),
        estimate_axis(normals),
        estimate_axis_from_positions(points),
    ];
    for axis_opt in &axis_candidates {
        let Some(axis) = axis_opt else { continue };
        if let Some(result) = try_fit_cone_with_axis(points, normals, *axis, tolerance) {
            return Some(result);
        }
    }
    None
}

/// Try fitting a cone with a given axis direction.
fn try_fit_cone_with_axis(
    points: &[Point3],
    normals: &[Vec3],
    axis: Vec3,
    tolerance: f64,
) -> Option<(Point3, Vec3, f64, f64)> {
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
    // Geometry-based half-angle from the linear fit slope, with normal-based fallback.
    let geo_ha = slope.abs().atan();
    let un: Vec<Vec3> = normals.iter().filter_map(|nn| nn.normalized()).collect();
    let mean_dot = un.iter().map(|nn| nn.dot(&axis).abs()).sum::<f64>() / (un.len().max(1) as f64);
    let normal_ha = mean_dot.clamp(0.0, 1.0).asin();
    let half_angle = if (1e-6..=std::f64::consts::FRAC_PI_2 - 1e-6).contains(&geo_ha) {
        geo_ha
    } else if (1e-6..=std::f64::consts::FRAC_PI_2 - 1e-6).contains(&normal_ha) {
        normal_ha
    } else {
        return None;
    };
    let tan_ha = half_angle.tan();
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

/// Compute variance of distances-to-major-circle for a given center candidate.
fn torus_fit_variance(points: &[Point3], center: Point3, axis: Vec3, major_r: f64) -> f64 {
    let dists: Vec<f64> = points
        .iter()
        .map(|p| dist_to_major_circle(*p - center, axis, major_r))
        .collect();
    let mean = dists.iter().sum::<f64>() / dists.len() as f64;
    dists.iter().map(|d| (d - mean) * (d - mean)).sum::<f64>() / dists.len() as f64
}

/// Estimate the torus center and major radius from per-point radial projections.
///
/// Projects points onto the plane perpendicular to the axis and iteratively
/// refines the center estimate. Then optimizes the height along the axis
/// to minimize variance in minor-radius distances.
fn estimate_torus_center(
    points: &[Point3],
    axis: Vec3,
    perp1: Vec3,
    perp2: Vec3,
) -> Option<(Point3, f64)> {
    let n = points.len() as f64;
    // Project all points to (perp1, perp2) coordinates.
    let pts_2d: Vec<(f64, f64)> = points
        .iter()
        .map(|p| {
            let pv = Vec3::new(p.x, p.y, p.z);
            (pv.dot(&perp1), pv.dot(&perp2))
        })
        .collect();
    // Initial center: centroid of 2D projections.
    let mut cx = pts_2d.iter().map(|(x, _)| x).sum::<f64>() / n;
    let mut cy = pts_2d.iter().map(|(_, y)| y).sum::<f64>() / n;
    let mut major_r = 0.0;
    // Iteratively refine in-plane center: each point votes for center at P - dir*R.
    for _ in 0..50 {
        let mut sum_r = 0.0_f64;
        let mut new_cx = 0.0_f64;
        let mut new_cy = 0.0_f64;
        let mut cnt = 0usize;
        for &(px, py) in &pts_2d {
            let dx = px - cx;
            let dy = py - cy;
            let r = (dx * dx + dy * dy).sqrt();
            if r > 1e-15 {
                sum_r += r;
                cnt += 1;
            }
        }
        if cnt == 0 {
            return None;
        }
        major_r = sum_r / cnt as f64;
        for &(px, py) in &pts_2d {
            let dx = px - cx;
            let dy = py - cy;
            let r = (dx * dx + dy * dy).sqrt();
            if r > 1e-15 {
                new_cx += px - dx / r * major_r;
                new_cy += py - dy / r * major_r;
            }
        }
        cx = new_cx / cnt as f64;
        cy = new_cy / cnt as f64;
    }
    // Refine height along axis to minimize variance of dist_to_major_circle.
    let heights: Vec<f64> = points
        .iter()
        .map(|p| Vec3::new(p.x, p.y, p.z).dot(&axis))
        .collect();
    let h_min = heights.iter().copied().fold(f64::INFINITY, f64::min);
    let h_max = heights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let base_center = Point3::new(0.0, 0.0, 0.0) + perp1 * cx + perp2 * cy;
    // Golden-section search for optimal height.
    let gr = 0.5 * (5.0_f64.sqrt() - 1.0);
    let (mut a, mut b) = (h_min - (h_max - h_min), h_max + (h_max - h_min));
    let mut c = b - gr * (b - a);
    let mut d = a + gr * (b - a);
    for _ in 0..30 {
        let fc = torus_fit_variance(points, base_center + axis * c, axis, major_r);
        let fd = torus_fit_variance(points, base_center + axis * d, axis, major_r);
        if fc < fd {
            b = d;
        } else {
            a = c;
        }
        c = b - gr * (b - a);
        d = a + gr * (b - a);
    }
    let ch = (a + b) * 0.5;
    let center = base_center + axis * ch;
    Some((center, major_r))
}

/// Try to detect a torus from sampled surface points and normals.
///
/// Estimates axis via position covariance (falling back to normal
/// cross-products), iteratively fits the major circle center and
/// radius, then derives minor_radius from distances to that circle.
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
    let axis = estimate_axis_from_positions(points)
        .or_else(|| estimate_cone_axis_from_normals(normals))
        .or_else(|| estimate_axis(normals))?;
    let (perp1, perp2) = orthonormal_basis(axis);
    let (center, major_radius) = estimate_torus_center(points, axis, perp1, perp2)?;
    if major_radius < tolerance {
        return None;
    }
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
