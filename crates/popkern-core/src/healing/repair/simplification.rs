//! Surface simplification: detect and replace over-parameterized NURBS
//! surfaces with exact analytical forms.
//!
//! Currently detects planar, cylindrical, and spherical NURBS surfaces.
//! Cone and torus detection are planned for future iterations.

use serde::{Deserialize, Serialize};

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::topology::{BRep, FaceHandle};

/// Kind of analytical surface detected during simplification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnalyticalKind {
    /// An infinite plane.
    Plane,
    /// A cylinder.
    Cylinder,
    /// A cone.
    Cone,
    /// A sphere.
    Sphere,
    /// A torus.
    Torus,
}

/// Result of simplifying a single face's surface.
#[derive(Debug, Clone)]
pub struct SimplificationResult {
    /// The face whose surface was replaced.
    pub face: FaceHandle,
    /// The analytical form that was detected.
    pub detected_kind: AnalyticalKind,
    /// Maximum deviation between the original and fitted surface.
    pub max_deviation: f64,
}

/// Attempt to detect an analytical surface form from a parametric surface.
///
/// Samples the surface on a 5x5 grid within `parameter_range` and tests
/// for planarity, cylindricality, and sphericality.
///
/// Returns `Some((kind, geometry, max_deviation))` if a match is found,
/// or `None` if the surface does not fit any known analytical form.
pub fn detect_analytical_form(
    surface: &dyn SurfaceEval,
    parameter_range: ((f64, f64), (f64, f64)),
    tolerance: f64,
) -> Option<(AnalyticalKind, SurfaceGeometry, f64)> {
    let ((u0, u1), (v0, v1)) = parameter_range;
    let n = 5;

    // Sample points and normals on a 5x5 grid.
    let mut points = Vec::with_capacity(n * n);
    let mut normals = Vec::with_capacity(n * n);

    for i in 0..n {
        let u = u0 + (u1 - u0) * (i as f64) / ((n - 1) as f64);
        for j in 0..n {
            let v = v0 + (v1 - v0) * (j as f64) / ((n - 1) as f64);
            points.push(surface.point_at(u, v));
            normals.push(surface.normal_at(u, v));
        }
    }

    // Try planarity detection.
    if let Some(result) = try_detect_plane(&points, &normals, tolerance) {
        return Some(result);
    }

    // Try cylinder detection.
    if let Some((geom, dev)) = try_detect_cylinder(&points, &normals, tolerance) {
        return Some((AnalyticalKind::Cylinder, geom, dev));
    }

    // Try sphere detection.
    if let Some((geom, dev)) = try_detect_sphere(&points, &normals, tolerance) {
        return Some((AnalyticalKind::Sphere, geom, dev));
    }

    None
}

/// Check if all sampled normals are parallel and all points lie on a plane.
fn try_detect_plane(
    points: &[Point3],
    normals: &[Vec3],
    tolerance: f64,
) -> Option<(AnalyticalKind, SurfaceGeometry, f64)> {
    if points.is_empty() || normals.is_empty() {
        return None;
    }

    // Find first non-degenerate normal as reference.
    let ref_normal = normals.iter().find_map(|n| n.normalized())?;

    // Check all normals are parallel (angle < 1e-6 rad).
    let max_angle_tol = 1e-6;
    for n in normals {
        if let Some(unit_n) = n.normalized() {
            let angle = ref_normal.angle_to(&unit_n);
            // Normals could be parallel or anti-parallel.
            let deviation = angle.min(std::f64::consts::PI - angle);
            if deviation > max_angle_tol {
                return None;
            }
        }
    }

    // Fit a plane: origin = centroid, normal = reference normal.
    let count = points.len() as f64;
    let centroid = Point3::new(
        points.iter().map(|p| p.x).sum::<f64>() / count,
        points.iter().map(|p| p.y).sum::<f64>() / count,
        points.iter().map(|p| p.z).sum::<f64>() / count,
    );

    // Compute max deviation from the fitted plane.
    let mut max_dev = 0.0_f64;
    for p in points {
        let diff = *p - centroid;
        let dist = diff.dot(&ref_normal).abs();
        max_dev = max_dev.max(dist);
    }

    if max_dev > tolerance {
        return None;
    }

    // Build orthonormal (u_axis, v_axis) for the plane.
    let (u_axis, v_axis) = build_plane_axes(&ref_normal);

    let geometry = SurfaceGeometry::Plane {
        origin: centroid,
        normal: ref_normal,
        u_axis,
        v_axis,
    };

    Some((AnalyticalKind::Plane, geometry, max_dev))
}

/// Build two orthonormal axes perpendicular to a given normal.
fn build_plane_axes(normal: &Vec3) -> (Vec3, Vec3) {
    // Choose a reference that is least parallel to the normal.
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
///
/// Estimates the cylinder axis from cross products of normal pairs, projects
/// points onto the perpendicular plane, fits a circle using the algebraic
/// Kasa method, and checks that all points lie within `tolerance`.
fn try_detect_cylinder(
    points: &[Point3],
    normals: &[Vec3],
    tolerance: f64,
) -> Option<(SurfaceGeometry, f64)> {
    if points.len() < 3 || normals.len() < 3 {
        return None;
    }

    // Step 1: estimate axis from cross products of normal pairs.
    let unit_normals: Vec<Vec3> = normals.iter().filter_map(|n| n.normalized()).collect();
    if unit_normals.len() < 2 {
        return None;
    }

    let mut axis_accum = Vec3::zero();
    let mut cross_count = 0usize;
    for i in 0..unit_normals.len() {
        for j in (i + 1)..unit_normals.len() {
            let c = unit_normals[i].cross(&unit_normals[j]);
            if c.length_squared() > 1e-12 {
                let to_add = if cross_count > 0 && c.dot(&axis_accum) < 0.0 {
                    c * -1.0
                } else {
                    c
                };
                axis_accum = Vec3::new(
                    axis_accum.x + to_add.x,
                    axis_accum.y + to_add.y,
                    axis_accum.z + to_add.z,
                );
                cross_count += 1;
            }
        }
    }

    let axis = axis_accum.normalized()?;

    // Step 2: project all points onto a plane perpendicular to the axis.
    let (perp1, perp2) = orthonormal_basis(axis);

    // Use absolute coordinates for projection (not relative to centroid).
    let projected_2d: Vec<(f64, f64)> = points
        .iter()
        .map(|p| {
            let v = Vec3::new(p.x, p.y, p.z);
            (v.dot(&perp1), v.dot(&perp2))
        })
        .collect();

    // Step 3: algebraic circle fit (Kasa method).
    // Solve for (a, b, c) in: x^2 + y^2 = a*x + b*y + c
    // via least squares: [sum_xx, sum_xy, sum_x] [a]   [sum_x3_xy2]
    //                    [sum_xy, sum_yy, sum_y] [b] = [sum_x2y_y3]
    //                    [sum_x,  sum_y,  n    ] [c]   [sum_x2_y2 ]
    let (cx, cy, radius) = fit_circle_2d(&projected_2d)?;

    if radius < tolerance {
        return None;
    }

    // Step 4: compute 3D origin on the axis line.
    let origin_3d = Point3::new(0.0, 0.0, 0.0) + perp1 * cx + perp2 * cy;

    // Step 5: max deviation = |distance_from_axis - radius| for each point.
    let mut max_dev = 0.0_f64;
    for p in points {
        let v = *p - origin_3d;
        let along = v.dot(&axis);
        let radial = Vec3::new(
            v.x - axis.x * along,
            v.y - axis.y * along,
            v.z - axis.z * along,
        );
        max_dev = max_dev.max((radial.length() - radius).abs());
    }

    if max_dev > tolerance {
        return None;
    }

    let geom = SurfaceGeometry::Cylinder {
        origin: origin_3d,
        axis,
        radius,
    };
    Some((geom, max_dev))
}

/// Fit a circle to 2D points using the algebraic Kasa method.
///
/// Returns `Some((cx, cy, radius))` or `None` if the system is degenerate.
fn fit_circle_2d(pts: &[(f64, f64)]) -> Option<(f64, f64, f64)> {
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

    // Right-hand side.
    let rhs0 = sxxx + sxyy;
    let rhs1 = sxxy + syyy;
    let rhs2 = sxx + syy;

    // Solve 3x3 system via Cramer's rule.
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

/// Attempt to detect a spherical surface from sampled points and normals.
///
/// Estimates the sphere center by finding the closest approach of
/// point+normal ray pairs, then checks that all points lie within
/// `tolerance` of the fitted sphere.
fn try_detect_sphere(
    points: &[Point3],
    normals: &[Vec3],
    tolerance: f64,
) -> Option<(SurfaceGeometry, f64)> {
    if points.len() < 3 || normals.len() < 3 {
        return None;
    }

    // Collect valid (point, unit_normal) pairs.
    let pairs: Vec<(Point3, Vec3)> = points
        .iter()
        .zip(normals.iter())
        .filter_map(|(p, n)| n.normalized().map(|un| (*p, un)))
        .collect();

    if pairs.len() < 3 {
        return None;
    }

    // Estimate center from closest approach of ray pairs.
    let mut center_accum = Vec3::zero();
    let mut pair_count = 0usize;
    let max_pairs = 50;

    'outer: for i in 0..pairs.len() {
        for j in (i + 1)..pairs.len() {
            if pair_count >= max_pairs {
                break 'outer;
            }
            let (pi, ni) = &pairs[i];
            let (pj, nj) = &pairs[j];

            let d = *pi - *pj;
            let a = ni.dot(ni);
            let b = ni.dot(nj);
            let c = nj.dot(nj);
            let e = ni.dot(&d);
            let f = nj.dot(&d);

            let denom = a * c - b * b;
            if denom.abs() < 1e-12 {
                continue; // Lines are parallel, skip.
            }

            let t = (b * f - c * e) / denom;
            let s = (a * f - b * e) / denom;

            let close_i = *pi + *ni * t;
            let close_j = *pj + *nj * s;
            let mid = Vec3::new(
                (close_i.x + close_j.x) * 0.5,
                (close_i.y + close_j.y) * 0.5,
                (close_i.z + close_j.z) * 0.5,
            );

            center_accum = Vec3::new(
                center_accum.x + mid.x,
                center_accum.y + mid.y,
                center_accum.z + mid.z,
            );
            pair_count += 1;
        }
    }

    if pair_count == 0 {
        return None;
    }

    let center = Point3::new(
        center_accum.x / pair_count as f64,
        center_accum.y / pair_count as f64,
        center_accum.z / pair_count as f64,
    );

    // Compute radius = average distance from center to all points.
    let count = points.len() as f64;
    let radius: f64 = points.iter().map(|p| (*p - center).length()).sum::<f64>() / count;

    if radius < tolerance {
        return None; // Degenerate.
    }

    // Max deviation = max |distance_to_center - radius|.
    let mut max_dev = 0.0_f64;
    for p in points {
        let dist = (*p - center).length();
        max_dev = max_dev.max((dist - radius).abs());
    }

    if max_dev > tolerance {
        return None;
    }

    let geom = SurfaceGeometry::Sphere { center, radius };
    Some((geom, max_dev))
}

/// Simplify a single face's surface if it is an over-parameterized NURBS.
///
/// If the face's surface is already analytical (Plane, Cylinder, etc.),
/// returns `None`.  If it is a NURBS surface that matches an analytical
/// form within `tolerance`, the face's surface is replaced in place and
/// a [`SimplificationResult`] is returned.
pub fn simplify_face(
    brep: &mut BRep,
    face: FaceHandle,
    tolerance: f64,
) -> Option<SimplificationResult> {
    let face_data = brep.faces.get(face)?;
    let sof = face_data.surface.as_ref()?;

    // Only attempt simplification on NURBS surfaces.
    let is_nurbs = matches!(sof.geometry, SurfaceGeometry::Nurbs(_));
    if !is_nurbs {
        return None;
    }

    let param_range = sof.geometry.parameter_range();
    let orientation = sof.orientation;

    let (kind, new_geom, max_dev) = detect_analytical_form(&sof.geometry, param_range, tolerance)?;

    // Replace the face's surface with the detected analytical form.
    if let Some(f) = brep.faces.get_mut(face) {
        f.surface = Some(SurfaceOnFace {
            geometry: new_geom,
            orientation,
        });
    }

    Some(SimplificationResult {
        face,
        detected_kind: kind,
        max_deviation: max_dev,
    })
}

/// Simplify all faces in the BRep whose NURBS surfaces match analytical forms.
///
/// Iterates every face and calls [`simplify_face`].  Returns a list of
/// all faces that were successfully simplified.
pub fn simplify_all_faces(brep: &mut BRep, tolerance: f64) -> Vec<SimplificationResult> {
    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    let mut results = Vec::new();
    for fh in face_handles {
        if let Some(result) = simplify_face(brep, fh, tolerance) {
            results.push(result);
        }
    }
    results
}

#[cfg(test)]
#[path = "simplification_tests.rs"]
mod tests;
