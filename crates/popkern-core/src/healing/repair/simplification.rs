//! Surface simplification: detect and replace over-parameterized NURBS
//! surfaces with exact analytical forms.
//!
//! Currently detects planar NURBS surfaces.  Cylinder, sphere, cone, and
//! torus detection are planned for future iterations.

use serde::{Deserialize, Serialize};

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
/// for planarity (all normals parallel within 1e-6 rad, all points within
/// `tolerance` of a fitted plane).
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

    // Future: try_detect_cylinder, try_detect_sphere, etc.
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
mod tests {
    use super::*;
    use crate::geometry::surface::SurfaceGeometry;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn detect_plane_from_planar_points() {
        // Create a plane surface and sample it.
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let result = detect_analytical_form(&plane, ((-1.0, 1.0), (-1.0, 1.0)), 1e-6);
        assert!(result.is_some(), "should detect plane");
        let (kind, _, dev) = result.unwrap();
        assert_eq!(kind, AnalyticalKind::Plane);
        assert!(dev < 1e-10, "deviation should be near zero for exact plane");
    }

    #[test]
    fn detect_fails_on_cylinder() {
        let cylinder = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        // A cylinder is not planar, so plane detection should fail.
        let result = detect_analytical_form(&cylinder, ((0.0, 1.0), (0.0, 1.0)), 1e-6);
        assert!(result.is_none(), "cylinder should not be detected as plane");
    }

    #[test]
    fn simplify_face_on_cube_is_noop() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        // Cube faces are already Plane geometry, so simplification should be a no-op.
        let result = simplify_face(&mut brep, face, 1e-6);
        assert!(
            result.is_none(),
            "already-analytical face should not be simplified"
        );
    }

    #[test]
    fn simplify_all_faces_on_cube_returns_empty() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let results = simplify_all_faces(&mut brep, 1e-6);
        assert!(results.is_empty(), "all cube faces are already planar");
    }

    #[test]
    fn analytical_kind_values() {
        assert_ne!(AnalyticalKind::Plane, AnalyticalKind::Cylinder);
        assert_ne!(AnalyticalKind::Sphere, AnalyticalKind::Torus);
    }
}
