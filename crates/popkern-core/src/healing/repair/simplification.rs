//! Surface simplification: detect and replace over-parameterized NURBS
//! surfaces with exact analytical forms.
//!
//! Detects planar, cylindrical, spherical, conical, and toroidal surfaces.
//! Detection logic lives in [`super::detect_analytical`]; this module provides
//! the high-level `simplify_face` / `simplify_all_faces` API.

use serde::{Deserialize, Serialize};

use crate::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use crate::topology::{BRep, FaceHandle};

use super::detect_analytical::{
    build_plane_axes, try_detect_cone, try_detect_cylinder, try_detect_plane, try_detect_sphere,
    try_detect_torus,
};

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
/// for planarity, cylindricality, sphericality, conicity, and toroidality.
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

    if let Some((ref_normal, centroid, max_dev)) = try_detect_plane(&points, &normals, tolerance) {
        let (u_axis, v_axis) = build_plane_axes(&ref_normal);
        let geometry = SurfaceGeometry::Plane {
            origin: centroid,
            normal: ref_normal,
            u_axis,
            v_axis,
        };
        return Some((AnalyticalKind::Plane, geometry, max_dev));
    }

    if let Some((axis, origin, radius, dev)) = try_detect_cylinder(&points, &normals, tolerance) {
        let geom = SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        };
        return Some((AnalyticalKind::Cylinder, geom, dev));
    }

    if let Some((center, radius, dev)) = try_detect_sphere(&points, &normals, tolerance) {
        let geom = SurfaceGeometry::Sphere {
            center,
            radius,
            axis: crate::geometry::Vec3::Z,
        };
        return Some((AnalyticalKind::Sphere, geom, dev));
    }

    if let Some((apex, axis, half_angle, dev)) =
        try_detect_cone(surface, &points, &normals, tolerance)
    {
        let geom = SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        };
        return Some((AnalyticalKind::Cone, geom, dev));
    }

    if let Some((center, axis, major_radius, minor_radius, dev)) =
        try_detect_torus(surface, &points, &normals, tolerance)
    {
        let geom = SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        };
        return Some((AnalyticalKind::Torus, geom, dev));
    }

    None
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

    let is_nurbs = matches!(sof.geometry, SurfaceGeometry::Nurbs(_));
    if !is_nurbs {
        return None;
    }

    let param_range = sof.geometry.parameter_range();
    let orientation = sof.orientation;

    let (kind, new_geom, max_dev) = detect_analytical_form(&sof.geometry, param_range, tolerance)?;

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
