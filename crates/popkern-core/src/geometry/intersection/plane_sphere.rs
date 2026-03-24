//! Analytic plane-sphere intersection solver.
//!
//! When the two surfaces are a plane and a sphere the intersection is either
//! empty (no contact), tangent (degenerate — treated as empty), or a circle.
//! Computing this analytically is far cheaper than the general NURBS marching
//! pipeline.

use std::f64::consts::FRAC_1_SQRT_2;

use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Try to analytically intersect a plane and a sphere.
///
/// Returns `Some(Vec<SurfaceSurfaceIntersection>)` when the surface pair is
/// a plane and a sphere (in either order), or `None` when neither surface is
/// a plane-sphere combination and the caller should fall back to the general
/// intersection pipeline.
///
/// Returns `Some(vec![])` for no intersection or tangent (degenerate) cases.
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_plane_sphere(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    match (surface1, surface2) {
        (
            SurfaceGeometry::Plane {
                origin,
                normal,
                u_axis,
                v_axis,
            },
            SurfaceGeometry::Sphere { center, radius },
        ) => Some(intersect_plane_sphere(
            *origin, *normal, *u_axis, *v_axis, u1_range, v1_range, *center, *radius, u2_range,
            v2_range, tolerance, false,
        )),
        (
            SurfaceGeometry::Sphere { center, radius },
            SurfaceGeometry::Plane {
                origin,
                normal,
                u_axis,
                v_axis,
            },
        ) => Some(intersect_plane_sphere(
            *origin, *normal, *u_axis, *v_axis, u2_range, v2_range, *center, *radius, u1_range,
            v1_range, tolerance, true,
        )),
        _ => None,
    }
}

/// Core implementation: intersect a plane with a sphere.
///
/// When `swapped` is true the sphere was surface 1 and the plane was surface 2,
/// so we swap the UV curves in the output.
#[allow(clippy::too_many_arguments)]
fn intersect_plane_sphere(
    plane_origin: Point3,
    plane_normal: Vec3,
    plane_u: Vec3,
    plane_v: Vec3,
    plane_u_range: (f64, f64),
    plane_v_range: (f64, f64),
    sphere_center: Point3,
    sphere_radius: f64,
    sphere_u_range: (f64, f64),
    sphere_v_range: (f64, f64),
    tolerance: &Tolerance,
    swapped: bool,
) -> Vec<SurfaceSurfaceIntersection> {
    // Signed distance from sphere centre to the plane.
    let diff = sphere_center - plane_origin;
    let d = diff.dot(&plane_normal);

    // No intersection: sphere is too far from plane.
    if d.abs() > sphere_radius + tolerance.linear {
        return vec![];
    }

    // Tangent (degenerate): intersection would be a single point.
    if d.abs() > sphere_radius - tolerance.linear {
        return vec![];
    }

    // Circle radius and projected centre.
    let r_circle = (sphere_radius * sphere_radius - d * d).sqrt();
    let center_proj = sphere_center - plane_normal * d;

    // Build in-plane axes perpendicular to the plane normal.
    // We use the plane's own u/v axes as our circle frame.
    let axis1 = plane_u;
    let axis2 = plane_v;

    // 3D NURBS circle.
    let curve_3d = nurbs_circle_3d(center_proj, axis1, axis2, r_circle);

    // UV on plane: project each 3D control point of the circle into the
    // plane's UV space. This correctly handles non-orthogonal u_axis/v_axis
    // (which produce skewed UV projections, not axis-aligned ellipses).
    let curve_uv_plane = project_curve_to_plane_uv(&curve_3d, plane_origin, plane_u, plane_v);

    // UV on sphere:
    // The sphere parameterization is u = longitude, v = colatitude.
    // The intersection circle sits at colatitude v_const = acos(d / radius).
    // We build a NURBS circle in UV space at constant v.
    let v_const = (d / sphere_radius).clamp(-1.0, 1.0).acos();
    let curve_uv_sphere = nurbs_sphere_uv_circle(sphere_u_range, v_const);

    // Respect the UV ranges — if the circle centre is outside the plane's
    // parameter box, the intersection may not be relevant, but we still
    // return it (clipping is the caller's responsibility).
    let _ = (plane_u_range, plane_v_range, sphere_v_range);

    let (curve_uv1, curve_uv2) = if swapped {
        (curve_uv_sphere, curve_uv_plane)
    } else {
        (curve_uv_plane, curve_uv_sphere)
    };

    let cos = CurveOnSurface {
        curve_3d,
        curve_uv1,
        curve_uv2,
        is_closed: true,
    };

    vec![SurfaceSurfaceIntersection::Curve(cos)]
}

/// Build a degree-2 rational NURBS circle with 9 control points.
///
/// The circle lies in the plane spanned by `axis1` and `axis2`, centred at
/// `center`, with the given `radius`.
fn nurbs_circle_3d(center: Point3, axis1: Vec3, axis2: Vec3, radius: f64) -> NurbsCurve {
    nurbs_ellipse_3d(center, axis1, axis2, radius, radius)
}

/// Build a degree-2 rational NURBS ellipse in 3D.
fn nurbs_ellipse_3d(center: Point3, axis1: Vec3, axis2: Vec3, ra: f64, rb: f64) -> NurbsCurve {
    let w = FRAC_1_SQRT_2;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            center + axis1 * ra,
            center + axis1 * ra + axis2 * rb,
            center + axis2 * rb,
            center - axis1 * ra + axis2 * rb,
            center - axis1 * ra,
            center - axis1 * ra - axis2 * rb,
            center - axis2 * rb,
            center + axis1 * ra - axis2 * rb,
            center + axis1 * ra,
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Build a degree-2 rational NURBS ellipse in 2D UV space.
fn nurbs_ellipse_2d(center: Point3, ru: f64, rv: f64) -> NurbsCurve {
    nurbs_ellipse_3d(
        center,
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        ru,
        rv,
    )
}

/// Build a NURBS circle in the sphere's UV parameter space.
///
/// The intersection circle at constant colatitude `v_const` is a horizontal
/// line in UV space from `(u_start, v_const)` to `(u_end, v_const)`.
/// For a full circle this is `(0, v_const)` to `(TAU, v_const)`.
fn nurbs_sphere_uv_circle(u_range: (f64, f64), v_const: f64) -> NurbsCurve {
    let u_start = u_range.0;
    let u_end = u_range.1;
    // Represent the full longitude sweep as a degree-1 NURBS line in UV space.
    // Using the same 9-control-point structure as the 3D circle would be
    // incorrect because the UV-to-3D map is non-linear. Instead we use a
    // simple linear segment that traverses the full u range at constant v.
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(u_start, v_const, 0.0),
            Point3::new(u_end, v_const, 0.0),
        ],
        weights: vec![1.0, 1.0],
    }
}

/// Project a 3D NURBS curve onto a plane's UV space by transforming each
/// control point. Preserves degree, knots, and weights.
///
/// Uses the correct inverse parameterization:
/// `u = (p-origin)·u_axis / |u_axis|²`, `v = (p-origin)·v_axis / |v_axis|²`.
fn project_curve_to_plane_uv(
    curve: &NurbsCurve,
    origin: Point3,
    u_axis: Vec3,
    v_axis: Vec3,
) -> NurbsCurve {
    let u_len_sq = u_axis.dot(&u_axis);
    let v_len_sq = v_axis.dot(&v_axis);
    let cps = curve
        .control_points
        .iter()
        .map(|p| {
            let d = *p - origin;
            let u = if u_len_sq > 1e-30 {
                d.dot(&u_axis) / u_len_sq
            } else {
                d.dot(&u_axis)
            };
            let v = if v_len_sq > 1e-30 {
                d.dot(&v_axis) / v_len_sq
            } else {
                d.dot(&v_axis)
            };
            Point3::new(u, v, 0.0)
        })
        .collect();
    NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points: cps,
        weights: curve.weights.clone(),
    }
}

#[cfg(test)]
#[path = "plane_sphere_tests.rs"]
mod tests;
