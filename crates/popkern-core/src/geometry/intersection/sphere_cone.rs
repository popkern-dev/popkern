//! Exact sphere × cone intersection for a sphere centred on the cone's
//! axis (PK-1133).
//!
//! With the sphere's centre at signed axial distance `c` from the apex
//! and the cone's points at `apex + a v + v tanα (cos u p1 + sin u p2)`,
//! a point of the cone lies on the sphere when
//! `(v − c)² + v² tan²α = R²`, a quadratic in `v` whose roots are the
//! axial stations of circles of radius `|v| tanα` (either nappe). The
//! pair had no arm and went to the marcher, whose trace of the circle
//! where an ice-cream cone touches its ball drifted 1.9e-6 off at a
//! thousand units from the origin (the tangent direction amplifies a
//! gap error δ to √(2δ/κ)), and with the ball turned at millimetre
//! scale did not finish.
//!
//! The sphere's distance to the cone is `|c| sinα`; within the linear
//! tolerance of `R` they touch along the circle at `v = c cos²α`, and
//! the arm emits nothing, as every arm does for a touch: contact
//! discovery types it beside the curves. A sphere off the axis keeps
//! the marching path.

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::cylinder_cone::uv_hline_of;
use crate::geometry::intersection::plane_sphere::{fit_sphere_uv_curve, nurbs_circle_3d};
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, PortableMath, Vec3};
use crate::tolerance::Tolerance;

/// The sphere × cone answer, or `None` when the pair is not a sphere
/// and a cone or the sphere's centre is off the cone's axis.
pub(crate) fn try_sphere_cone(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    cone_v_window: ((f64, f64), (f64, f64)),
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let (sphere, cone, cone_v, swapped) = match (surface1, surface2) {
        (SurfaceGeometry::Sphere { .. }, SurfaceGeometry::Cone { .. }) => {
            (surface1, surface2, cone_v_window.1, false)
        }
        (SurfaceGeometry::Cone { .. }, SurfaceGeometry::Sphere { .. }) => {
            (surface2, surface1, cone_v_window.0, true)
        }
        _ => return None,
    };
    let SurfaceGeometry::Sphere {
        center,
        radius,
        axis: sphere_axis,
    } = sphere
    else {
        return None;
    };
    let SurfaceGeometry::Cone {
        apex,
        axis,
        half_angle,
    } = cone
    else {
        return None;
    };
    let lin = tolerance.linear();
    let a = axis.normalized()?;
    let to_centre = *center - *apex;
    let c = to_centre.dot(&a);
    if (to_centre - a * c).length() > lin {
        return None;
    }
    let (sin, cos) = (half_angle.psin(), half_angle.pcos());
    if cos.abs() <= tolerance.angular() || sin.abs() <= tolerance.angular() {
        // A cone degenerate to a plane or a line: no circle family.
        return None;
    }
    let tan = sin / cos;
    // A touch along a circle (or none at all): nothing to emit.
    let gap = c.abs() * sin.abs() - radius;
    if gap.abs() <= lin {
        return Some(Vec::new());
    }
    // (1 + tan²) v² − 2 c v + (c² − R²) = 0.
    let k = 1.0 + tan * tan;
    let disc = c * c - k * (c * c - radius * radius);
    if disc < 0.0 {
        return Some(Vec::new());
    }
    let root = disc.sqrt();
    let (v0, v1) = (cone_v.0.min(cone_v.1), cone_v.0.max(cone_v.1));
    let mut out = Vec::new();
    for v in [(c - root) / k, (c + root) / k] {
        // A circle no wider than the tolerance is the apex on the
        // sphere: a point contact, discovery's to type.
        let rho = v.abs() * tan.abs();
        if rho <= lin || v < v0 - lin || v > v1 + lin {
            continue;
        }
        if let Some(curve) = circle_at(
            v,
            rho,
            (*apex, a),
            (*center, *radius, *sphere_axis),
            cone,
            swapped,
            tolerance,
        ) {
            out.push(curve);
        }
    }
    Some(out)
}

/// The circle of radius `rho` at station `v` of the cone, on both
/// charts.
fn circle_at(
    v: f64,
    rho: f64,
    (apex, a): (Point3, Vec3),
    sphere: (Point3, f64, Vec3),
    cone: &SurfaceGeometry,
    swapped: bool,
    tolerance: &Tolerance,
) -> Option<SurfaceSurfaceIntersection> {
    let (p1, p2) = orthonormal_basis(a);
    let curve_3d = nurbs_circle_3d(apex + a * v, p1, p2, rho);
    let on_cone = uv_hline_of(cone, &curve_3d)?;
    let on_sphere = fit_sphere_uv_curve(&curve_3d, sphere, tolerance);
    let (curve_uv1, curve_uv2) = if swapped {
        (on_cone, on_sphere)
    } else {
        (on_sphere, on_cone)
    };
    Some(SurfaceSurfaceIntersection::Curve(CurveOnSurface {
        curve_3d,
        curve_uv1,
        curve_uv2,
        is_closed: true,
    }))
}

#[cfg(test)]
#[path = "sphere_cone_tests.rs"]
mod tests;
