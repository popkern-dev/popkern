//! Analytic cylinder–sphere intersection solver (coaxial case: the
//! sphere's centre lies on the cylinder axis).
//!
//! Three exact regimes, by sphere radius `rs` against cylinder radius
//! `rc`:
//!
//! - `rs < rc`: the sphere sits inside the cylinder — no intersection.
//! - `rs == rc` (within the linear tolerance): the sphere touches the
//!   cylinder along one tangent circle at the centre's axial station.
//!   Under the tangency-silent contract this emits nothing — and that
//!   is precisely the case that motivated this arm: a straight-edge
//!   fillet cylinder meeting its own corner sphere (same rolling-ball
//!   radius, centre on the axis) is tangent along a full circle, so the
//!   marching path seeded ~2100 tangent-classified points and marched
//!   each one to the 20 001-step cap without a crossing certificate
//!   (measured 26 ms per seed, ~55 s per pair call, called three times
//!   per chamfer reconcile pass — the PK-390 "hang").
//! - `rs > rc`: two transversal circles at axial offset
//!   `±sqrt(rs² − rc²)` from the centre. Emitted exactly when the
//!   sphere chart's pole axis is parallel to the cylinder axis (the
//!   circles are then constant-`v` lines in both charts); otherwise
//!   the pair stays on the marching path, which handles transversal
//!   seeds correctly.
//!
//! Non-coaxial pairs return `None` and keep the marching path.

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::curve::CurveEval;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::uv_projection::project_point_to_uv;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// Try to analytically intersect a cylinder and a sphere.
///
/// Returns `Some(intersections)` for a COAXIAL pair (sphere centre on
/// the cylinder axis, in either operand order) whose regime is exact
/// (see the module doc); `None` otherwise.
pub(crate) fn try_cylinder_sphere(
    surface1: &SurfaceGeometry,
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    match (surface1, surface2) {
        (
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            },
            SurfaceGeometry::Sphere {
                center,
                radius: sphere_radius,
                ..
            },
        ) => intersect_impl(
            *origin,
            *axis,
            *radius,
            v1_range,
            *center,
            *sphere_radius,
            surface1,
            surface2,
            tolerance,
            false,
        ),
        (
            SurfaceGeometry::Sphere {
                center,
                radius: sphere_radius,
                ..
            },
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            },
        ) => intersect_impl(
            *origin,
            *axis,
            *radius,
            v2_range,
            *center,
            *sphere_radius,
            surface2,
            surface1,
            tolerance,
            true,
        ),
        _ => None,
    }
}

/// Coaxial-case implementation. `swapped` = surface1 was the sphere.
#[allow(clippy::too_many_arguments)]
fn intersect_impl(
    origin: Point3,
    axis: Vec3,
    rc: f64,
    cyl_v_range: (f64, f64),
    center: Point3,
    rs: f64,
    cyl_surface: &SurfaceGeometry,
    sphere_surface: &SurfaceGeometry,
    tolerance: &Tolerance,
    swapped: bool,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let axis_n = axis.normalized()?;
    let lin = tolerance.linear();
    // Centre on the axis?
    let d = center - origin;
    let h_center = d.dot(&axis_n);
    if (d - axis_n * h_center).length() > lin {
        return None; // off-axis → marching
    }
    let rc = rc.abs();
    let rs = rs.abs();
    if (rs - rc).abs() <= lin {
        // Tangent circle: tangency-silent.
        return Some(Vec::new());
    }
    if rs < rc {
        return Some(Vec::new());
    }
    // Transversal circles need the sphere chart's pole axis along the
    // cylinder axis for an exact constant-v image.
    let pole = match sphere_surface {
        SurfaceGeometry::Sphere { axis, .. } => *axis,
        _ => return None,
    };
    if axis_n.cross(&pole).length() > tolerance.angular() {
        return None;
    }
    let h = (rs * rs - rc * rc).sqrt();
    let (p1, p2) = orthonormal_basis(axis_n);
    let mut curves = Vec::new();
    for station in [h_center + h, h_center - h] {
        if station < cyl_v_range.0 - lin || station > cyl_v_range.1 + lin {
            continue;
        }
        let circle_center = origin + axis_n * station;
        let curve_3d = nurbs_circle(circle_center, p1, p2, rc);
        let uv_cyl = uv_hline_of(cyl_surface, &curve_3d)?;
        let uv_sphere = uv_hline_of(sphere_surface, &curve_3d)?;
        let (uv1, uv2) = if swapped {
            (uv_sphere, uv_cyl)
        } else {
            (uv_cyl, uv_sphere)
        };
        curves.push(SurfaceSurfaceIntersection::Curve(CurveOnSurface {
            curve_3d,
            curve_uv1: uv1,
            curve_uv2: uv2,
            is_closed: true,
        }));
    }
    Some(curves)
}

/// Exact rational-quadratic full circle (parameter 0..1).
fn nurbs_circle(c: Point3, a1: Vec3, a2: Vec3, r: f64) -> NurbsCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            c + a1 * r,
            c + a1 * r + a2 * r,
            c + a2 * r,
            c - a1 * r + a2 * r,
            c - a1 * r,
            c - a1 * r - a2 * r,
            c - a2 * r,
            c + a1 * r - a2 * r,
            c + a1 * r,
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// The UV image of a full circle on a periodic chart whose angular
/// direction follows the circle: a degree-1 constant-v line sweeping
/// one signed period in u, anchored by projecting the 0 and quarter
/// points.
fn uv_hline_of(surface: &SurfaceGeometry, circle: &NurbsCurve) -> Option<NurbsCurve> {
    let p0 = circle.point_at(0.0);
    let pq = circle.point_at(0.25);
    let (u0, v0) = project_point_to_uv(surface, &p0)?;
    let (uq, _) = project_point_to_uv(surface, &pq)?;
    let mut du = uq - u0;
    while du > std::f64::consts::PI {
        du -= std::f64::consts::TAU;
    }
    while du < -std::f64::consts::PI {
        du += std::f64::consts::TAU;
    }
    let sweep = if du >= 0.0 {
        std::f64::consts::TAU
    } else {
        -std::f64::consts::TAU
    };
    Some(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(u0, v0, 0.0), Point3::new(u0 + sweep, v0, 0.0)],
        weights: vec![1.0, 1.0],
    })
}

#[cfg(test)]
#[path = "cylinder_sphere_tests.rs"]
mod tests;
