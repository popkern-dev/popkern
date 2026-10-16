//! Analytic pair dispatch for surface×curve intersection — split
//! from `surface_curve.rs` under the 500-line budget. Line arms hit
//! exact solvers; everything else takes the general NURBS solver
//! (with the dense multi-wrap grid when requested — see
//! `intersect_surface_curve_dense`).

use super::super::{
    line_cone, line_cylinder, line_plane, line_sphere, line_torus, surface_curve_general,
};
use super::SurfaceCurveIntersection;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

/// Route one surface×curve query to its solver.
#[allow(clippy::too_many_arguments)]
pub(super) fn dispatch(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    curve: &CurveGeometry,
    t_range: (f64, f64),
    tolerance: &Tolerance,
    dense_grid: bool,
) -> Vec<SurfaceCurveIntersection> {
    match (surface, curve) {
        (
            SurfaceGeometry::Plane {
                origin,
                normal,
                u_axis,
                v_axis,
            },
            CurveGeometry::Line {
                origin: lo,
                direction,
            },
        ) => line_plane::intersect_line_plane(
            lo, direction, t_range, origin, normal, u_axis, v_axis, u_range, v_range, tolerance,
        ),

        (
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            },
            CurveGeometry::Line {
                origin: lo,
                direction,
            },
        ) => line_cylinder::intersect_line_cylinder(
            lo, direction, t_range, origin, axis, *radius, u_range, v_range, tolerance,
        ),

        (
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            },
            CurveGeometry::Line {
                origin: lo,
                direction,
            },
        ) => line_cone::intersect_line_cone(
            lo,
            direction,
            t_range,
            apex,
            axis,
            *half_angle,
            u_range,
            v_range,
            tolerance,
        ),

        (
            SurfaceGeometry::Sphere {
                center,
                radius,
                axis,
            },
            CurveGeometry::Line {
                origin: lo,
                direction,
            },
        ) => line_sphere::intersect_line_sphere(
            lo, direction, t_range, center, *radius, axis, u_range, v_range, tolerance,
        ),

        (
            SurfaceGeometry::Torus {
                center,
                axis,
                major_radius,
                minor_radius,
            },
            CurveGeometry::Line {
                origin: lo,
                direction,
            },
        ) => line_torus::intersect_line_torus(
            lo,
            direction,
            t_range,
            center,
            axis,
            *major_radius,
            *minor_radius,
            u_range,
            v_range,
            tolerance,
        ),

        (
            SurfaceGeometry::Plane {
                origin,
                normal,
                u_axis,
                v_axis,
            },
            CurveGeometry::Circle { .. } | CurveGeometry::Arc { .. },
        ) => super::super::plane_circle::intersect_plane_circle(
            super::super::plane_circle::PlaneFrame {
                origin,
                normal,
                u_axis,
                v_axis,
            },
            curve,
            t_range,
            (u_range, v_range),
            tolerance,
        )
        .unwrap_or_default(),

        _ => surface_curve_general::intersect_surface_curve_general(
            surface, u_range, v_range, curve, t_range, tolerance, dense_grid,
        ),
    }
}
