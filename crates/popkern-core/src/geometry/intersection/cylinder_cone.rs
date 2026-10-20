//! Analytic cylinder–cone intersection solver (coaxial case).
//!
//! A cylinder and a cone sharing an axis intersect in the exact
//! CIRCLE(s) at the axial station(s) where the cone's radius equals
//! the cylinder's. The marching solver fragmented this everyday
//! configuration (countersunk bore): each seed traced a partial arc
//! that never closed, the per-pair curve cap kept 20 rotated copies,
//! and consolidation left a gapped loop the downstream chart walk
//! could not balance (PK-147). Non-coaxial pairs return `None` and
//! keep the marching path.

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::curve::CurveEval;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::uv_projection::project_point_to_uv;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

/// Try to analytically intersect a cylinder and a cone.
///
/// Returns `Some(intersections)` when the pair is a COAXIAL cylinder
/// and cone (in either order); `None` for other combinations —
/// including parallel-offset and skew axes, which stay on the
/// marching path.
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_cylinder_cone(
    surface1: &SurfaceGeometry,
    _u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    _u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    match (surface1, surface2) {
        (
            SurfaceGeometry::Cylinder {
                origin: co,
                axis: ca,
                radius: cr,
            },
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            },
        ) => intersect_impl(
            *co,
            *ca,
            *cr,
            v1_range,
            *apex,
            *axis,
            *half_angle,
            v2_range,
            surface1,
            surface2,
            tolerance,
            false,
        ),
        (
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            },
            SurfaceGeometry::Cylinder {
                origin: co,
                axis: ca,
                radius: cr,
            },
        ) => intersect_impl(
            *co,
            *ca,
            *cr,
            v2_range,
            *apex,
            *axis,
            *half_angle,
            v1_range,
            surface2,
            surface1,
            tolerance,
            true,
        ),
        _ => None,
    }
}

/// Whether a cylinder of radius `cr` and a cone whose axis is parallel
/// to it, `offset` away, cannot meet over the cone's axial window: the
/// cone's circles there (radii |v| tan h) all lie strictly inside the
/// cylinder, all strictly outside it, or all strictly around it, by 10
/// linear tolerances. A cylinder window does not bound the test, so
/// it is conservative.
fn parallel_offset_disjoint(
    cr: f64,
    offset: f64,
    half_angle: f64,
    cone_v_range: (f64, f64),
    tolerance: &Tolerance,
) -> bool {
    let (v0, v1) = cone_v_range;
    if !(v0.is_finite() && v1.is_finite() && v0 <= v1) {
        return false;
    }
    let tan_h = half_angle.ptan().abs();
    let r_max = v0.abs().max(v1.abs()) * tan_h;
    // The window crossing the apex holds the zero-radius circle.
    let r_min = if v0 <= 0.0 && v1 >= 0.0 {
        0.0
    } else {
        v0.abs().min(v1.abs()) * tan_h
    };
    let margin = 10.0 * tolerance.linear();
    let inside = offset + r_max < cr - margin;
    let outside = offset - r_max > cr + margin;
    let around = offset + cr < r_min - margin;
    inside || outside || around
}

/// Coaxial-case implementation. `swapped` = surface1 was the cone.
#[allow(clippy::too_many_arguments)]
fn intersect_impl(
    co: Point3,
    ca: Vec3,
    cr: f64,
    cyl_v_range: (f64, f64),
    apex: Point3,
    axis: Vec3,
    half_angle: f64,
    cone_v_range: (f64, f64),
    cyl_surface: &SurfaceGeometry,
    cone_surface: &SurfaceGeometry,
    tolerance: &Tolerance,
    swapped: bool,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let ca_n = ca.normalized()?;
    let ax_n = axis.normalized()?;
    // Parallel axes?
    if ca_n.cross(&ax_n).length() > tolerance.angular() {
        return None; // skew → marching
    }
    // Apex on the cylinder axis?
    let d = apex - co;
    let offset = (d - ca_n * d.dot(&ca_n)).length();
    if offset > tolerance.linear() {
        // Parallel-offset axes: exactly empty when the cone's window
        // (a frustum, radius |v| tan h over the window's v) lies wholly
        // inside the solid cylinder, wholly outside it, or wholly
        // around it; otherwise marching. Every tapped hole's entry
        // chamfer cone beside a flange rim paid 80-150 ms of marching
        // to find nothing (PK-876).
        return parallel_offset_disjoint(cr, offset, half_angle, cone_v_range, tolerance)
            .then(Vec::new);
    }
    let tan_h = half_angle.ptan();
    if tan_h.abs() < tolerance.angular() {
        // Degenerate cone (≈ cylinder): no transversal circle; the
        // coincident detector owns the tangent/coincident cases.
        return None;
    }
    // Stations where the cone's radius reaches the cylinder's, one
    // per nappe (cone v is the SIGNED axial distance from the apex).
    let v_mag = cr / tan_h.abs();
    let mut curves = Vec::new();
    for v_cone in [v_mag, -v_mag] {
        if v_cone < cone_v_range.0 || v_cone > cone_v_range.1 {
            continue;
        }
        let center = apex + ax_n * v_cone;
        let v_cyl = (center - co).dot(&ca_n);
        if v_cyl < cyl_v_range.0 || v_cyl > cyl_v_range.1 {
            continue;
        }
        let (p1, p2) = orthonormal_basis(ca_n);
        let curve_3d = nurbs_circle(center, p1, p2, cr);
        // UV images: horizontal lines in each chart, anchored by
        // PROJECTION so the phase and direction are exact for any
        // axis-orientation combination (the two surfaces' angular
        // bases need not agree).
        let uv_cyl = uv_hline_of(cyl_surface, &curve_3d)?;
        let uv_cone = uv_hline_of(cone_surface, &curve_3d)?;
        let (uv1, uv2) = if swapped {
            (uv_cone.clone(), uv_cyl.clone())
        } else {
            (uv_cyl.clone(), uv_cone.clone())
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

/// The UV image of a full circle on a coaxial periodic chart: a
/// degree-1 line at constant v sweeping one signed period in u,
/// anchored and oriented by projecting the curve's 0 and quarter
/// points into the chart.
fn uv_hline_of(surface: &SurfaceGeometry, circle: &NurbsCurve) -> Option<NurbsCurve> {
    let p0 = circle.point_at(0.0);
    let pq = circle.point_at(0.25);
    let (u0, v0) = project_point_to_uv(surface, &p0)?;
    let (uq, _) = project_point_to_uv(surface, &pq)?;
    // Quarter-turn unwrap: the angular step is ±π/2.
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
#[path = "cylinder_cone_offset_tests.rs"]
mod offset_tests;
