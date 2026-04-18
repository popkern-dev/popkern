//! §5.3 item #6 — surface-aware bridge-curve construction for the
//! blend-boundary face-edge segments.
//!
//! Previously the hex/N-gon boundary's face-edge bridges were
//! straight-line NURBS on planar faces. For corners where one or
//! more adjacent faces are curved (cylinder, sphere, cone, torus,
//! NURBS, Gregory), a straight 3D line leaves the surface and the
//! three-quad Gregory dispatch degrades to the polar-NURBS fallback.
//!
//! This module emits a cubic-Bezier bridge that **lies on the
//! supporting surface** by interpolating in the surface's UV
//! parameter space. The resulting NURBS curve is C² on analytic
//! supports (Plane, Cylinder, Sphere, Cone, Torus) and an excellent
//! approximation on NURBS / Gregory supports (parametric linear in
//! UV, which is exact only for Plane but close enough on smooth
//! surfaces near a small blend region).
//!
//! # Periodic unwrap
//!
//! Cylinders, spheres, cones, and tori have periodic `u`-parameter
//! spaces (wrap-around angles). Naive linear interpolation in UV
//! can pick the *long* way around — we unwrap `u_b` modulo the
//! period so the interpolation takes the shorter arc.
//!
//! # Degenerate fallback
//!
//! When projection fails or the surface is otherwise ill-defined
//! at `pa` / `pb`, we fall back to a straight-line cubic Bezier
//! with the 3D midpoint control points. This matches the pre-item
//! #6 behaviour — still a usable bridge, just not surface-hugging.

use std::f64::consts::TAU;

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3};

use super::boundary_hex::cubic_nurbs_from_control_points;

/// Short human-readable name for a surface variant — used for UI /
/// log display and as the `surface_kind` field on the
/// [`BlendWarning::SurfaceBridgeFallback`](crate::error::BlendWarning::SurfaceBridgeFallback)
/// variant.
pub(crate) fn surface_kind_name(surface: &SurfaceGeometry) -> &'static str {
    match surface {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cylinder",
        SurfaceGeometry::Sphere { .. } => "Sphere",
        SurfaceGeometry::Cone { .. } => "Cone",
        SurfaceGeometry::Torus { .. } => "Torus",
        SurfaceGeometry::Nurbs(_) => "NURBS",
        SurfaceGeometry::Gregory(_) => "Gregory",
    }
}

/// Distance from a 3D point to a surface, defined as
/// `|p − surface.point_at(u, v)|` where `(u, v)` is the UV
/// projection of `p` onto the surface. Returns `f64::INFINITY` if
/// projection fails.
///
/// Used to filter which candidate face a blend-boundary bridge
/// should live on when multiple faces are incident to the blend
/// region.
pub(crate) fn point_on_surface_distance(p: &Point3, surface: &SurfaceGeometry) -> f64 {
    let Some((u, v)) = project_point_to_uv(surface, p) else {
        return f64::INFINITY;
    };
    let q = surface.point_at(u, v);
    p.distance_to(&q)
}

/// Build a cubic-Bezier bridge from `pa` to `pb` that lies on the
/// given surface. Returns `None` if `pa` or `pb` projects too far
/// from the surface (`> tol_linear`).
///
/// # Algorithm
///
/// 1. Project `pa` and `pb` to UV via
///    [`project_point_to_uv`].
/// 2. Unwrap `u` (and for a torus, `v`) to the shortest periodic
///    branch between the two endpoints.
/// 3. Sample the surface at `t ∈ {0, 1/3, 2/3, 1}` in UV space.
/// 4. Solve for the cubic Bezier's interior control points so the
///    curve passes through all four samples (same 4-point
///    interpolation used in `fit_cubic_bezier_through`).
///
/// The resulting curve has:
/// - Endpoints at exactly `pa_projected` and `pb_projected` (which
///   equal `pa`, `pb` within `tol_linear`).
/// - Two intermediate control points positioned so that the
///   Bezier samples at `t = 1/3` and `t = 2/3` match the surface's
///   UV-linear midpoints.
pub(crate) fn build_surface_bridge(
    pa: &Point3,
    pb: &Point3,
    surface: &SurfaceGeometry,
    tol_linear: f64,
) -> Option<NurbsCurve> {
    let (u_a, v_a) = project_point_to_uv(surface, pa)?;
    let (u_b, v_b) = project_point_to_uv(surface, pb)?;
    // Verify the projections land close to the input points —
    // rules out off-face endpoints.
    if pa.distance_to(&surface.point_at(u_a, v_a)) > tol_linear
        || pb.distance_to(&surface.point_at(u_b, v_b)) > tol_linear
    {
        return None;
    }

    let (u_a, u_b, v_a, v_b) = unwrap_periodic_uv(surface, u_a, v_a, u_b, v_b);

    let lerp = |a: f64, b: f64, t: f64| (1.0 - t) * a + t * b;
    let s0 = surface.point_at(u_a, v_a);
    let s1 = surface.point_at(lerp(u_a, u_b, 1.0 / 3.0), lerp(v_a, v_b, 1.0 / 3.0));
    let s2 = surface.point_at(lerp(u_a, u_b, 2.0 / 3.0), lerp(v_a, v_b, 2.0 / 3.0));
    let s3 = surface.point_at(u_b, v_b);

    let cps = cubic_bezier_through_4_samples(s0, s1, s2, s3);
    cubic_nurbs_from_control_points(cps).ok()
}

/// Solve for a cubic Bezier `P_0..P_3` so the curve passes through
/// the four sample points `S_0..S_3` at parameters `t = 0, 1/3,
/// 2/3, 1`. `P_0 = S_0`, `P_3 = S_3`; the interior control points
/// come from the 2×2 linear system (same formula as
/// `fit_cubic_bezier_through` — determinant 108).
fn cubic_bezier_through_4_samples(s0: Point3, s1: Point3, s2: Point3, s3: Point3) -> [Point3; 4] {
    let axis = |q1c: f64, q2c: f64, p0c: f64, p3c: f64| -> (f64, f64) {
        let a = 27.0 * q1c - 8.0 * p0c - p3c;
        let b = 27.0 * q2c - p0c - 8.0 * p3c;
        ((2.0 * a - b) / 18.0, (2.0 * b - a) / 18.0)
    };
    let (p1x, p2x) = axis(s1.x, s2.x, s0.x, s3.x);
    let (p1y, p2y) = axis(s1.y, s2.y, s0.y, s3.y);
    let (p1z, p2z) = axis(s1.z, s2.z, s0.z, s3.z);
    [
        s0,
        Point3::new(p1x, p1y, p1z),
        Point3::new(p2x, p2y, p2z),
        s3,
    ]
}

/// Unwrap UV parameters on periodic surfaces to minimise the
/// parameter-space distance between `(u_a, v_a)` and `(u_b, v_b)`.
///
/// - `Cylinder`, `Cone`, `Sphere`: `u` period = `2π`; `v` is
///   non-periodic.
/// - `Torus`: both `u` and `v` are periodic with period `2π`.
/// - `Plane`, `NURBS`, `Gregory`: no periodicity.
fn unwrap_periodic_uv(
    surface: &SurfaceGeometry,
    u_a: f64,
    v_a: f64,
    u_b: f64,
    v_b: f64,
) -> (f64, f64, f64, f64) {
    let (u_periodic, v_periodic) = surface.is_periodic();
    let u_b = if u_periodic {
        unwrap_to_nearest(u_a, u_b, TAU)
    } else {
        u_b
    };
    let v_b = if v_periodic {
        unwrap_to_nearest(v_a, v_b, TAU)
    } else {
        v_b
    };
    (u_a, u_b, v_a, v_b)
}

/// Shift `target` by multiples of `period` so that it is as close
/// as possible to `ref_val` on the real line. Idempotent for
/// non-periodic inputs (caller should pass only periodic coords).
fn unwrap_to_nearest(ref_val: f64, target: f64, period: f64) -> f64 {
    let diff = target - ref_val;
    let shifts = (diff / period).round();
    target - shifts * period
}

/// Straight-line cubic-Bezier fallback through 2 explicit endpoints
/// (used when surface projection fails and the caller wants a
/// best-effort bridge that at least interpolates `pa` / `pb`).
pub(super) fn build_straight_line_bridge(pa: &Point3, pb: &Point3) -> Option<NurbsCurve> {
    let lerp = |t: f64| {
        Point3::new(
            (1.0 - t) * pa.x + t * pb.x,
            (1.0 - t) * pa.y + t * pb.y,
            (1.0 - t) * pa.z + t * pb.z,
        )
    };
    let cps = [*pa, lerp(1.0 / 3.0), lerp(2.0 / 3.0), *pb];
    cubic_nurbs_from_control_points(cps).ok()
}

#[cfg(test)]
#[path = "surface_bridge_tests.rs"]
mod tests;
