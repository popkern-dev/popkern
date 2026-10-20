//! Closed sweep walls in periodic form (PK-892,
//! `docs/architecture/closed-nurbs-rings.md`).
//!
//! A one-span closed profile (a wire circle) sweeps to a wall that
//! closes on itself in u. Written clamped, the wall is a chart holding
//! its seam twice; written in periodic form it is a ring like a
//! cylinder's lateral face. Only a seam the profile crosses smoothly
//! is converted: a profile that closes at a corner keeps the clamped
//! chart, whose seam is then a real crease.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsSurface, Point3};
use popkern_core::tolerance::Tolerance;

/// Sine bound for the profile's two end tangents to count as one line
/// across the seam: the rows a producer writes for a smooth closed
/// profile agree to roundoff.
const SEAM_TANGENT_SINE: f64 = 1e-9;

/// The wall geometry of `surface`: in periodic form in u when it
/// closes there within the linear tolerance and every column crosses
/// the seam tangent continuously; otherwise `surface` unchanged.
pub(crate) fn ring_wall(surface: &NurbsSurface, tolerance: Tolerance) -> SurfaceGeometry {
    if smooth_u_seam(surface)
        && let Some(ring) = surface.periodic_u_from_closed_clamped(tolerance.linear())
    {
        return SurfaceGeometry::Nurbs(ring);
    }
    SurfaceGeometry::Nurbs(surface.clone())
}

/// Whether each u-column leaves its first row along the line it
/// arrives at its last row on (the end control legs are parallel and
/// point the same way).
fn smooth_u_seam(surface: &NurbsSurface) -> bool {
    let (n, m) = (surface.num_u, surface.num_v);
    if n < 3 || surface.control_points.len() != n * m {
        return false;
    }
    let at = |i: usize, j: usize| -> Point3 { surface.control_points[i * m + j] };
    (0..m).all(|j| {
        let leave = at(1, j) - at(0, j);
        let arrive = at(n - 1, j) - at(n - 2, j);
        match (leave.normalized(), arrive.normalized()) {
            (Some(a), Some(b)) => a.dot(&b) > 0.0 && a.cross(&b).length() <= SEAM_TANGENT_SINE,
            _ => false,
        }
    })
}
