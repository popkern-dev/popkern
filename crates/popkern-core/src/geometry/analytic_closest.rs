//! Exact closest points on analytic surfaces (PK-737).
//!
//! `SurfaceGeometry::closest_point` routed every analytic surface
//! through the generic bounded Gauss-Newton multi-start. For a point one
//! unit from a cylinder's axis (radius 8) that returned a foot 0.18 off
//! the true radial foot (a rounded-rim fillet junction at r = 7, found
//! by PK-733). Each analytic surface has a closed-form inverse of its
//! own parameterisation:
//! - plane: projection onto the (u, v) axes;
//! - cylinder: the axial coordinate and the radial angle;
//! - sphere: longitude and colatitude in the world frame the
//!   parameterisation uses;
//! - cone: the foot on the generator in the axial half-plane (the apex
//!   when the point lies behind it);
//! - torus: the radial angle, then the tube angle about the tube
//!   centre.
//!
//! **Planes** solve the 2 × 2 Gram system. That is the Gauss-Newton
//! foot to roundoff, in closed form instead of a 29-start search.
//! Enabling it first moved setback-rim end layers by up to 41%, an
//! end-layer ill-conditioning (PK-776). The rim gate splits interior and
//! end layer with the perturbation band measured (PK-773).
//!
//! **Branches.** Periodic parameters are returned in the principal
//! branch `[0, 2π)` of `parameter_range`. **Degenerate points**, where
//! every value of a periodic parameter is equally close (on the axis, at
//! a pole, at a tube centre), take that parameter as 0: a defined choice
//! rather than a solver's accident.

use std::f64::consts::TAU;

use super::analytic_surfaces::orthonormal_basis;
use super::surface::SurfaceGeometry;
use super::{Point3, Vec3};
use popkern_math::PortableMath;

/// Below this (relative to the surface's own scale) a radial direction
/// is degenerate and its angle is taken as 0.
const DEGENERATE: f64 = 1e-14;

/// An angle in the principal branch `[0, 2π)`.
fn principal(angle: f64) -> f64 {
    let a = angle.rem_euclid(TAU);
    if a >= TAU { 0.0 } else { a }
}

/// The angle of `radial` in the frame `orthonormal_basis(axis)`, or 0
/// when `radial` is degenerate against `scale`.
fn angle_about(axis: Vec3, radial: Vec3, scale: f64) -> f64 {
    if radial.length() <= DEGENERATE * scale.max(1.0) {
        return 0.0;
    }
    let (perp1, perp2) = orthonormal_basis(axis);
    principal(radial.dot(&perp2).patan2(radial.dot(&perp1)))
}

/// The exact closest point's `(u, v)` on an analytic surface; `None`
/// for free-form surfaces (NURBS, Gregory, polyhedral spline) and a
/// degenerate frame (a zero axis or plane basis).
pub(crate) fn analytic_closest_point(surface: &SurfaceGeometry, p: &Point3) -> Option<(f64, f64)> {
    match surface {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let aa = axis.dot(axis);
            if aa <= 1e-300 {
                return None;
            }
            let d = *p - *origin;
            let v = d.dot(axis) / aa;
            let radial = d - *axis * v;
            Some((angle_about(*axis, radial, *radius), v))
        }
        SurfaceGeometry::Sphere { center, radius } => {
            let d = *p - *center;
            let len = d.length();
            if len <= DEGENERATE * radius.max(1.0) {
                return Some((0.0, 0.0));
            }
            let v = (d.z / len).clamp(-1.0, 1.0).pacos();
            let rho = (d.x * d.x + d.y * d.y).sqrt();
            let u = if rho <= DEGENERATE * radius.max(1.0) {
                0.0
            } else {
                principal(d.y.patan2(d.x))
            };
            Some((u, v))
        }
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            let a = axis.normalized()?;
            let d = *p - *apex;
            let axial = d.dot(&a);
            let radial = d - a * axial;
            let rho = radial.length();
            // In the (axial, radial) half-plane the generator runs
            // along (cos α, sin α); its foot from (axial, rho), clamped
            // to the apex (the parameter range is v ≥ 0).
            let (s, c) = (half_angle.psin(), half_angle.pcos());
            let along = (axial * c + rho * s).max(0.0);
            let v = along * c;
            Some((angle_about(*axis, radial, rho.max(axial.abs())), v))
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let a = axis.normalized()?;
            let d = *p - *center;
            let axial = d.dot(&a);
            let radial = d - a * axial;
            let rho = radial.length();
            let u = angle_about(*axis, radial, *major_radius);
            // The tube angle about the tube centre at radius `major`.
            let (w_r, w_a) = (rho - major_radius, axial);
            let v = if w_r.phypot(w_a) <= DEGENERATE * minor_radius.max(1.0) {
                0.0
            } else {
                principal(w_a.patan2(w_r))
            };
            Some((u, v))
        }
        SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        } => {
            // The 2 × 2 Gram system (the axes need not be unit or
            // orthogonal).
            let d = *p - *origin;
            let (uu, uv, vv) = (u_axis.dot(u_axis), u_axis.dot(v_axis), v_axis.dot(v_axis));
            let (du, dv) = (d.dot(u_axis), d.dot(v_axis));
            let det = uu * vv - uv * uv;
            if det.abs() <= 1e-300 {
                return None;
            }
            Some(((du * vv - dv * uv) / det, (dv * uu - du * uv) / det))
        }
        SurfaceGeometry::Nurbs(_)
        | SurfaceGeometry::Gregory(_)
        | SurfaceGeometry::PolyhedralSpline(_) => None,
    }
}

#[cfg(test)]
#[path = "analytic_closest_tests.rs"]
mod tests;
