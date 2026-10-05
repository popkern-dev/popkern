//! Analytic plane–torus intersection solver (PK-279).
//!
//! Results by configuration:
//! - **Axis-perpendicular plane** — up to two exact circles (outer and
//!   inner tube crossings); tangency at the tube extremes is silent.
//! - **Axis-containing plane** — two exact tube circles at ±R along
//!   the in-plane direction perpendicular to the axis.
//! - **Villarceau bitangent plane** (through the centre, inclined so
//!   `|cos| = √(R²−r²)/R` to the axis) — the two exact Villarceau
//!   circles of radius R, verified against the torus before emission;
//!   a failed verification falls through to the general arm.
//! - **General position** — the quartic spiric section, solved by the
//!   sampled-analytic branch solver in [`super::plane_torus_spiric`]:
//!   closed-form tube-angle solutions per axis angle, complete branch
//!   enumeration from the trig-quartic feasibility function, loops
//!   fitted through the same machinery as the marching path.
//!
//! The pair is claimed in all orientations; `None` is returned only
//! for non-plane/torus pairs or when a numeric guard (non-finite
//! input, chart failure, fit failure) declines — in which case the
//! marching path remains the documented fallback boundary.

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::circle_of_revolution::{
    circle_curve_on, station_in_window, torus_tube_station,
};
use crate::geometry::intersection::plane_torus_spiric::try_spiric_sections;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;
use std::f64::consts::TAU;

/// Try to analytically intersect a plane and a torus.
///
/// Returns `Some(intersections)` when the surface pair is a plane and
/// a torus (in either order) and the solver produced a verdict —
/// including the empty and tangency-silent outcomes. Returns `None`
/// for other surface combinations, or when a numeric guard declines
/// and the marching fallback should run.
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_plane_torus(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    match (surface1, surface2) {
        (SurfaceGeometry::Plane { .. }, SurfaceGeometry::Torus { .. }) => intersect_impl(
            surface1,
            (u1_range, v1_range),
            surface2,
            (u2_range, v2_range),
            tolerance,
            false,
        ),
        (SurfaceGeometry::Torus { .. }, SurfaceGeometry::Plane { .. }) => intersect_impl(
            surface2,
            (u2_range, v2_range),
            surface1,
            (u1_range, v1_range),
            tolerance,
            true,
        ),
        _ => None,
    }
}

/// `swapped` = the torus was operand 1.
fn intersect_impl(
    plane: &SurfaceGeometry,
    plane_windows: ((f64, f64), (f64, f64)),
    torus: &SurfaceGeometry,
    torus_windows: ((f64, f64), (f64, f64)),
    tolerance: &Tolerance,
    swapped: bool,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let torus_v_range = torus_windows.1;
    let SurfaceGeometry::Plane {
        origin: po, normal, ..
    } = plane
    else {
        return None;
    };
    let SurfaceGeometry::Torus {
        center,
        axis,
        major_radius,
        minor_radius,
    } = torus
    else {
        return None;
    };
    let n_hat = normal.normalized()?;
    let axis_n = axis.normalized()?;
    let (big_r, r) = (*major_radius, *minor_radius);
    if !big_r.is_finite() || !r.is_finite() || r <= 0.0 {
        return None;
    }
    let lin = tolerance.linear();
    let ang = tolerance.angular();
    let cos_ax = n_hat.dot(&axis_n);
    // Signed distance of the torus centre from the plane.
    let dc = (*center - *po).dot(&n_hat);

    let (s1, w1, s2, w2) = if swapped {
        (torus, torus_windows, plane, plane_windows)
    } else {
        (plane, plane_windows, torus, torus_windows)
    };

    // Snap gates are GEOMETRIC, not angular: an angular gate of 1e-4
    // admits tilts whose snapped circles sit R·tilt ≈ 5e-3 off the
    // true intersection (measured on a 2e-3-tilted plane over the
    // R=2 fixture) — far outside linear tolerance — while the spiric
    // arm handles those planes exactly. A configuration snaps only
    // when the snap's worst-case geometric error stays within the
    // linear tolerance.
    let sin_ax = (1.0 - cos_ax * cos_ax).max(0.0).sqrt();

    // ── Axis-perpendicular plane: circles at fixed tube stations ──
    if sin_ax * (big_r.abs() + r) <= lin {
        // Axial offset of the plane from the torus centre.
        let h = -dc / cos_ax;
        if (h.abs() - r).abs() <= lin {
            return Some(Vec::new()); // tangent at tube extreme: silent
        }
        if h.abs() > r {
            return Some(Vec::new());
        }
        let s = (r * r - h * h).sqrt();
        let (a1, a2) = orthonormal_basis(axis_n);
        let mut curves = Vec::new();
        for rho in [big_r + s, big_r - s] {
            if rho.abs() <= lin {
                continue;
            }
            let station = torus_tube_station(rho.abs(), h, big_r);
            if !station_in_window(station, torus_v_range, Some(TAU)) {
                continue;
            }
            let (f1, f2) = if rho >= 0.0 { (a1, a2) } else { (-a1, -a2) };
            let curve = circle_curve_on(
                s1,
                w1,
                s2,
                w2,
                *center + axis_n * h,
                f1,
                f2,
                rho.abs(),
                tolerance,
            )?;
            curves.push(SurfaceSurfaceIntersection::Curve(curve));
        }
        return Some(curves);
    }

    // ── Axis-containing plane: two exact tube circles ──
    if cos_ax.abs() * r <= lin && dc.abs() <= lin {
        let m_hat = axis_n.cross(&n_hat).normalized()?;
        let mut curves = Vec::new();
        for sign in [1.0f64, -1.0] {
            let curve = circle_curve_on(
                s1,
                w1,
                s2,
                w2,
                *center + m_hat * (sign * big_r),
                m_hat * sign,
                axis_n,
                r,
                tolerance,
            )?;
            curves.push(SurfaceSurfaceIntersection::Curve(curve));
        }
        return Some(curves);
    }

    // ── Villarceau bitangent plane: two exact circles of radius R ──
    if dc.abs() <= lin && big_r > r + lin {
        let villarceau_cos = ((big_r * big_r - r * r).sqrt() / big_r).clamp(0.0, 1.0);
        if (cos_ax.abs() - villarceau_cos).abs() <= ang
            && let Some(curves) =
                villarceau_circles(s1, w1, s2, w2, *center, axis_n, n_hat, big_r, r, tolerance)
        {
            return Some(curves);
        }
        // Verification failed (near-miss inside the angular gate):
        // fall through to the exact general arm.
    }

    // ── General position: sampled-analytic spiric sections ──
    try_spiric_sections(plane, torus, torus_windows, swapped, tolerance)
}

/// Construct and VERIFY the two Villarceau circles. The detection gate
/// upstream is angular (loose); acceptance here is the honest bar —
/// every sample must lie on the torus within the linear tolerance.
/// Returns `None` when verification fails so the caller can use the
/// general arm instead.
#[allow(clippy::too_many_arguments)]
fn villarceau_circles(
    s1: &SurfaceGeometry,
    w1: ((f64, f64), (f64, f64)),
    s2: &SurfaceGeometry,
    w2: ((f64, f64), (f64, f64)),
    center: Point3,
    axis_n: Vec3,
    n_hat: Vec3,
    big_r: f64,
    r: f64,
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    // Tilt-line direction (intersection of the cutting plane with the
    // equatorial plane) — the circle centres sit at ±r along it.
    let m_hat = n_hat.cross(&axis_n).normalized()?;
    // In-plane companion, perpendicular to the tilt line.
    let w_hat = n_hat.cross(&m_hat).normalized()?;
    let mut curves = Vec::new();
    for sign in [1.0f64, -1.0] {
        let c = center + m_hat * (sign * r);
        // Verify the candidate lies on the torus before accepting.
        let samples = 16;
        for i in 0..samples {
            let theta = TAU * (i as f64) / (samples as f64);
            let p = c + m_hat * (big_r * theta.pcos()) + w_hat * (big_r * theta.psin());
            let d = p - center;
            let h = d.dot(&axis_n);
            let rho = (d - axis_n * h).length();
            let residual = (((rho - big_r) * (rho - big_r) + h * h).sqrt() - r).abs();
            if residual > tolerance.linear() {
                return None;
            }
        }
        let curve = circle_curve_on(s1, w1, s2, w2, c, m_hat, w_hat, big_r, tolerance)?;
        // The pcurve fit's mid-knot verification bounds chart error;
        // additionally require the 3D curve to close on itself as a
        // sanity check on the frame.
        let p0 = curve.curve_3d.point_at(0.0);
        let p1 = curve.curve_3d.point_at(1.0);
        if p0.distance_to(&p1) > tolerance.linear() {
            return None;
        }
        curves.push(SurfaceSurfaceIntersection::Curve(curve));
    }
    Some(curves)
}

#[cfg(test)]
#[path = "plane_torus_tests.rs"]
mod tests;
