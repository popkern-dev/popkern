//! Tangent-match end cap (PK-536): the Class-A step beyond the
//! rolling-ball sphere.
//!
//! A [`CapKind::Rounded`](crate::types::CapKind::Rounded) sphere is G1
//! to the strip along its end section and to both supports at the
//! contacts. However, it jumps in curvature across the section: the
//! strip's rulings are straight, while the sphere's meridians curve at
//! 1/r.
//!
//! The tangent-match cap is a surface of revolution about the end
//! tangent line through the ball centre, swept over the section arc's
//! angular span. Its meridian profile is a quartic Bézier in (axial,
//! radial) units of r:
//!
//! ```text
//! (0, 1) (1/2, 1) (1, 1) (1, 1/2) (1, 0)
//! ```
//!
//! - The first three control points are collinear on the section
//!   circle's radius. The meridian therefore leaves the section along
//!   the ruling with zero curvature, and the cap is G2 to a cylinder
//!   strip across the section. It is G1 to any rolling-ball strip,
//!   since the section lies in the ball's normal plane.
//! - At the contacts, the cap's tangent plane is the strip's, which is
//!   the support's.
//! - The last three control points are collinear across the axis, so
//!   the pole is flat and smooth.
//! - The pole sits r beyond the station, where the Rounded pole does.
//!
//! The surface is exact NURBS: rational quadratic in the angle u (one
//! arc segment up to 90°, two beyond), and polynomial quartic in the
//! profile v, with v = 0 at the section and v = 1 at the pole. Its
//! natural normal points away from the axis, which is the sphere's
//! convention. Its meridians are the iso-curves u = 0 and u = 1.

use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::{Point3, Vec3};

use crate::cap::CapEndpoint;

/// Meridian profile control points, as (axial, radial) multiples of r.
const PROFILE: [(f64, f64); 5] = [(0.0, 1.0), (0.5, 1.0), (1.0, 1.0), (1.0, 0.5), (1.0, 0.0)];

/// Arc spans below this (radians) or within this of π are degenerate
/// section arcs: the contacts coincide or sit antipodal on the ball.
const MIN_ARC_SPAN: f64 = 1e-9;

/// The tangent-match cap surface for `endpoint`, or `None` when the
/// endpoint is degenerate: zero radius or tangent, contacts on the
/// axis, or a section arc of zero or half-turn span.
pub(crate) fn tangent_match_surface(endpoint: &CapEndpoint) -> Option<NurbsSurface> {
    let t = endpoint.spine_tangent_outward.normalized()?;
    let r = endpoint.radius;
    if !(r.is_finite() && r > 0.0) {
        return None;
    }
    let radial = |c: Point3| {
        let d = c - endpoint.center;
        (d - t * d.dot(&t)).normalized()
    };
    let (mut e_a, mut e_b) = (radial(endpoint.contact_1)?, radial(endpoint.contact_2)?);
    // Sweep counter-clockwise about +t, so that ∂u × ∂v points away
    // from the axis.
    if e_a.cross(&e_b).dot(&t) < 0.0 {
        std::mem::swap(&mut e_a, &mut e_b);
    }
    let span = e_a.dot(&e_b).clamp(-1.0, 1.0).acos();
    if !(MIN_ARC_SPAN..std::f64::consts::PI - MIN_ARC_SPAN).contains(&span) {
        return None;
    }
    let (directions, weights, knots_u) = unit_arc(e_a, t, span);
    let num_v = PROFILE.len();
    let mut control_points = Vec::with_capacity(directions.len() * num_v);
    let mut all_weights = Vec::with_capacity(directions.len() * num_v);
    for (direction, weight) in directions.iter().zip(&weights) {
        for (axial, radius) in PROFILE {
            control_points.push(endpoint.center + t * (axial * r) + *direction * (radius * r));
            all_weights.push(*weight);
        }
    }
    Some(NurbsSurface {
        degree_u: 2,
        degree_v: PROFILE.len() - 1,
        knots_u,
        knots_v: vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        num_u: directions.len(),
        num_v,
        control_points,
        weights: all_weights,
    })
}

/// The rational quadratic unit arc from `e_a`, `span` radians
/// counter-clockwise about `axis`: its control directions, weights and
/// knots. One segment up to 90°, two beyond.
fn unit_arc(e_a: Vec3, axis: Vec3, span: f64) -> (Vec<Vec3>, Vec<f64>, Vec<f64>) {
    let rotate = |angle: f64| e_a * angle.cos() + axis.cross(&e_a) * angle.sin();
    let segments = if span > std::f64::consts::FRAC_PI_2 {
        2
    } else {
        1
    };
    let step = span / segments as f64;
    let mut directions = vec![rotate(0.0)];
    let mut weights = vec![1.0];
    for k in 0..segments {
        let (d0, d1) = (rotate(k as f64 * step), rotate((k + 1) as f64 * step));
        directions.push((d0 + d1) * (1.0 / (1.0 + step.cos())));
        weights.push((0.5 * step).cos());
        directions.push(d1);
        weights.push(1.0);
    }
    let knots = if segments == 1 {
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    } else {
        vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0]
    };
    (directions, weights, knots)
}

/// The meridian iso-curve of `surface` at `u = 0` (`at_end == false`)
/// or `u = 1`, from the section (t = 0) to the pole (t = 1). The end
/// rows of a rational arc have unit weight, so the meridian is the
/// polynomial quartic through that row's control points.
pub(crate) fn meridian(surface: &NurbsSurface, at_end: bool) -> NurbsCurve {
    let row = if at_end { surface.num_u - 1 } else { 0 };
    let start = row * surface.num_v;
    NurbsCurve {
        degree: surface.degree_v,
        knots: surface.knots_v.clone(),
        control_points: surface.control_points[start..start + surface.num_v].to_vec(),
        weights: vec![1.0; surface.num_v],
    }
}

#[cfg(test)]
#[path = "cap_tangent_match_tests.rs"]
mod tests;
