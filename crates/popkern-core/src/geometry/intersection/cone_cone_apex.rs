//! Two cones sharing their apex meet in rulings (PK-1127).
//!
//! Every ruling of a cone passes through its apex, so two cones with one
//! apex meet in the lines through it that are rulings of both: the unit
//! directions `d` at angle α₁ to the first axis and α₂ to the second (on
//! either nappe, `d · aᵢ = ±cos αᵢ`). Written `d = x a₁ + y a₂ + z n`
//! with `n = a₁ × a₂`, the two conditions fix `x` and `y` and the unit
//! length fixes `z²`: up to two lines per pair of nappe signs. The
//! marcher traced them as fitted curves of over a thousand control
//! points that bent at the apex, a chart singularity of both surfaces;
//! a ruling is a line of constant `u` on each cone's chart, through the
//! apex on both nappes alike, so each is exact.

use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, PortableMath, Vec3};
use crate::tolerance::Tolerance;

/// A parameter window `((u0, u1), (v0, v1))`.
type Window = ((f64, f64), (f64, f64));

/// The rulings two cones with one apex share, clipped to both windows;
/// `None` when the pair is not two cones with one apex and distinct axes.
pub(crate) fn try_cone_cone_apex(
    s1: &SurfaceGeometry,
    w1: Window,
    s2: &SurfaceGeometry,
    w2: Window,
    tolerance: &Tolerance,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let (
        SurfaceGeometry::Cone {
            apex: p1,
            axis: a1,
            half_angle: h1,
        },
        SurfaceGeometry::Cone {
            apex: p2,
            axis: a2,
            half_angle: h2,
        },
    ) = (s1, s2)
    else {
        return None;
    };
    if p1.distance_to(p2) > tolerance.linear() {
        return None;
    }
    let (a1, a2) = (a1.normalized()?, a2.normalized()?);
    let n = a1.cross(&a2);
    let s = n.length();
    if s <= tolerance.angular() {
        return None;
    }
    let c = a1.dot(&a2);
    let mut out = Vec::new();
    // One line per (s1, s2) up to its reversal: (+, ±) covers all four.
    for sign2 in [1.0, -1.0] {
        let (c1, c2) = (h1.pcos(), sign2 * h2.pcos());
        // [1 c; c 1] (x, y) = (c1, c2).
        let det = 1.0 - c * c;
        let x = (c1 - c * c2) / det;
        let y = (c2 - c * c1) / det;
        let in_plane = a1 * x + a2 * y;
        let z_sq = (1.0 - in_plane.length() * in_plane.length()) / (s * s);
        if z_sq < -1e-12 {
            continue;
        }
        let z = z_sq.max(0.0).sqrt();
        let zs: &[f64] = if z <= 1e-12 { &[0.0] } else { &[1.0, -1.0] };
        for &zsign in zs {
            let d = in_plane + n * (zsign * z);
            if let Some(cos) = ruling(*p1, d, (s1, w1, a1), (s2, w2, a2), tolerance) {
                out.push(cos);
            }
        }
    }
    Some(out)
}

/// The ruling through `apex` along `d`, clipped to both cones' v
/// windows, as an exact line on both charts.
fn ruling(
    apex: Point3,
    d: Vec3,
    (s1, w1, a1): (&SurfaceGeometry, Window, Vec3),
    (s2, w2, a2): (&SurfaceGeometry, Window, Vec3),
    tolerance: &Tolerance,
) -> Option<SurfaceSurfaceIntersection> {
    let d = d.normalized()?;
    // t along d gives v = t (d · a) on each cone.
    let range = |w: Window, a: Vec3| -> Option<(f64, f64)> {
        let k = d.dot(&a);
        if k.abs() <= 1e-12 {
            return None;
        }
        let (v0, v1) = (w.1.0.min(w.1.1), w.1.0.max(w.1.1));
        let (t0, t1) = (v0 / k, v1 / k);
        Some((t0.min(t1), t0.max(t1)))
    };
    let (r1, r2) = (range(w1, a1)?, range(w2, a2)?);
    let (t0, t1) = (r1.0.max(r2.0), r1.1.min(r2.1));
    if t1 - t0 <= tolerance.linear() {
        return None;
    }
    let line = |a: Point3, b: Point3| NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    };
    let curve_3d = line(apex + d * t0, apex + d * t1);
    let uv = |s: &SurfaceGeometry, a: Vec3| -> Option<NurbsCurve> {
        let u = ruling_u(s, d, a)?;
        Some(line(
            Point3::new(u, t0 * d.dot(&a), 0.0),
            Point3::new(u, t1 * d.dot(&a), 0.0),
        ))
    };
    Some(SurfaceSurfaceIntersection::Curve(CurveOnSurface {
        curve_3d,
        curve_uv1: uv(s1, a1)?,
        curve_uv2: uv(s2, a2)?,
        is_closed: false,
    }))
}

/// The chart angle `u` of the ruling along `d` (the same on both
/// nappes: the ruling at `u` passes through the apex).
fn ruling_u(s: &SurfaceGeometry, d: Vec3, a: Vec3) -> Option<f64> {
    use crate::geometry::analytic_surfaces::orthonormal_basis;
    let v = d.dot(&a);
    let radial = d - a * v;
    let (p1, p2) = orthonormal_basis(a);
    // point = apex + a v + v tanα (cos u p1 + sin u p2): the radial part
    // carries the sign of v.
    let sign = v.signum();
    let u = (sign * radial.dot(&p2)).patan2(sign * radial.dot(&p1));
    let _ = s;
    Some(u.rem_euclid(std::f64::consts::TAU))
}

#[cfg(test)]
#[path = "cone_cone_apex_tests.rs"]
mod tests;
