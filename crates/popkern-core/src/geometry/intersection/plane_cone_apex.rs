//! Plane through a cone's apex: the intersection is the cone's
//! generator lines in the plane (PK-931).
//!
//! The degenerate case used to fall through to the marching solver,
//! which traced each generator as a fitted curve. On a cone whose chart
//! seam (u = 0) lies on one of those generators the fit's samples
//! straddled the 0 / 2π branches and its uv curve crossed the whole
//! chart (measured: a cone minus a box whose corner edge runs down the
//! cone axis, rotated about y by −0.4 — the generator record's uv
//! midpoint sat at u 5.32, 0.32 off its curve in 3D, and the cut
//! quadrant inside the box sampled the cone's far side and was kept).
//! The generators are exact lines: u constant on the cone, a straight
//! segment on the plane.
//!
//! Cone chart: `apex + axis·v + (p1·cos u + p2·sin u)·v·tan(α)`, v ≥ 0.
//! A generator at azimuth u lies in the plane when
//! `n·axis + tan(α)·(n·p1·cos u + n·p2·sin u) = 0`, i.e.
//! `amp·cos(u − base) = −(n·axis)/tan(α)`.

use popkern_math::PortableMath;

use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

/// The generator lines of the cone (`apex`, axis `ca`, `half_angle`)
/// in the plane through its apex (`po`, normal `pn`, axes `pu`/`pv`),
/// spanning the cone's v window `cvr` with headroom from the apex. A
/// plane meeting the cone only at its apex, or tangent along one
/// generator, yields nothing (as plane × cylinder treats tangency).
#[allow(clippy::too_many_arguments)]
pub(super) fn case_through_apex(
    (po, pn, pu, pv): (Point3, Vec3, Vec3, Vec3),
    apex: Point3,
    ca: Vec3,
    half_angle: f64,
    cvr: (f64, f64),
    tol: &Tolerance,
    swapped: bool,
) -> Vec<SurfaceSurfaceIntersection> {
    let slope = half_angle.ptan();
    let (p1, p2) = orthonormal_basis(ca);
    let (n1, n2) = (pn.dot(&p1), pn.dot(&p2));
    let amp = (n1 * n1 + n2 * n2).sqrt() * slope;
    if amp < tol.angular() {
        return Vec::new(); // perpendicular to the axis: the apex alone
    }
    let ratio = -pn.dot(&ca) / amp;
    if ratio.abs() > 1.0 - tol.angular() {
        return Vec::new(); // the apex alone, or tangent along a generator
    }
    let base = n2.patan2(n1);
    let delta = ratio.pacos();
    let tau = std::f64::consts::TAU;
    let span = if cvr.1.is_finite() {
        cvr.1.max(0.0)
    } else {
        0.0
    };
    let v_hi = span + span.max(1.0);
    [base + delta, base - delta]
        .iter()
        .map(|&a| {
            let u = a.rem_euclid(tau);
            let dir = ca + (p1 * u.pcos() + p2 * u.psin()) * slope;
            let tip = apex + dir * v_hi;
            let uv_cone = line(Point3::new(u, 0.0, 0.0), Point3::new(u, v_hi, 0.0));
            let uv_plane = line(plane_uv(apex, po, pu, pv), plane_uv(tip, po, pu, pv));
            let (curve_uv1, curve_uv2) = if swapped {
                (uv_cone, uv_plane)
            } else {
                (uv_plane, uv_cone)
            };
            SurfaceSurfaceIntersection::Curve(CurveOnSurface {
                curve_3d: line(apex, tip),
                curve_uv1,
                curve_uv2,
                is_closed: false,
            })
        })
        .collect()
}

fn line(a: Point3, b: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    }
}

/// `pt` in the plane chart (origin `o`, axes `pu`, `pv`, possibly
/// non-orthogonal: solved through the axes' Gram matrix).
fn plane_uv(pt: Point3, o: Point3, pu: Vec3, pv: Vec3) -> Point3 {
    let d = pt - o;
    let (a, b, c) = (pu.dot(&pu), pu.dot(&pv), pv.dot(&pv));
    let (x, y) = (d.dot(&pu), d.dot(&pv));
    let det = a * c - b * b;
    if det.abs() < 1e-30 {
        return Point3::new(x, y, 0.0);
    }
    Point3::new((c * x - b * y) / det, (a * y - b * x) / det, 0.0)
}

#[cfg(test)]
#[path = "plane_cone_apex_tests.rs"]
mod tests;
