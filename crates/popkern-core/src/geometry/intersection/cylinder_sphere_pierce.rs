//! Cylinder×sphere with the sphere's centre off the cylinder's axis:
//! the ruling scan (`ruling_pierce`). Along the ruling at angle θ,
//! `o(θ) + t·a` (unit `a`), the sphere is met where
//! `t² + 2(a·d)t + |d|² − R² = 0`, `d = o(θ) − c`; sweeping θ yields
//! the exact branch structure. The marching path traced a ball bored
//! off its centre as one loop distorted and the other collapsed to a
//! point, and the boolean served the ball unchanged (brep-bench
//! boolean/offset-hole-sphere: volume 33510 for 29241).
//!
//! The samples are exact surface points; the fit through them is
//! densified until it lies on both surfaces within a tenth of the
//! linear tolerance. A loop sweeping the sphere chart's longitude near
//! a pole (a rod beside the pole axis) declines to the marching path.

use popkern_math::PortableMath;
use std::f64::consts::TAU;

use super::cylinder_cylinder_singular::singular_crossings;
use crate::error::CoreError;
use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::intersection::marching::MarchSample;
use crate::geometry::intersection::surface_surface::{
    SurfaceSurfaceIntersection, fit_march_samples,
};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::uv_projection::project_point_to_uv;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

const LOOP_SAMPLES: usize = 512;
const MAX_LOOP_SAMPLES: usize = 16_384;
const FIT_BUDGET: f64 = 0.1;
/// The largest step in sphere-chart longitude between consecutive
/// samples the scan fits (radians); a larger one is a pole sweep.
const POLE_SWEEP: f64 = 0.25;

/// The pierce of the cylinder `(origin, axis, rc)` and the sphere
/// `(center, rs)`; `swapped` when the caller's first surface is the
/// sphere. `cyl_v` is the cylinder's requested height window.
#[allow(clippy::too_many_arguments)]
pub(super) fn pierce(
    origin: Point3,
    axis: Vec3,
    rc: f64,
    cyl_v: (f64, f64),
    center: Point3,
    rs: f64,
    (cyl_surface, sphere_surface): (&SurfaceGeometry, &SurfaceGeometry),
    tolerance: &Tolerance,
    swapped: bool,
) -> Option<Vec<SurfaceSurfaceIntersection>> {
    let a = axis.normalized()?;
    let SurfaceGeometry::Sphere { axis: pole, .. } = sphere_surface else {
        return None;
    };
    let pole_axis = pole.normalized()?;
    let (e1, e2) = orthonormal_basis(a);
    let (rc, rs) = (rc.abs(), rs.abs());
    let foot = |theta: f64| origin + e1 * (rc * theta.pcos()) + e2 * (rc * theta.psin());
    // t² + b t + c with b = 2 a·d, c = |d|² − R².
    let quad = |theta: f64| -> (f64, f64, f64) {
        let d = foot(theta) - center;
        let b = 2.0 * a.dot(&d);
        let c = d.dot(&d) - rs * rs;
        (b * b - 4.0 * c, b, c)
    };
    let discs: Vec<f64> = (0..LOOP_SAMPLES)
        .map(|i| quad(i as f64 / LOOP_SAMPLES as f64 * TAU).0)
        .collect();
    if discs.iter().all(|&d| d <= 0.0) {
        return Some(Vec::new());
    }
    let crossings = singular_crossings(&quad, &discs, 1.0, tolerance.linear());
    // The axial position `t` of a ruling point is the cylinder's chart
    // v up to the origin's own offset (v = (p − origin)·axis / |axis|²).
    let cyl_chart = |p: Point3, prev: Option<f64>| -> (f64, f64) {
        let rel = p - origin;
        let v = rel.dot(&axis) / axis.length_squared();
        let u = unwrap(rel.dot(&e2).patan2(rel.dot(&e1)), prev);
        (u, v)
    };
    let sphere_chart = |p: Point3, prev: Option<f64>| -> Option<(f64, f64)> {
        let (u, v) = project_point_to_uv(sphere_surface, &p)?;
        Some((unwrap(u, prev), v))
    };
    let pad = (cyl_v.1 - cyl_v.0).abs() * 0.5 + tolerance.linear();
    let mut density = LOOP_SAMPLES;
    loop {
        let (loops, open) =
            super::ruling_pierce::ruling_loops(&quad, &discs, 1.0, &crossings, density);
        let mut worst = 0.0_f64;
        let mut results = Vec::new();
        for pts in loops {
            if pts.len() < 4 {
                continue;
            }
            let mut samples = Vec::with_capacity(pts.len());
            let (mut uc, mut us) = (None, None);
            let points: Vec<Point3> = pts.iter().map(|&(theta, t)| foot(theta) + a * t).collect();
            let longitudes = pole_safe_longitudes(&points, (center, pole_axis), sphere_chart)?;
            for (k, &p) in points.iter().enumerate() {
                let (cu, cv) = cyl_chart(p, uc);
                let (raw_u, sv) = sphere_chart(p, us)?;
                let su = longitudes[k].map_or(raw_u, |u| unwrap(u, us));
                // Near the sphere chart's pole a loop sweeps longitude
                // faster than θ steps: the fit there is not this
                // scan's to make (a rod beside a sphere's poles lost
                // its faces' areas). Decline to the marching path,
                // which owns pole sweeps (PK-1055).
                if us.is_some_and(|prev: f64| (su - prev).abs() > POLE_SWEEP) {
                    return None;
                }
                uc = Some(cu);
                us = Some(su);
                let ((u1, v1), (u2, v2)) = if swapped {
                    ((su, sv), (cu, cv))
                } else {
                    ((cu, cv), (su, sv))
                };
                samples.push(MarchSample {
                    u1,
                    v1,
                    u2,
                    v2,
                    point: p,
                });
            }
            // A loop round the chart's pole winds its longitude a whole
            // turn: the marching path's too (PK-1055 traces it).
            let sphere_u = |s: &MarchSample| if swapped { s.u1 } else { s.u2 };
            if let (Some(first), Some(last)) = (samples.first(), samples.last())
                && (sphere_u(last) - sphere_u(first)).abs() > std::f64::consts::PI
            {
                return None;
            }
            let cyl_vs = samples.iter().map(|s| if swapped { s.v2 } else { s.v1 });
            if !cyl_vs
                .into_iter()
                .any(|v| v > cyl_v.0 - pad && v < cyl_v.1 + pad)
            {
                continue;
            }
            match fit_march_samples(&samples, !open) {
                Ok(curve) => {
                    worst = worst.max(super::ruling_pierce::fit_deviation(
                        &curve,
                        cyl_surface,
                        sphere_surface,
                        samples.len(),
                    ));
                    results.push(SurfaceSurfaceIntersection::Curve(curve));
                }
                Err(CoreError::PreconditionFailed(_)) | Err(_) => continue,
            }
        }
        if worst <= FIT_BUDGET * tolerance.linear() || density >= MAX_LOOP_SAMPLES {
            return Some(results);
        }
        density *= 2;
    }
}

/// For each point, its longitude when it lies on a pole of the sphere,
/// taken from the nearest point off it; `None` elsewhere. At a pole
/// the chart's longitude is `atan2(0, 0)`, arbitrary: a branch ending
/// there (Viviani's curve with the ball's pole at its double point)
/// jumped by it, read as a sweep round the pole and declined to the
/// marcher, which stopped 7.5e-5 short of the point and served
/// pcurves 3.1e-4 off (PK-1139). Approaching a pole along a branch the
/// longitude has a limit, its neighbour's. "On the pole" is relative
/// to the branch's own motion: within a thousandth of the step to the
/// next sample of the axis (a crossing refined by golden section lies
/// 1e-7 off it, its longitude noise).
fn pole_safe_longitudes(
    points: &[Point3],
    axis: (Point3, Vec3),
    chart: impl Fn(Point3, Option<f64>) -> Option<(f64, f64)>,
) -> Option<Vec<Option<f64>>> {
    let (centre, pole) = axis;
    let off_axis = |p: &Point3| {
        let d = *p - centre;
        (d - pole * d.dot(&pole)).length()
    };
    let step = |k: usize| {
        let next = if k + 1 < points.len() {
            k + 1
        } else {
            k.saturating_sub(1)
        };
        points[k].distance_to(&points[next])
    };
    let on_pole: Vec<bool> = (0..points.len())
        .map(|k| off_axis(&points[k]) <= 1e-3 * step(k))
        .collect();
    if on_pole.iter().all(|&b| b) {
        return None;
    }
    let mut out = vec![None; points.len()];
    for k in (0..points.len()).filter(|&k| on_pole[k]) {
        let neighbour = (1..points.len()).find_map(|d| {
            [k.checked_sub(d), Some(k + d)]
                .into_iter()
                .flatten()
                .find(|&j| j < points.len() && !on_pole[j])
        })?;
        out[k] = Some(chart(points[neighbour], None)?.0);
    }
    Some(out)
}

/// `u` on the branch nearest `prev` (charts are periodic in u).
fn unwrap(mut u: f64, prev: Option<f64>) -> f64 {
    if let Some(pu) = prev {
        while u - pu > std::f64::consts::PI {
            u -= TAU;
        }
        while pu - u > std::f64::consts::PI {
            u += TAU;
        }
    }
    u
}

#[cfg(test)]
#[path = "cylinder_sphere_pierce_tests.rs"]
mod tests;
