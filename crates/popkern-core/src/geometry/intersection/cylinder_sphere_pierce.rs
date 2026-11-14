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
//! linear tolerance. A branch through a pole of the sphere's chart is
//! cut there, each piece ending at the pole at its own limit longitude
//! (`split_at_poles`); a loop sweeping the longitude near a pole (a rod
//! beside the pole axis) declines to the marching path.

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
/// The largest step in sphere-chart longitude the refinement beside a
/// pole leaves (radians): the fit's cubic through coarser steps
/// overshoots the longitude between them, 8e-5 off the ball at points
/// 1e-3 from the axis under quarter-radian steps.
const REFINED_SWEEP: f64 = 0.005;
/// The longitude step of a scan gap (radians) that marks it a pole's
/// sweep, to be refined.
const POLE_SWEEP_GAP: f64 = 0.0625;
/// Bisections of one gap between samples beside a pole.
const MAX_POLE_DEPTH: usize = 48;

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
    // The sphere's poles lying on the cylinder: a branch through one
    // is split there (`split_at_poles`).
    let poles: Vec<Point3> = [1.0, -1.0]
        .into_iter()
        .map(|s| center + pole_axis * (s * rs))
        .filter(|p| {
            let rel = *p - origin;
            let off = (rel - a * rel.dot(&a)).length();
            (off - rc).abs() <= tolerance.linear()
        })
        .collect();
    let pad = (cyl_v.1 - cyl_v.0).abs() * 0.5 + tolerance.linear();
    let mut density = LOOP_SAMPLES;
    loop {
        let (loops, open) =
            super::ruling_pierce::ruling_loops(&quad, &discs, 1.0, &crossings, density);
        let mut worst = 0.0_f64;
        let mut results = Vec::new();
        let (p1, p2) = orthonormal_basis(pole_axis);
        let longitude = |p: Point3| {
            let d = p - center;
            d.dot(&p2).patan2(d.dot(&p1))
        };
        let floor = 1e-3 * tolerance.linear();
        let off_axis = |p: Point3| {
            let d = p - center;
            (d - pole_axis * d.dot(&pole_axis)).length()
        };
        let root_near = |theta: f64, near: f64| {
            let (disc, b, _) = quad(theta);
            let s = disc.max(0.0).sqrt();
            let (lo, hi) = (0.5 * (-b - s), 0.5 * (-b + s));
            if (lo - near).abs() <= (hi - near).abs() {
                lo
            } else {
                hi
            }
        };
        let pieces = loops.iter().flat_map(|pts| {
            let before = pts.len();
            let pts = refine_pole_sweeps(
                pts,
                &|theta, t| foot(theta) + a * t,
                &longitude,
                &root_near,
                (floor, tolerance.linear(), &off_axis),
            );
            let refined = pts.len() > before;
            let points: Vec<Point3> = pts.iter().map(|&(theta, t)| foot(theta) + a * t).collect();
            split_at_poles(points, !open, &poles)
                .into_iter()
                .map(move |(piece, closed)| (piece, closed, refined))
        });
        for (points, closed, refined) in pieces {
            if points.len() < 4 {
                continue;
            }
            let mut samples = Vec::with_capacity(points.len());
            let (mut uc, mut us) = (None, None);
            // The branch's tangent at a pole: across both normals,
            // the sphere's the pole axis and the cylinder's radial.
            let tangent_at = |p: Point3| -> Option<Vec3> {
                let rel = p - origin;
                let radial = (rel - a * rel.dot(&a)).normalized()?;
                let t = pole_axis.cross(&radial);
                (t.length() > 1e-3).then_some(t)
            };
            let longitudes =
                pole_safe_longitudes(&points, (center, pole_axis), sphere_chart, tangent_at)?;
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
            // turn, its sphere pcurve closing a period on as a latitude
            // circle's does. Passing close enough to the pole that its
            // steps there were refined (`refine_pole_sweeps`), it is
            // fitted here: left to the marching path, a rod whose wall
            // passed 2e-6 beside a ball's pole lost half its loop as a
            // duplicate record and the union served the rod alone
            // (PK-1176). Farther off, the marching path's records (cut
            // at the chart's seam) keep a rod holding the poles at one
            // volume through every turn to 3e-9, where the fitted ring,
            // cut at the face's seam and joined again, read 2.7e-6 off.
            let sphere_u = |s: &MarchSample| if swapped { s.u1 } else { s.u2 };
            if !refined
                && let (Some(first), Some(last)) = (samples.first(), samples.last())
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
            match fit_march_samples(&samples, closed) {
                Ok(curve) => {
                    // The fit's first pcurve is the caller's first
                    // surface's: the sphere's when `swapped`.
                    let (first, second) = if swapped {
                        (sphere_surface, cyl_surface)
                    } else {
                        (cyl_surface, sphere_surface)
                    };
                    worst = worst.max(super::ruling_pierce::fit_deviation(
                        &curve,
                        first,
                        second,
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

/// For each point, its longitude when it lies on a pole of the sphere;
/// `None` elsewhere. At a pole the chart's longitude is `atan2(0, 0)`,
/// arbitrary: a branch ending there (Viviani's curve with the ball's
/// pole at its double point) jumped by it, read as a sweep round the
/// pole and declined to the marcher, which stopped 7.5e-5 short of the
/// point and served pcurves 3.1e-4 off (PK-1139). Approaching a pole
/// along a branch the longitude has a limit, the direction the branch
/// leaves the pole in: `tangent_at` the pole, turned toward the nearest
/// point off it, where the surfaces cross there; that point's own
/// longitude where they are tangent (no tangent across the normals).
/// The neighbour's longitude is a step's turn off the limit, a kink in
/// the pcurve its fit converged to at second order (PK-1169: 1.7e-7 at
/// 16384 samples a turn). "On the pole" is relative to the branch's
/// own motion: within a thousandth of the step to the next sample of
/// the axis (a crossing refined by golden section lies 1e-7 off it,
/// its longitude noise).
fn pole_safe_longitudes(
    points: &[Point3],
    axis: (Point3, Vec3),
    chart: impl Fn(Point3, Option<f64>) -> Option<(f64, f64)>,
    tangent_at: impl Fn(Point3) -> Option<Vec3>,
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
        let toward = points[neighbour] - points[k];
        let along = match tangent_at(points[k]) {
            Some(t) if t.dot(&toward) < 0.0 => points[k] - t * toward.length(),
            Some(t) => points[k] + t * toward.length(),
            None => points[neighbour],
        };
        out[k] = Some(chart(along, None)?.0);
    }
    Some(out)
}

/// The branch's `(θ, t)` samples with more between any two whose
/// sphere longitudes differ by more than `POLE_SWEEP`, while their
/// chord exceeds `floor`: a branch passing beside a pole of the ball
/// sweeps half a turn of longitude over a stretch as short as its
/// distance from the pole, which the scan's uniform θ steps jump, and
/// the scan declined it to the marching path, whose records a rod
/// 3e-4 or less beside a ball's pole failed to assemble (PK-1176).
/// Each sample added at the middle θ takes the ruling root nearer the
/// two it falls between (the branch's own). A sample within `tol` of
/// the pole axis is at the pole, its longitude the pole handling's
/// (`split_at_poles`, `pole_safe_longitudes`): a gap reaching one is
/// not refined.
fn refine_pole_sweeps(
    pts: &[(f64, f64)],
    at: &dyn Fn(f64, f64) -> Point3,
    longitude: &dyn Fn(Point3) -> f64,
    root_near: &dyn Fn(f64, f64) -> f64,
    (floor, tol, off_axis): (f64, f64, &dyn Fn(Point3) -> f64),
) -> Vec<(f64, f64)> {
    let step = |x: (f64, f64), y: (f64, f64)| {
        let (p, q) = (at(x.0, x.1), at(y.0, y.1));
        let d = longitude(q) - longitude(p);
        let d = d - TAU * (d / TAU).round();
        let open = p.distance_to(&q) > floor && off_axis(p) > tol && off_axis(q) > tol;
        (d.abs(), open)
    };
    let sweep = |x: (f64, f64), y: (f64, f64)| {
        let (d, open) = step(x, y);
        d > REFINED_SWEEP && open
    };
    let mut out = Vec::with_capacity(pts.len());
    for w in pts.windows(2) {
        out.push(w[0]);
        // Only a gap stepping a pole's sweep (a sixteenth of a radian or
        // more: a loop far from the poles steps a hundredth at 512 a
        // turn) is refined, down to `REFINED_SWEEP`.
        if !matches!(step(w[0], w[1]), (d, true) if d > POLE_SWEEP_GAP) {
            continue;
        }
        // Depth-first bisection of this gap, emitted in order.
        let mut stack = vec![(w[0], w[1], 0usize)];
        let mut pending: Vec<(f64, f64)> = Vec::new();
        while let Some((x, y, depth)) = stack.pop() {
            if depth >= MAX_POLE_DEPTH || !sweep(x, y) {
                pending.push(y);
                continue;
            }
            let theta = 0.5 * (x.0 + y.0);
            let m = (theta, root_near(theta, 0.5 * (x.1 + y.1)));
            stack.push((m, y, depth + 1));
            stack.push((x, m, depth + 1));
        }
        // The gap's end is pushed last of all; the samples before it
        // came in order.
        pending.pop();
        out.extend(pending);
    }
    if let Some(&last) = pts.last() {
        out.push(last);
    }
    out
}

/// The branch `points` (closed when `closed`) cut at each of `poles` it
/// passes through, the pole itself ending both pieces: a branch through
/// a pole of the sphere's chart leaves it at the longitude half a turn
/// from the one it arrived at, a jump no single pcurve makes, and the
/// scan declined it to the marching path as a sweep round the pole,
/// which traced Viviani's curve with the ball's poles on it as nine
/// records, near-duplicate arcs 2e-3 apart and loops at the double
/// point among them (PK-1169). Each piece's end at the pole takes the
/// longitude of its own side (`pole_safe_longitudes`). A branch passes
/// through a pole where the pole lies within a sample step of a chord
/// of the branch; the sample nearest it within a thousandth of the step
/// is replaced, otherwise the pole is put in between.
fn split_at_poles(
    mut points: Vec<Point3>,
    closed: bool,
    poles: &[Point3],
) -> Vec<(Vec<Point3>, bool)> {
    let mut cuts: Vec<usize> = Vec::new();
    for pole in poles {
        let n = points.len();
        if n < 2 {
            break;
        }
        let Some((k, d)) = (0..n - 1)
            .map(|k| {
                let (p, q) = (points[k], points[k + 1]);
                let w = q - p;
                let s = if w.length_squared() > 0.0 {
                    ((*pole - p).dot(&w) / w.length_squared()).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                (k, pole.distance_to(&(p + w * s)))
            })
            .min_by(|x, y| x.1.total_cmp(&y.1))
        else {
            continue;
        };
        let step = points[k].distance_to(&points[k + 1]);
        if d > step {
            continue;
        }
        let at = if points[k].distance_to(pole) <= 1e-3 * step {
            points[k] = *pole;
            k
        } else if points[k + 1].distance_to(pole) <= 1e-3 * step {
            points[k + 1] = *pole;
            k + 1
        } else {
            points.insert(k + 1, *pole);
            for c in &mut cuts {
                if *c > k {
                    *c += 1;
                }
            }
            k + 1
        };
        cuts.push(at);
    }
    let n = points.len();
    // A cut at an end of an open branch, or at the closing sample of a
    // closed one, cuts nothing.
    cuts.retain(|&c| c > 0 && c + 1 < n);
    cuts.sort_unstable();
    cuts.dedup();
    if cuts.is_empty() {
        return vec![(points, closed)];
    }
    if closed {
        // The closed branch's samples run round once, the last
        // repeating the first: start it at the first cut.
        let first = cuts[0];
        let mut ring: Vec<Point3> = points[first..n - 1].to_vec();
        ring.extend_from_slice(&points[..=first]);
        let shifted: Vec<usize> = cuts.iter().map(|&c| c - first).collect();
        return pieces_between(&ring, &shifted[1..]);
    }
    pieces_between(&points, &cuts)
}

/// The open pieces of `points` between consecutive `cuts` (interior
/// indices, ascending), each cut sample ending one piece and starting
/// the next.
fn pieces_between(points: &[Point3], cuts: &[usize]) -> Vec<(Vec<Point3>, bool)> {
    let mut out = Vec::new();
    let mut start = 0;
    for &c in cuts.iter().chain(std::iter::once(&(points.len() - 1))) {
        if c > start {
            out.push((points[start..=c].to_vec(), false));
        }
        start = c;
    }
    out
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
