//! Local NURBS / projection helpers used by the analytic dispatch.
//! Kept private to the `edge_triple` module — no public surface.

use crate::geometry::curve::CurveEval;
use crate::geometry::curve_on_surface::CurveOnSurface;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, EdgeHint};
use crate::geometry::edge_triple::options::{BranchPolicy, EdgeTripleOptions};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};

// ─── Branch disambiguation (ADR §4) ───────────────────────────────

pub(super) fn pick_branch(
    curves: &[&CurveOnSurface],
    anchors: &EdgeAnchors,
    hint: Option<&EdgeHint>,
    opts: &EdgeTripleOptions,
) -> Result<usize, EdgeTripleError> {
    if curves.len() == 1 {
        return Ok(0);
    }
    if let Some(EdgeHint::BranchIndex(i)) = hint {
        let idx = *i as usize;
        if idx >= curves.len() {
            return Err(EdgeTripleError::AmbiguousBranch {
                count: curves.len(),
            });
        }
        return Ok(idx);
    }
    if let Some(EdgeHint::Curve3D(c)) = hint {
        let mid = curve_midpoint(c);
        return Ok(min_by_key(curves, |b| {
            curve_midpoint(&b.curve_3d).distance_to(&mid)
        }));
    }
    match opts.prefer_branches {
        BranchPolicy::First => Ok(0),
        BranchPolicy::Reject => Err(EdgeTripleError::AmbiguousBranch {
            count: curves.len(),
        }),
        BranchPolicy::ClosestToAnchors => {
            let target = anchor_midpoint(anchors);
            let scored: Vec<(usize, f64)> = curves
                .iter()
                .enumerate()
                .map(|(i, c)| (i, curve_midpoint(&c.curve_3d).distance_to(&target)))
                .collect();
            let (i, best) = scored
                .iter()
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .copied()
                .unwrap();
            let _ = best;
            Ok(i)
        }
    }
}

pub(super) fn anchor_midpoint(anchors: &EdgeAnchors) -> Point3 {
    match anchors {
        EdgeAnchors::Open { start, end } => Point3::new(
            0.5 * (start.x + end.x),
            0.5 * (start.y + end.y),
            0.5 * (start.z + end.z),
        ),
        EdgeAnchors::Closed { sample, .. } | EdgeAnchors::Seam { sample, .. } => *sample,
    }
}

pub(super) fn curve_midpoint(c: &NurbsCurve) -> Point3 {
    let (t0, t1) = c.parameter_range();
    c.point_at(0.5 * (t0 + t1))
}

fn min_by_key<T, F>(slice: &[T], key: F) -> usize
where
    F: Fn(&T) -> f64,
{
    slice
        .iter()
        .enumerate()
        .map(|(i, t)| (i, key(t)))
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0)
}

// ─── NURBS slice / reverse ────────────────────────────────────────

pub(super) fn slice_curve(curve: &NurbsCurve, lo: f64, hi: f64) -> Option<NurbsCurve> {
    let (full_min, full_max) = curve.parameter_range();
    let eps = 1.0e-12 * (full_max - full_min).abs().max(1.0);
    if lo < full_min - eps || hi > full_max + eps {
        return None;
    }
    if hi <= lo + eps {
        return None;
    }
    // `split_at` requires strict interior. Clamp endpoints to the
    // curve's range so an anchor at exactly the start/end boundary
    // still trims successfully.
    let lo = lo.max(full_min);
    let hi = hi.min(full_max);
    let at_start = (lo - full_min).abs() <= eps;
    let at_end = (hi - full_max).abs() <= eps;

    let mut current = curve.clone();
    if !at_start {
        current = current.split_at(lo).ok()?.1;
    }
    if !at_end {
        let (cur_min, _) = current.parameter_range();
        let cut = if at_start { hi } else { cur_min + (hi - lo) };
        current = current.split_at(cut).ok()?.0;
    }
    Some(current)
}

pub(super) fn reverse_nurbs(curve: &NurbsCurve) -> NurbsCurve {
    let (t_min, t_max) = curve.parameter_range();
    let mut reversed_knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_min + t_max - k)
        .collect();
    reversed_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut reversed_cps = curve.control_points.clone();
    reversed_cps.reverse();
    let mut reversed_weights = curve.weights.clone();
    reversed_weights.reverse();
    NurbsCurve {
        degree: curve.degree,
        knots: reversed_knots,
        control_points: reversed_cps,
        weights: reversed_weights,
    }
}

// ─── Endpoint snapping (G1) ───────────────────────────────────────

pub(super) fn snap_endpoints_3d(curve: &mut NurbsCurve, start: Point3, end: Point3) {
    if let Some(first) = curve.control_points.first_mut() {
        *first = start;
    }
    if let Some(last) = curve.control_points.last_mut() {
        *last = end;
    }
}

pub(super) fn snap_endpoints_pcurve(
    pcurve: &mut NurbsCurve,
    surf: &SurfaceGeometry,
    start: Point3,
    end: Point3,
) {
    if let Some(first) = pcurve.control_points.first_mut() {
        let (u, v) = project_to_uv_clamped(surf, start);
        *first = Point3::new(u, v, 0.0);
    }
    if let Some(last) = pcurve.control_points.last_mut() {
        let (u, v) = project_to_uv_clamped(surf, end);
        *last = Point3::new(u, v, 0.0);
    }
}

/// Closed-form UV projection for analytic surfaces. Stage 2's
/// dispatch only reaches this for {Plane, Cylinder, Sphere, Cone};
/// later stages add Torus and NURBS via dedicated paths.
pub(super) fn project_to_uv_clamped(surf: &SurfaceGeometry, p: Point3) -> (f64, f64) {
    use SurfaceGeometry::*;
    let (u, v) = match surf {
        Plane {
            origin,
            u_axis,
            v_axis,
            ..
        } => {
            let d = p - *origin;
            (d.dot(u_axis), d.dot(v_axis))
        }
        Cylinder {
            origin,
            axis,
            radius,
        } => {
            let d = p - *origin;
            let v = d.dot(axis);
            let radial = d - (*axis) * v;
            let u_basis = pick_u_basis(axis);
            let v_basis = axis.cross(&u_basis);
            let u = radial.dot(&v_basis).atan2(radial.dot(&u_basis));
            let two_pi = std::f64::consts::TAU;
            let u = ((u % two_pi) + two_pi) % two_pi;
            let _ = radius;
            (u, v)
        }
        Sphere { center, radius } => {
            let d = p - *center;
            let v = (d.z / *radius).clamp(-1.0, 1.0).asin();
            let u = d.y.atan2(d.x);
            let two_pi = std::f64::consts::TAU;
            let u = ((u % two_pi) + two_pi) % two_pi;
            (u, v)
        }
        Cone {
            apex,
            axis,
            half_angle,
        } => {
            let d = p - *apex;
            let v = d.dot(axis);
            let u_basis = pick_u_basis(axis);
            let v_basis = axis.cross(&u_basis);
            let radial = d - (*axis) * v;
            let u = radial.dot(&v_basis).atan2(radial.dot(&u_basis));
            let two_pi = std::f64::consts::TAU;
            let u = ((u % two_pi) + two_pi) % two_pi;
            let _ = half_angle;
            (u, v)
        }
        // Out of scope for stage 2 — the classifier rejects these
        // before reaching this function.
        Torus { .. } | Nurbs(_) | Gregory(_) => (p.x, p.y),
    };
    let ((u_lo, u_hi), (v_lo, v_hi)) = surf.parameter_range();
    let clamp = |x: f64, lo: f64, hi: f64| {
        if !lo.is_finite() || !hi.is_finite() {
            x
        } else {
            x.clamp(lo, hi)
        }
    };
    (clamp(u, u_lo, u_hi), clamp(v, v_lo, v_hi))
}

fn pick_u_basis(axis: &Vec3) -> Vec3 {
    let candidate = if axis.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let n = candidate - (*axis) * candidate.dot(axis);
    n.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0))
}

// ─── pcurve rebuild for consistent parameterisation (G3) ──────────

/// Rebuild a [`CurveOnSurface`] so both pcurves share the 3-D
/// curve's degree, knots, and weights — projecting curve_3d's CPs
/// to each surface's UV space individually.
///
/// Necessary because some `intersect_surfaces` fast paths return
/// pcurves with different parameterisations than `curve_3d`
/// (e.g. `case_perpendicular` returns a degree-1 hline UV against
/// a degree-2 rational NURBS circle). The factory's G2/G3
/// guarantees require all three curves to share `(degree, knots,
/// weights)`.
pub(super) fn rebuild_consistent_cos(
    cos: &CurveOnSurface,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
) -> CurveOnSurface {
    let curve_3d = cos.curve_3d.clone();
    let curve_uv1 = pcurve_from_curve_3d(&curve_3d, surf_a);
    let curve_uv2 = pcurve_from_curve_3d(&curve_3d, surf_b);
    CurveOnSurface {
        curve_3d,
        curve_uv1,
        curve_uv2,
        is_closed: cos.is_closed,
    }
}

fn pcurve_from_curve_3d(curve: &NurbsCurve, surf: &SurfaceGeometry) -> NurbsCurve {
    let mut cps: Vec<Point3> = curve
        .control_points
        .iter()
        .map(|cp| {
            let (u, v) = project_to_uv_clamped(surf, *cp);
            Point3::new(u, v, 0.0)
        })
        .collect();
    // Unwrap u so the pcurve is continuous across any 2π / period seam.
    let two_pi = std::f64::consts::TAU;
    if matches!(
        surf,
        SurfaceGeometry::Cylinder { .. }
            | SurfaceGeometry::Sphere { .. }
            | SurfaceGeometry::Cone { .. }
    ) {
        for i in 1..cps.len() {
            let diff = cps[i].x - cps[i - 1].x;
            if diff > two_pi * 0.5 {
                cps[i].x -= two_pi;
            } else if diff < -two_pi * 0.5 {
                cps[i].x += two_pi;
            }
        }
    }
    NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points: cps,
        weights: curve.weights.clone(),
    }
}

// ─── Search-box derivation ────────────────────────────────────────

/// Compute a `(u_range, v_range)` search box for `surf`. For each
/// UV axis, the range is the surface's intrinsic range (for bounded
/// surfaces like Sphere / Torus / Gregory) or — for unbounded
/// directions (plane axes, cylinder v, cone v) — derived from the
/// anchor bounding box projected onto the surface's UV with a 10×
/// padding so the analytic intersection has headroom and avoids
/// catastrophic cancellation at huge parameter values.
pub(super) fn derive_search_box(
    surf: &SurfaceGeometry,
    anchors: &EdgeAnchors,
) -> ((f64, f64), (f64, f64)) {
    let ((u_lo, u_hi), (v_lo, v_hi)) = surf.parameter_range();
    let (anchor_u, anchor_v) = anchor_uv_extent(surf, anchors);
    let resolve = |intrinsic_lo: f64, intrinsic_hi: f64, derived: (f64, f64)| -> (f64, f64) {
        let lo = if intrinsic_lo.is_finite() {
            intrinsic_lo
        } else {
            derived.0
        };
        let hi = if intrinsic_hi.is_finite() {
            intrinsic_hi
        } else {
            derived.1
        };
        (lo, hi)
    };
    (resolve(u_lo, u_hi, anchor_u), resolve(v_lo, v_hi, anchor_v))
}

fn anchor_uv_extent(surf: &SurfaceGeometry, anchors: &EdgeAnchors) -> ((f64, f64), (f64, f64)) {
    let pts = anchor_points(anchors);
    let mut u_lo = f64::INFINITY;
    let mut u_hi = f64::NEG_INFINITY;
    let mut v_lo = f64::INFINITY;
    let mut v_hi = f64::NEG_INFINITY;
    for p in &pts {
        let (u, v) = project_to_uv_clamped(surf, *p);
        u_lo = u_lo.min(u);
        u_hi = u_hi.max(u);
        v_lo = v_lo.min(v);
        v_hi = v_hi.max(v);
    }
    let u_span = (u_hi - u_lo).abs().max(1.0);
    let v_span = (v_hi - v_lo).abs().max(1.0);
    let pad_u = 10.0 * u_span;
    let pad_v = 10.0 * v_span;
    ((u_lo - pad_u, u_hi + pad_u), (v_lo - pad_v, v_hi + pad_v))
}

fn anchor_points(anchors: &EdgeAnchors) -> Vec<Point3> {
    match anchors {
        EdgeAnchors::Open { start, end } => vec![*start, *end],
        EdgeAnchors::Closed { sample, .. } => vec![*sample],
        EdgeAnchors::Seam { sample, .. } => vec![*sample],
    }
}
