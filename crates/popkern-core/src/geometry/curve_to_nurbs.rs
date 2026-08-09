//! Analytic-to-NURBS conversion helpers for [`CurveGeometry`].
//!
//! Extracted from `curve.rs` to keep that file under the 500-line
//! structural cap. Builds degree-2 rational NURBS approximations of
//! analytic primitives (full circles, arcs, full ellipses, elliptical
//! arcs) using the standard 9-control-point construction with weight
//! `1/√2` on the alternating tangent-intersection points.

use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Vec3};

/// Parameter on a circle-derived rational NURBS whose POINT sits at
/// the requested absolute angle (measured in the `(u, v)` basis about
/// `center`).
///
/// A rational quadratic circle's parameter is NOT the angle: inside
/// each 90° knot span the two diverge by up to ~4% of the span, and
/// only agree at the span boundaries. Splitting such a curve at
/// `t = angle` therefore cuts at the wrong geometric point for any
/// bound that is not a multiple of π/2 (measured: a chamfer rim on
/// the notched flange asked to end at its notch-plane triple point
/// at angle 4.5025 was cut 0.0153 rad away — 0.37 of arc — and every
/// downstream stage inherited the shortfall). Callers must split at
/// the parameter this returns instead.
///
/// The unwrapped angle-from-start is monotone along the curve, so
/// bisection is exact and unconditionally convergent. Returns `None`
/// when the target lies at (or beyond) either end of the curve's
/// span — no split is needed there.
fn param_at_angle(
    curve: &NurbsCurve,
    center: Point3,
    u: Vec3,
    v: Vec3,
    target_angle: f64,
) -> Option<f64> {
    use crate::geometry::curve::CurveEval;
    use std::f64::consts::TAU;
    let (t_lo, t_hi) = curve.parameter_range();
    let raw_angle = |t: f64| -> f64 {
        let p = curve.point_at(t);
        let d = p - center;
        d.dot(&v).atan2(d.dot(&u))
    };
    let a_start = raw_angle(t_lo);
    // Angle relative to the curve start, unwrapped into [0, TAU).
    let rel = |t: f64| -> f64 {
        let mut a = raw_angle(t) - a_start;
        while a < -1e-9 {
            a += TAU;
        }
        while a >= TAU - 1e-9 {
            a -= TAU;
        }
        a.max(0.0)
    };
    let span = {
        // The end angle needs unwrapping too: a closed curve reports
        // rel(t_hi) = 0.
        let r = rel(t_hi);
        if r <= 1e-9 { TAU } else { r }
    };
    let mut target = target_angle - a_start;
    while target < 0.0 {
        target += TAU;
    }
    while target >= TAU {
        target -= TAU;
    }
    if target <= 1e-12 || target >= span - 1e-12 {
        return None;
    }
    let (mut lo, mut hi) = (t_lo, t_hi);
    for _ in 0..64 {
        let mid = 0.5 * (lo + hi);
        if rel(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

/// Build a full-circle NURBS (degree 2, 9 control points, 4 quarter arcs).
pub(super) fn circle_to_nurbs(center: Point3, normal: Vec3, radius: f64) -> NurbsCurve {
    use crate::geometry::analytic_curves::circle_basis;
    use std::f64::consts::{FRAC_PI_2, PI, TAU};

    let (u, v) = circle_basis(normal);
    let w = 1.0_f64 / 2.0_f64.sqrt(); // weight for midpoints

    let pts = [
        center + u * radius,        // angle 0
        center + (u + v) * radius,  // control pt (tangent intersection)
        center + v * radius,        // angle π/2
        center + (v - u) * radius,  // control pt
        center - u * radius,        // angle π
        center + (-u - v) * radius, // control pt
        center - v * radius,        // angle 3π/2
        center + (u - v) * radius,  // control pt
        center + u * radius,        // angle 2π
    ];

    let weights = vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
    let knots = vec![
        0.0,
        0.0,
        0.0,
        FRAC_PI_2,
        FRAC_PI_2,
        PI,
        PI,
        PI + FRAC_PI_2,
        PI + FRAC_PI_2,
        TAU,
        TAU,
        TAU,
    ];

    NurbsCurve {
        degree: 2,
        knots,
        control_points: pts.to_vec(),
        weights,
    }
}

/// Build a degree-2 rational NURBS for a full ellipse using
/// `(u_basis, v_basis)` as the two semi-axis vectors. Same 9-CP /
/// w=1/√2 structure as the circle, just with `radius·u_axis` and
/// `radius·v_axis` replaced by `u_basis` and `v_basis` directly.
pub(super) fn ellipse_full_to_nurbs(center: Point3, u_basis: Vec3, v_basis: Vec3) -> NurbsCurve {
    use std::f64::consts::{FRAC_PI_2, PI, TAU};
    let w = 1.0_f64 / 2.0_f64.sqrt();
    let pts = [
        center + u_basis,           // angle 0
        center + u_basis + v_basis, // tangent intersection at 0
        center + v_basis,           // angle π/2
        center + v_basis - u_basis, // tangent intersection
        center - u_basis,           // angle π
        center - u_basis - v_basis,
        center - v_basis, // angle 3π/2
        center + u_basis - v_basis,
        center + u_basis, // angle 2π
    ];
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0,
            0.0,
            0.0,
            FRAC_PI_2,
            FRAC_PI_2,
            PI,
            PI,
            PI + FRAC_PI_2,
            PI + FRAC_PI_2,
            TAU,
            TAU,
            TAU,
        ],
        control_points: pts.to_vec(),
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Build a NURBS elliptical arc by constructing the full ellipse
/// and splitting to the arc range.
///
/// Terminal-sliver guard (PK-257): when an endpoint lands within
/// atan2-noise of the full ellipse's quarter grid, the split strands
/// a sub-tolerance nano piece with a full-multiplicity break (the
/// PK-255 class — oblique fillet-end sections are ellipse arcs and
/// carried it exactly like the circular ones). Those arcs take the
/// uniform segmentation; the rational-quadratic shoulder construction
/// is exact on ellipses by affine invariance.
pub(super) fn ellipse_to_nurbs(
    center: Point3,
    u_basis: Vec3,
    v_basis: Vec3,
    start_angle: f64,
    end_angle: f64,
) -> NurbsCurve {
    use std::f64::consts::TAU;
    let quarter = TAU / 4.0;
    let collides = (1..4).any(|index| {
        let knot = index as f64 * quarter;
        knot > start_angle
            && knot < end_angle
            && ((knot - start_angle).abs() <= QUARTER_COLLISION_WINDOW
                || (end_angle - knot).abs() <= QUARTER_COLLISION_WINDOW)
    });
    if collides {
        return uniform_segmented_conic_arc(center, u_basis, v_basis, start_angle, end_angle);
    }
    let full = ellipse_full_to_nurbs(center, u_basis, v_basis);
    let right = if start_angle > 0.0 {
        match full.split_at(start_angle) {
            Ok((_l, r)) => r,
            Err(_) => full,
        }
    } else {
        full
    };
    if end_angle < std::f64::consts::TAU {
        match right.split_at(end_angle) {
            Ok((l, _r)) => l,
            Err(_) => right,
        }
    } else {
        right
    }
}

/// Uniformly segmented rational-quadratic conic arc: `n = ceil(span /
/// (π/2))` pieces of `span/n` each over exactly `[a0, a1]`, knots in
/// absolute angle, shoulder construction per piece (exact on circles
/// and — by affine invariance — ellipses). Sliver spans are
/// impossible by construction: every piece has width `span/n`.
pub(super) fn uniform_segmented_conic_arc(
    center: Point3,
    u_basis: Vec3,
    v_basis: Vec3,
    a0: f64,
    a1: f64,
) -> NurbsCurve {
    use std::f64::consts::TAU;
    let quarter = TAU / 4.0;
    let span = a1 - a0;
    let segments = ((span / quarter).ceil() as usize).max(1);
    let piece = span / segments as f64;
    let point_at_angle = |angle: f64, scale: f64| -> Point3 {
        center + u_basis * (scale * angle.cos()) + v_basis * (scale * angle.sin())
    };
    let shoulder_scale = {
        let half = 0.5 * piece;
        if half.cos().abs() > 1.0e-300 {
            1.0 / half.cos()
        } else {
            1.0
        }
    };
    let weight = (0.5 * piece).cos();
    let mut control_points = Vec::with_capacity(2 * segments + 1);
    let mut weights = Vec::with_capacity(2 * segments + 1);
    let mut knots = Vec::with_capacity(2 * segments + 4);
    knots.push(a0);
    for segment in 0..segments {
        let seg_start = a0 + piece * segment as f64;
        knots.push(seg_start);
        knots.push(seg_start);
        control_points.push(point_at_angle(seg_start, 1.0));
        weights.push(1.0);
        control_points.push(point_at_angle(seg_start + 0.5 * piece, shoulder_scale));
        weights.push(weight);
    }
    control_points.push(point_at_angle(a1, 1.0));
    weights.push(1.0);
    knots.push(a1);
    knots.push(a1);
    knots.push(a1);
    NurbsCurve {
        degree: 2,
        knots,
        control_points,
        weights,
    }
}

/// Build an exact degree-2 Bézier-equivalent NURBS for a parabolic
/// arc `vertex + t·linear + t²·quadratic` over `[t_start, t_end]`.
///
/// A parabola is exactly representable as a non-rational degree-2
/// Bézier curve. Endpoints are evaluations at `t_start` and `t_end`;
/// the middle control point sits at the intersection of the two
/// endpoint tangent lines. Equivalently, `b1 = b0 + tangent(t_start)
/// · (t_end − t_start) / 2`.
pub(super) fn parabola_to_nurbs(
    vertex: Point3,
    linear: Vec3,
    quadratic: Vec3,
    t_start: f64,
    t_end: f64,
) -> NurbsCurve {
    let span = t_end - t_start;
    let b0 = vertex + linear * t_start + quadratic * (t_start * t_start);
    let b2 = vertex + linear * t_end + quadratic * (t_end * t_end);
    let tangent_at_start = linear + quadratic * (2.0 * t_start);
    let b1 = b0 + tangent_at_start * (span * 0.5);
    NurbsCurve {
        degree: 2,
        knots: vec![t_start, t_start, t_start, t_end, t_end, t_end],
        control_points: vec![b0, b1, b2],
        weights: vec![1.0, 1.0, 1.0],
    }
}

/// Build an exact rational degree-2 NURBS approximation of a
/// hyperbolic arc `center + cosh(t)·semi_major + sinh(t)·semi_minor`
/// over `[t_start, t_end]`.
///
/// The hyperbola is exactly representable as a rational degree-2
/// Bézier via the substitution `u = tanh(t/2)`, which gives
/// `cosh(t) = (1+u²)/(1−u²)` and `sinh(t) = 2u/(1−u²)`. In
/// homogeneous form: `(b0 = center + sm_maj, b1 = center, b2 =
/// center + sm_maj)` doesn't quite work directly — the cleanest
/// closed form is to evaluate three sample points at the
/// endpoints and the parametric midpoint, then build the unique
/// rational quadratic Bézier through them. We fall back to a
/// sample-and-fit since the rational-quadratic construction needs
/// careful derivation; degree-5 with 33 samples gives ~1e-12
/// residual on typical sweeps.
pub(super) fn hyperbola_to_nurbs(
    center: Point3,
    semi_major: Vec3,
    semi_minor: Vec3,
    t_start: f64,
    t_end: f64,
) -> NurbsCurve {
    let n_samples = 33;
    let mut points: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let frac = i as f64 / (n_samples as f64 - 1.0);
        let t = t_start + frac * (t_end - t_start);
        let p = center + semi_major * t.cosh() + semi_minor * t.sinh();
        points.push(p);
        params.push(frac);
    }
    let nurbs = match crate::geometry::fitting::fit_curve(&points, &params, 5, n_samples) {
        Ok(c) => c,
        Err(_) => {
            // Degenerate: fall back to a degree-1 spline through endpoints.
            let Some(&first) = points.first() else {
                let origin = Point3::new(0.0, 0.0, 0.0);
                return NurbsCurve {
                    degree: 1,
                    knots: vec![t_start, t_start, t_end, t_end],
                    control_points: vec![origin, origin],
                    weights: vec![1.0, 1.0],
                };
            };
            let last = points.last().copied().unwrap_or(first);
            return NurbsCurve {
                degree: 1,
                knots: vec![t_start, t_start, t_end, t_end],
                control_points: vec![first, last],
                weights: vec![1.0, 1.0],
            };
        }
    };
    // Reknot the [0, 1] parameter range to [t_start, t_end] so
    // downstream consumers see the natural hyperbola parameter.
    let mut nurbs = nurbs;
    let k_lo = *nurbs.knots.first().unwrap_or(&0.0);
    let k_hi = *nurbs.knots.last().unwrap_or(&1.0);
    let span_old = k_hi - k_lo;
    let span_new = t_end - t_start;
    if span_old.abs() > 1.0e-30 {
        for k in &mut nurbs.knots {
            *k = t_start + (*k - k_lo) / span_old * span_new;
        }
    }
    nurbs
}

/// Build a NURBS arc covering the arc's OWN angular range, knots in
/// ABSOLUTE angle so the NURBS parameter matches the
/// [`CurveGeometry::Arc`] evaluator's parameter.
///
/// Arcs within the principal `[0, 2π]` range keep the historical
/// construction (splitting the standard full-circle NURBS) with its
/// parameterization byte-identical — the chamfer surgery pairs
/// exact-arc pcurves with these curves parametrically. Ranges
/// OUTSIDE the principal branch previously CLAMPED silently — a
/// closed rim authored as `Arc[−π/2, 3π/2]` (seam at −π/2) converted
/// to a curve covering only `[0, 3π/2]` with its start collapsed to
/// angle 0, and every consumer (edge splits, clipping, inversion)
/// inherited a rim missing a quarter. Out-of-principal ranges now
/// build on a basis ROTATED to the range start with knots shifted to
/// absolute angle.
///
/// [`CurveGeometry::Arc`]: crate::geometry::curve::CurveGeometry::Arc
pub(super) fn arc_to_nurbs(
    center: Point3,
    normal: Vec3,
    radius: f64,
    start_angle: f64,
    end_angle: f64,
) -> NurbsCurve {
    use std::f64::consts::TAU;

    let a0 = start_angle.min(end_angle);
    let mut a1 = start_angle.max(end_angle);
    if a1 - a0 > TAU {
        a1 = a0 + TAU;
    }
    if std::env::var_os("PK_ARC_TRACE").is_some() && !((0.0..=TAU).contains(&a0) && a1 <= TAU) {
        eprintln!("  [arc2nurbs] OUT-OF-PRINCIPAL range ({a0:.4},{a1:.4})");
    }

    if (0.0..=TAU).contains(&a0) && a1 <= TAU {
        // Principal range: historical construction, parameterization
        // preserved exactly on the retained span. Splits happen at
        // the PARAMETER whose point sits at the requested angle
        // (`param_at_angle`), not at `t = angle` — the rational
        // quadratic's parameter and the angle only agree at π/2
        // multiples.
        //
        // Terminal-sliver guard (PK-257): when an endpoint lands
        // within atan2-noise of the standard circle's quarter grid,
        // the historical split strands a sub-tolerance nano piece
        // whose full-multiplicity break poisons downstream knot
        // unions (the PK-255 class). Those arcs take the uniform
        // segmentation below — the historical byte-identical
        // parameterization is unreachable for them anyway (it IS the
        // sliver). Non-colliding arcs (every exactly-authored
        // chamfer/fillet arc) keep the historical construction.
        let quarter = TAU / 4.0;
        let collides = (1..4).any(|index| {
            let knot = index as f64 * quarter;
            knot > a0
                && knot < a1
                && ((knot - a0).abs() <= QUARTER_COLLISION_WINDOW
                    || (a1 - knot).abs() <= QUARTER_COLLISION_WINDOW)
        });
        if !collides {
            use crate::geometry::analytic_curves::circle_basis;
            let (u, v) = circle_basis(normal);
            let full = circle_to_nurbs(center, normal, radius);
            let right_of_start = match param_at_angle(&full, center, u, v, a0) {
                Some(t) => match full.split_at(t) {
                    Ok((_left, right)) => right,
                    Err(_) => full,
                },
                None => full,
            };
            return match param_at_angle(&right_of_start, center, u, v, a1) {
                Some(t) => match right_of_start.split_at(t) {
                    Ok((left, _right)) => left,
                    Err(_) => right_of_start,
                },
                None => right_of_start,
            };
        }
    }

    // Out-of-principal range: build the arc DIRECTLY as uniformly
    // segmented rational quadratics (n = ceil(span / (π/2)) pieces of
    // span/n each, shoulder construction per piece), with knots in
    // absolute angle over exactly [a0, a1].
    //
    // Terminal-sliver guard (PK-257): the previous construction split
    // an a0-anchored full circle down to the range, so an endpoint
    // carrying upstream atan2 noise (~1e-8 rad off an exact quarter
    // angle) stranded a sub-tolerance nano piece against the quarter-
    // knot grid — a full-multiplicity break that poisoned every
    // downstream knot union (measured PK-255: fillet-end seam arcs of
    // width π/2 ± 2e-8 emitted [main, 2e-8] composites that collapsed
    // the exact-boundary derivative rows). Uniform segmentation makes
    // sliver spans impossible by construction: every piece has width
    // span/n ∈ (π/4, π/2]. For spans that are exact multiples of π/2
    // the segment boundaries coincide with the old quarter grid.
    use crate::geometry::analytic_curves::circle_basis;
    let (u, v) = circle_basis(normal);
    uniform_segmented_conic_arc(center, u * radius, v * radius, a0, a1)
}

/// Angular window (radians) for the quarter-grid collision tests in
/// [`arc_to_nurbs`] and [`ellipse_to_nurbs`]. Well above the ~1e-8
/// atan2 noise class carried by upstream angle authors, well below
/// any legitimate arc feature; a stranded piece this size is
/// sub-tolerance geometry for kernel-scale radii.
const QUARTER_COLLISION_WINDOW: f64 = 1.0e-6;
