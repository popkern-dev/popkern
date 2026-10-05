//! Exact rational-quadratic spellings of circular arcs in a 2-D frame.
//!
//! Every arc here is a run of rational quadratic Bézier spans, each
//! spanning at most 90° so the middle weight `cos(Δ/2)` stays well
//! conditioned. The spelling is exact — no sampling, no fitting — and
//! it is the form the sweep kernel's canonical recognition maps to a
//! true `Cylinder`. Knot values are ANGLE labels: each span's knot
//! interval is its angular sub-range, so a wire's parameter agrees
//! with the arc's angle exactly at every span joint (and within a few
//! percent mid-span, where a rational quadratic's parameter is not
//! angle-linear).

use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_math::PortableMath;

/// One quadratic span `(P0, Pm, P2, w_mid)` with its `[t0, t1]` knot
/// labels; end weights are `1`.
#[derive(Debug, Clone, Copy)]
pub(crate) struct QuadraticSpan {
    /// Start control point.
    pub p0: Point3,
    /// Middle (shoulder) control point.
    pub pm: Point3,
    /// End control point.
    pub p2: Point3,
    /// Middle weight (`1` for a straight span).
    pub w: f64,
    /// Knot label at the span start.
    pub t0: f64,
    /// Knot label at the span end.
    pub t1: f64,
}

/// A straight span from `a` to `b` with knot labels `[t0, t1]`.
pub(crate) fn line_span(a: Point3, b: Point3, t0: f64, t1: f64) -> QuadraticSpan {
    QuadraticSpan {
        p0: a,
        pm: Point3::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0, (a.z + b.z) / 2.0),
        p2: b,
        w: 1.0,
        t0,
        t1,
    }
}

/// Exact spans for the arc of centre `(cx, cy)`, radius `r`, running
/// counter-clockwise from angle `a0` to `a1` (`a1 > a0`, at most one
/// full turn), placed at height `z`. Knot labels run from `a0` to
/// `a1`. Arcs wider than 90° are subdivided into equal sub-spans.
pub(crate) fn arc_spans(cx: f64, cy: f64, r: f64, a0: f64, a1: f64, z: f64) -> Vec<QuadraticSpan> {
    let sweep = a1 - a0;
    let n_sub = (sweep / std::f64::consts::FRAC_PI_2 - 1.0e-12)
        .ceil()
        .max(1.0) as usize;
    let step = sweep / n_sub as f64;
    let w = (step / 2.0).pcos();
    let on = |ang: f64, scale: f64| {
        Point3::new(cx + scale * r * ang.pcos(), cy + scale * r * ang.psin(), z)
    };
    (0..n_sub)
        .map(|i| {
            let sa = a0 + step * i as f64;
            let sb = if i + 1 == n_sub { a1 } else { sa + step };
            QuadraticSpan {
                p0: on(sa, 1.0),
                pm: on((sa + sb) / 2.0, 1.0 / w),
                p2: on(sb, 1.0),
                w,
                t0: sa,
                t1: sb,
            }
        })
        .collect()
}

/// Assemble consecutive spans (each starting where the previous ended)
/// into one clamped degree-2 rational NURBS with FULL-multiplicity
/// interior knots at every joint, so every span is an exact
/// independent piece.
pub(crate) fn assemble_quadratic(spans: &[QuadraticSpan]) -> NurbsCurve {
    let n = spans.len();
    let mut control_points = Vec::with_capacity(2 * n + 1);
    let mut weights = Vec::with_capacity(2 * n + 1);
    let mut knots = Vec::with_capacity(2 * n + 4);
    let t_first = spans.first().map(|s| s.t0).unwrap_or(0.0);
    knots.extend_from_slice(&[t_first; 3]);
    for (i, s) in spans.iter().enumerate() {
        if i == 0 {
            control_points.push(s.p0);
            weights.push(1.0);
        } else {
            knots.push(s.t0);
            knots.push(s.t0);
        }
        control_points.push(s.pm);
        weights.push(s.w);
        control_points.push(s.p2);
        weights.push(1.0);
    }
    let t_last = spans.last().map(|s| s.t1).unwrap_or(1.0);
    knots.extend_from_slice(&[t_last; 3]);
    NurbsCurve {
        degree: 2,
        knots,
        control_points,
        weights,
    }
}

/// Exact full circle in the standard 9-control-point form (four
/// quadrant spans, knots `0, ¼, ½, ¾, 1`) — the spelling the sweep
/// kernel's canonical recognition maps to a true `Cylinder`.
pub(crate) fn circle_wire(cx: f64, cy: f64, r: f64, z: f64) -> NurbsCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let pts = vec![
        Point3::new(cx + r, cy, z),
        Point3::new(cx + r, cy + r, z),
        Point3::new(cx, cy + r, z),
        Point3::new(cx - r, cy + r, z),
        Point3::new(cx - r, cy, z),
        Point3::new(cx - r, cy - r, z),
        Point3::new(cx, cy - r, z),
        Point3::new(cx + r, cy - r, z),
        Point3::new(cx + r, cy, z),
    ];
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: pts,
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Closed degree-1 wire through `points` (first point NOT repeated)
/// at height `z`, knots uniform on `[0, 1]`.
pub(crate) fn polygon_wire(points: &[(f64, f64)], z: f64) -> NurbsCurve {
    let mut pts: Vec<Point3> = points.iter().map(|&(x, y)| Point3::new(x, y, z)).collect();
    pts.push(pts[0]);
    let n = pts.len();
    let mut knots = vec![0.0];
    for i in 0..n {
        knots.push(i as f64 / (n - 1) as f64);
    }
    knots.push(1.0);
    NurbsCurve {
        degree: 1,
        knots,
        control_points: pts,
        weights: vec![1.0; n],
    }
}

#[cfg(test)]
#[path = "exact_arc_tests.rs"]
mod tests;
