//! Seed-consumption test for the surface–surface marching driver:
//! distance from a candidate seed to the polylines already traced.

use crate::geometry::Point3;

/// A traced polyline registered for seed consumption.
///
/// A CLOSED march claims the whole loop, including the end→start
/// stretch the step loop never walked (closure fires when a step
/// lands within the closure tolerance of the start, leaving a gap of
/// up to that tolerance in arc). Consumption must claim the same
/// set: the closing segment is tested with the measured gap as its
/// radius — the same tolerance class the closure decision itself
/// asserted — because an on-curve point inside the gap can bulge up
/// to the gap's own sagitta off the closing chord (measured: a
/// rotated-torus lens tip seed sat 1.1e-4 off a 6.3e-4 closing
/// chord while the base radius was 1e-5, re-marched the identical
/// 5.37-arc loop, and the duplicate records stranded the boolean
/// assembly, PK-484).
pub(super) struct ConsumedTrace {
    pub points: Vec<Point3>,
    /// `Some(gap)` when the march reported closure; `gap` is the
    /// end→start distance left unwalked.
    pub closure_gap: Option<f64>,
}

impl ConsumedTrace {
    pub(super) fn open(points: Vec<Point3>) -> Self {
        Self {
            points,
            closure_gap: None,
        }
    }

    pub(super) fn closed(points: Vec<Point3>) -> Self {
        let gap = match (points.first(), points.last()) {
            (Some(a), Some(b)) => Some(a.distance_to(b)),
            _ => None,
        };
        Self {
            points,
            closure_gap: gap,
        }
    }
}

fn point_near_segment(p: &Point3, a: &Point3, b: &Point3, radius: f64) -> bool {
    let ab = *b - *a;
    let len_sq = ab.dot(&ab);
    let t = if len_sq > 0.0 {
        ((*p - *a).dot(&ab) / len_sq).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let closest = *a + ab * t;
    let d = *p - closest;
    d.dot(&d) < radius * radius
}

/// True when `p` lies within `radius` of any segment of any traced
/// polyline. The point-to-segment metric (rather than point-to-sample)
/// is what makes consumption honest for adaptively-stepped traces —
/// sample spacing may exceed the consumption radius by orders of
/// magnitude while the chord stays within tolerance of the curve.
/// Closed traces additionally test their closing segment at the
/// closure-gap radius (see [`ConsumedTrace`]).
pub(super) fn point_near_polylines(p: &Point3, traces: &[ConsumedTrace], radius: f64) -> bool {
    traces.iter().any(|trace| {
        let line = &trace.points;
        if line.len() == 1 {
            return line[0].distance_to(p) < radius;
        }
        if line
            .windows(2)
            .any(|seg| point_near_segment(p, &seg[0], &seg[1], radius))
        {
            return true;
        }
        match (trace.closure_gap, line.first(), line.last()) {
            (Some(gap), Some(first), Some(last)) => {
                point_near_segment(p, last, first, radius.max(gap))
            }
            _ => false,
        }
    })
}
