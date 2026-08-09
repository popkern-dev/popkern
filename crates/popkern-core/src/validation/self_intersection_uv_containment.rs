//! UV trim-polygon containment for the self-intersection validator:
//! period-branch unwrapping, polygon assembly helpers, and the
//! point-in-loop classification (periodic-aware).
//!
//! Split out of `self_intersection_trim.rs` to keep it under the
//! 500-line structural cap.

use std::f64::consts::TAU;

#[derive(Debug, Clone, Copy)]
pub(in crate::validation) struct UvPoint {
    pub(in crate::validation) u: f64,
    pub(in crate::validation) v: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(in crate::validation) enum TrimContainment {
    Outside,
    Boundary,
    Interior,
}

pub(super) fn unwrap_uv_for_query(uv: UvPoint, query: UvPoint, periodic_u: bool) -> UvPoint {
    if !periodic_u {
        return uv;
    }

    let mut u = uv.u;
    while u - query.u > TAU * 0.5 {
        u -= TAU;
    }
    while query.u - u > TAU * 0.5 {
        u += TAU;
    }
    UvPoint { u, v: uv.v }
}

pub(super) fn push_distinct_uv(polygon: &mut Vec<UvPoint>, uv: UvPoint) {
    let duplicate = polygon
        .last()
        .is_some_and(|last| (last.u - uv.u).hypot(last.v - uv.v) < 1.0e-12);
    if !duplicate {
        polygon.push(uv);
    }
}

/// Periodic-aware containment: a point on a u-periodic surface is
/// covered by the trim if ANY period image of it lies inside the
/// chain-unwrapped loop polygon. The polygon is branch-coherent (see
/// `loop_uv_polygon`), so testing the query shifted by the whole
/// periods that land inside the polygon's u-extent is exact.
pub(super) fn point_in_uv_loop_periodic(
    point: UvPoint,
    loop_polygon: &[UvPoint],
    tolerance: f64,
    periodic_u: bool,
) -> TrimContainment {
    if !periodic_u {
        return point_in_uv_loop(point, loop_polygon, tolerance);
    }
    let (mut min_u, mut max_u) = (f64::INFINITY, f64::NEG_INFINITY);
    for p in loop_polygon {
        min_u = min_u.min(p.u);
        max_u = max_u.max(p.u);
    }
    if !min_u.is_finite() || !max_u.is_finite() {
        return point_in_uv_loop(point, loop_polygon, tolerance);
    }
    let k_lo = ((min_u - tolerance - point.u) / TAU).floor() as i64;
    let k_hi = ((max_u + tolerance - point.u) / TAU).ceil() as i64;
    let mut strongest = TrimContainment::Outside;
    for k in k_lo..=k_hi {
        let shifted = UvPoint {
            u: point.u + k as f64 * TAU,
            v: point.v,
        };
        match point_in_uv_loop(shifted, loop_polygon, tolerance) {
            TrimContainment::Interior => return TrimContainment::Interior,
            TrimContainment::Boundary => strongest = TrimContainment::Boundary,
            TrimContainment::Outside => {}
        }
    }
    strongest
}

fn point_in_uv_loop(point: UvPoint, loop_polygon: &[UvPoint], tolerance: f64) -> TrimContainment {
    if loop_polygon.len() < 3 {
        return TrimContainment::Outside;
    }

    let Some(&last) = loop_polygon.last() else {
        return TrimContainment::Outside;
    };
    let mut previous = last;
    for &current in loop_polygon {
        if distance_to_uv_segment(point, previous, current) <= tolerance {
            return TrimContainment::Boundary;
        }
        previous = current;
    }

    let mut inside = false;
    previous = last;
    for &current in loop_polygon {
        let crosses_v = (current.v > point.v) != (previous.v > point.v);
        if crosses_v {
            let u_at_v = (previous.u - current.u) * (point.v - current.v)
                / (previous.v - current.v)
                + current.u;
            if point.u < u_at_v {
                inside = !inside;
            }
        }
        previous = current;
    }
    if inside {
        TrimContainment::Interior
    } else {
        TrimContainment::Outside
    }
}

fn distance_to_uv_segment(point: UvPoint, a: UvPoint, b: UvPoint) -> f64 {
    let ab_u = b.u - a.u;
    let ab_v = b.v - a.v;
    let ap_u = point.u - a.u;
    let ap_v = point.v - a.v;
    let len_sq = ab_u * ab_u + ab_v * ab_v;
    if len_sq <= f64::EPSILON {
        return ap_u.hypot(ap_v);
    }

    let t = ((ap_u * ab_u + ap_v * ab_v) / len_sq).clamp(0.0, 1.0);
    let closest = UvPoint {
        u: a.u + ab_u * t,
        v: a.v + ab_v * t,
    };
    (point.u - closest.u).hypot(point.v - closest.v)
}
