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

pub(super) fn unwrap_uv_for_query(
    uv: UvPoint,
    query: UvPoint,
    periodic_u: bool,
    periodic_v: bool,
) -> UvPoint {
    let mut u = uv.u;
    if periodic_u {
        while u - query.u > TAU * 0.5 {
            u -= TAU;
        }
        while query.u - u > TAU * 0.5 {
            u += TAU;
        }
    }
    let mut v = uv.v;
    if periodic_v {
        while v - query.v > TAU * 0.5 {
            v -= TAU;
        }
        while query.v - v > TAU * 0.5 {
            v += TAU;
        }
    }
    UvPoint { u, v }
}

pub(super) fn push_distinct_uv(polygon: &mut Vec<UvPoint>, uv: UvPoint) {
    let duplicate = polygon
        .last()
        .is_some_and(|last| (last.u - uv.u).hypot(last.v - uv.v) < 1.0e-12);
    if !duplicate {
        polygon.push(uv);
    }
}

/// Periodic-aware containment: a point on a periodic surface is
/// covered by the trim if ANY period image of it lies inside the
/// chain-unwrapped loop polygon. The polygon is branch-coherent per
/// axis (see `loop_uv_polygon`), so testing the query shifted by the
/// whole periods that land inside the polygon's extent is exact.
/// Both axes are handled — a torus is periodic in u AND v, and
/// u-only image search left v-seam-crossing loops classified against
/// a torn polygon (measured: rotated torus×torus union, phantom
/// SelfIntersection between an A-face and a B-face 0.267 from the
/// B-face's true boundary; 18.75% of the B-chart read as
/// double-covered under a 48×48 sweep).
pub(super) fn point_in_uv_loop_periodic(
    point: UvPoint,
    loop_polygon: &[UvPoint],
    tolerance: f64,
    periodic_u: bool,
    periodic_v: bool,
) -> TrimContainment {
    let image_range = |periodic: bool, coord: f64, extent: Option<(f64, f64)>| -> (i64, i64) {
        if !periodic {
            return (0, 0);
        }
        match extent {
            Some((lo, hi)) => (
                ((lo - tolerance - coord) / TAU).floor() as i64,
                ((hi + tolerance - coord) / TAU).ceil() as i64,
            ),
            None => (0, 0),
        }
    };
    let extent = |pick: fn(&UvPoint) -> f64| -> Option<(f64, f64)> {
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for p in loop_polygon {
            lo = lo.min(pick(p));
            hi = hi.max(pick(p));
        }
        (lo.is_finite() && hi.is_finite()).then_some((lo, hi))
    };
    let (ku_lo, ku_hi) = image_range(periodic_u, point.u, extent(|p| p.u));
    let (kv_lo, kv_hi) = image_range(periodic_v, point.v, extent(|p| p.v));
    let mut strongest = TrimContainment::Outside;
    for ku in ku_lo..=ku_hi {
        for kv in kv_lo..=kv_hi {
            let shifted = UvPoint {
                u: point.u + ku as f64 * TAU,
                v: point.v + kv as f64 * TAU,
            };
            match point_in_uv_loop(shifted, loop_polygon, tolerance) {
                TrimContainment::Interior => return TrimContainment::Interior,
                TrimContainment::Boundary => strongest = TrimContainment::Boundary,
                TrimContainment::Outside => {}
            }
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

/// A raw (un-unwrapped) loop boundary sample, flagged when the loop
/// traversal reaches it by crossing a parametrically degenerate row
/// with a period-sized branch jump (sphere pole, cone apex).
pub(super) struct RawLoopSample {
    pub(super) uv: UvPoint,
    pub(super) pole_jump_before: bool,
}

/// Assemble a loop's UV trim polygon from raw samples.
///
/// Default assembly is chain-continuous unwrapping: each sample folds
/// to the branch of the previous point (the first anchors against the
/// query). That is correct for every loop that stays on one branch —
/// unwrapping against the query instead folds half-period loops into
/// self-overlapping polygons (measured: laptop-stand phantom
/// self-intersection), and tearing at pole crossings breaks
/// pole-touching octant caps (measured: rounded-cube corner sphere
/// classified far-side azimuths Interior).
///
/// The exception is a loop whose chains connect ONLY through
/// degenerate rows and net a full period — the 1-face sphere's
/// meridian seam loop. Pure folding collapses that polygon to a
/// zero-area strip (every query classifies Outside, which made the
/// material-sense resolver unsampleable). For those loops the chains
/// after a pole crossing keep their RAW branch, with the implied
/// segments along the degenerate rows serving as the (3D-zero-length)
/// bridges. The raw variant is used ONLY when the folded polygon is
/// area-degenerate and the closure jump is licensed by degenerate
/// rows, so every previously-served polygon assembles bit-identically.
pub(super) fn assemble_periodic_polygon(
    raw: &[RawLoopSample],
    query: UvPoint,
    periodic_u: bool,
    periodic_v: bool,
    closure_on_degenerate_rows: bool,
) -> Vec<UvPoint> {
    let chain_starts: Vec<usize> = raw
        .iter()
        .enumerate()
        .filter_map(|(i, s)| s.pole_jump_before.then_some(i))
        .collect();

    let assemble = |raw_branch_chains: &[usize]| -> Vec<UvPoint> {
        let mut polygon: Vec<UvPoint> = Vec::with_capacity(raw.len());
        let mut anchor_k = 0.0_f64;
        for (i, s) in raw.iter().enumerate() {
            let reference = polygon.last().copied().unwrap_or(query);
            let folded = unwrap_uv_for_query(s.uv, reference, periodic_u, periodic_v);
            if i == 0 {
                anchor_k = ((folded.u - s.uv.u) / TAU).round();
                push_distinct_uv(&mut polygon, folded);
                continue;
            }
            if s.pole_jump_before && raw_branch_chains.contains(&i) {
                // Bridge along the degenerate row: enter the row at the
                // previous chain's branch, exit at this chain's raw
                // branch (shifted by the loop's query anchoring).
                push_distinct_uv(&mut polygon, folded);
                let v = unwrap_uv_for_query(s.uv, reference, false, periodic_v).v;
                push_distinct_uv(
                    &mut polygon,
                    UvPoint {
                        u: s.uv.u + anchor_k * TAU,
                        v,
                    },
                );
                continue;
            }
            push_distinct_uv(&mut polygon, folded);
        }
        polygon
    };

    let folded_polygon = assemble(&[]);
    if chain_starts.is_empty() || chain_starts.len() > 4 || !periodic_u {
        return folded_polygon;
    }

    let closure_ok = |polygon: &[UvPoint]| -> bool {
        match (polygon.first(), polygon.last()) {
            (Some(a), Some(b)) => {
                (a.u - b.u).abs() <= TAU * 0.5 + 1.0e-9 || closure_on_degenerate_rows
            }
            _ => false,
        }
    };
    let area = |polygon: &[UvPoint]| -> f64 {
        let mut sum = 0.0;
        for (i, a) in polygon.iter().enumerate() {
            let b = &polygon[(i + 1) % polygon.len()];
            sum += a.u * b.v - b.u * a.v;
        }
        sum * 0.5
    };

    if closure_ok(&folded_polygon) && area(&folded_polygon).abs() > 1.0e-9 {
        return folded_polygon;
    }

    // Folded polygon is area-degenerate (or tears its closure):
    // enumerate raw-branch assignments for the pole-separated chains,
    // all-fold first, and serve the first valid non-degenerate variant.
    let m = chain_starts.len();
    for mask in 1u32..(1 << m) {
        let selected: Vec<usize> = chain_starts
            .iter()
            .enumerate()
            .filter(|(bit, _)| (mask >> bit) & 1 == 1)
            .map(|(_, &idx)| idx)
            .collect();
        let candidate = assemble(&selected);
        if candidate.len() >= 3 && closure_ok(&candidate) && area(&candidate).abs() > 1.0e-9 {
            return candidate;
        }
    }
    folded_polygon
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A v-seam-crossing torus trim: the loop's coherent (chain-
    /// unwrapped) polygon spans v ∈ [5.8, 6.9] — beyond the canonical
    /// cell. A query at the canonical image v = 0.35 lies inside only
    /// via its +τ period image; u-only image search misses it (the
    /// measured pk305 phantom-SelfIntersection class, inverted: here
    /// the miss direction is a false OUTSIDE, the boolean fixture hit
    /// the paired false INTERIOR from the torn polygon).
    #[test]
    fn v_periodic_image_found_for_seam_crossing_loop() {
        let polygon = [
            UvPoint { u: 1.0, v: 5.8 },
            UvPoint { u: 2.0, v: 5.8 },
            UvPoint { u: 2.0, v: 6.9 },
            UvPoint { u: 1.0, v: 6.9 },
        ];
        let query = UvPoint { u: 1.5, v: 0.35 };
        assert_eq!(
            point_in_uv_loop_periodic(query, &polygon, 1.0e-9, true, true),
            TrimContainment::Interior,
            "query's +tau v-image lies inside the seam-crossing loop"
        );
        assert_eq!(
            point_in_uv_loop_periodic(query, &polygon, 1.0e-9, true, false),
            TrimContainment::Outside,
            "without v-periodicity the canonical query is outside"
        );
    }

    /// Bi-axis chain unwrapping: a sample one period away in v snaps
    /// to the branch of its reference point.
    #[test]
    fn unwrap_resolves_v_branch_against_reference() {
        let reference = UvPoint { u: 1.0, v: 6.0 };
        let sample = UvPoint { u: 1.0, v: 0.2 };
        let unwrapped = unwrap_uv_for_query(sample, reference, true, true);
        assert!((unwrapped.v - (0.2 + TAU)).abs() < 1.0e-12);
        let clamped = unwrap_uv_for_query(sample, reference, true, false);
        assert!((clamped.v - 0.2).abs() < 1.0e-12);
    }
}
