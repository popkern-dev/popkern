//! CHART-EXTERIOR PROMOTION (PK-705): the 3D-Newell election files a
//! loop as a hole by SIZE (under half the largest loop's area) when it
//! shares no edge with the elected outer — a proxy that is blind to
//! where the loop lies. On a face the chart walk partitioned, every
//! loop carries its exact chart polyline, and a hole must lie INSIDE
//! its outer there. A loop whose chart samples all lie outside the
//! outer, on every whole-period lift, and whose outer lies outside it
//! in turn, is a disjoint region of the same face, never a hole.
//!
//! Measured: a pin cylinder unioned through three disjoint knuckles
//! cut each half-cylinder face into five strips; the strip inside the
//! last knuckle measured 47.5 against a threshold of 48.7 (its 5.94
//! length vs the 12.18 exposed strips), was filed as a hole of an
//! exposed strip, and classified with it — the buried band served
//! as the exposed face's inner loop and the mesh tore (125 open
//! edges). The PK-343/PK-438 promotions cannot see this class: the
//! loops are not turnarounds and the elected outer is outside the
//! opposite operand.
//!
//! Positive proof only: a loop that wraps the period (its unrolled
//! chart polyline does not close), a sample riding the outer's
//! boundary band, or fewer than three decisive samples keeps the loop
//! filed where the election put it.

use super::fragment_chart::OrbitPolylines;
use crate::periodic_chart_parity_torus::ChartLoop;

/// Whole-period lifts of the outer tried for containment: a face whose
/// trim fits within one period meets the hole on one of these.
const LIFTS: std::ops::RangeInclusive<i32> = -2..=2;

/// Promote, from `inner_loops_out`, every loop whose chart polyline is
/// strictly exterior to the outer's (see the module docs). `view` is
/// the face's chart-walk orbit polylines (absent when the face took
/// the 3D-Newell walk: then nothing is promoted); `u_period` the
/// chart's u-period.
pub(super) fn promote_chart_exterior_holes(
    inner_loops_out: &mut Vec<Vec<usize>>,
    outer_loop: &[usize],
    view: Option<&OrbitPolylines>,
    u_period: Option<f64>,
) -> Vec<Vec<usize>> {
    let Some(view) = view else {
        return Vec::new();
    };
    let polygon = |lp: &[usize]| -> Option<Vec<(f64, f64)>> {
        unrolled_closed(view.get(lp.iter().min()?)?, u_period)
    };
    let Some(outer) = polygon(outer_loop) else {
        return Vec::new();
    };
    let mut promoted = Vec::new();
    inner_loops_out.retain(|hole| {
        let exterior = polygon(hole).is_some_and(|h| {
            strictly_exterior(&h, &outer, u_period) && strictly_exterior(&outer, &h, u_period)
        });
        if exterior {
            if popkern_core::trace::trace_enabled("fragment-walk") {
                eprintln!(
                    "    [elect] chart-exterior hole {:?}{} promoted to its own region",
                    &hole[..hole.len().min(6)],
                    if hole.len() > 6 { "…" } else { "" }
                );
            }
            promoted.push(hole.clone());
        }
        !exterior
    });
    promoted
}

/// The loop's chart points unrolled continuously in u (each step taken
/// as its nearest whole-period representative), or `None` when the
/// unrolled loop does not close — it winds round the period and has no
/// interior polygon in the chart.
fn unrolled_closed(lp: &ChartLoop, u_period: Option<f64>) -> Option<Vec<(f64, f64)>> {
    let points = lp.points();
    let first = *points.first()?;
    let mut out = vec![first];
    let mut prev = first;
    for &(u, v) in &points[1..] {
        let mut du = u - prev.0;
        if let Some(p) = u_period {
            du -= p * (du / p).round();
        }
        prev = (prev.0 + du, v);
        out.push(prev);
    }
    let raw = first.0 - prev.0;
    let turns = u_period.map_or(0.0, |p| (raw / p).round());
    (out.len() >= 3 && turns == 0.0).then_some(out)
}

fn bbox_extent(points: &[(f64, f64)]) -> f64 {
    let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
    for &(x, y) in points {
        lo = (lo.0.min(x), lo.1.min(y));
        hi = (hi.0.max(x), hi.1.max(y));
    }
    (hi.0 - lo.0) + (hi.1 - lo.1)
}

/// True when every sample of `lp` away from `ring`'s boundary band lies
/// outside `ring` on every whole-period lift, with at least three such
/// decisive samples.
fn strictly_exterior(lp: &[(f64, f64)], ring: &[(f64, f64)], u_period: Option<f64>) -> bool {
    let band = 1.0e-3 * bbox_extent(ring).max(1e-9);
    let lifts: Vec<f64> = match u_period {
        Some(p) => LIFTS.map(|k| f64::from(k) * p).collect(),
        None => vec![0.0],
    };
    let mut decisive = 0usize;
    for &(pu, pv) in lp {
        let mut near = false;
        for &shift in &lifts {
            let q = (pu - shift, pv);
            if boundary_distance(q, ring) < band {
                near = true;
                break;
            }
            if inside(q, ring) {
                return false;
            }
        }
        if !near {
            decisive += 1;
        }
    }
    decisive >= 3
}

fn boundary_distance((px, py): (f64, f64), ring: &[(f64, f64)]) -> f64 {
    let mut best = f64::MAX;
    for i in 0..ring.len() {
        let (x0, y0) = ring[i];
        let (x1, y1) = ring[(i + 1) % ring.len()];
        let (dx, dy) = (x1 - x0, y1 - y0);
        let l2 = (dx * dx + dy * dy).max(1e-300);
        let t = (((px - x0) * dx + (py - y0) * dy) / l2).clamp(0.0, 1.0);
        best = best.min(((px - x0 - dx * t).powi(2) + (py - y0 - dy * t).powi(2)).sqrt());
    }
    best
}

fn inside((px, py): (f64, f64), ring: &[(f64, f64)]) -> bool {
    let mut crossings = 0;
    for i in 0..ring.len() {
        let (x0, y0) = ring[i];
        let (x1, y1) = ring[(i + 1) % ring.len()];
        if (y0 > py) != (y1 > py) && x0 + (py - y0) / (y1 - y0) * (x1 - x0) > px {
            crossings += 1;
        }
    }
    crossings % 2 == 1
}

#[cfg(test)]
mod tests {
    use super::{strictly_exterior, unrolled_closed};
    use crate::periodic_chart_parity_torus::ChartLoop;
    use std::f64::consts::{PI, TAU};

    /// A chart rectangle sampled densely along each side, as the walk's
    /// embedded runs are.
    fn rect(u0: f64, u1: f64, v0: f64, v1: f64) -> Vec<(f64, f64)> {
        let corners = [(u0, v0), (u1, v0), (u1, v1), (u0, v1)];
        let mut out = Vec::new();
        for i in 0..4 {
            let (a, b) = (corners[i], corners[(i + 1) % 4]);
            for k in 0..8 {
                let f = f64::from(k) / 8.0;
                out.push((a.0 + (b.0 - a.0) * f, a.1 + (b.1 - a.1) * f));
            }
        }
        out
    }

    #[test]
    fn stacked_strips_are_exterior_and_nested_loops_are_not() {
        // Two strips of a half-cylinder chart, one above the other.
        let exposed = rect(0.0, PI, 35.94, 48.12);
        let buried = rect(0.0, PI, 48.12, 54.06);
        assert!(strictly_exterior(&buried, &exposed, Some(TAU)));
        assert!(strictly_exterior(&exposed, &buried, Some(TAU)));
        // A genuine hole inside the strip is not exterior.
        let hole = rect(1.0, 2.0, 40.0, 44.0);
        assert!(!strictly_exterior(&hole, &exposed, Some(TAU)));
        // A loop one period over is the same region — not exterior.
        let lifted = rect(TAU + 1.0, TAU + 2.0, 40.0, 44.0);
        assert!(!strictly_exterior(&lifted, &exposed, Some(TAU)));
    }

    #[test]
    fn a_loop_winding_the_period_has_no_chart_polygon() {
        // A rim ring: once round the period at v = 0 — its unrolled
        // end sits one period from its start.
        let rim: Vec<(f64, f64)> = (0..8).map(|k| (TAU * f64::from(k) / 8.0, 0.0)).collect();
        let band = ChartLoop::from_runs([rim.as_slice()]);
        assert!(unrolled_closed(&band, Some(TAU)).is_none());
        let strip = ChartLoop::from_runs([rect(0.0, PI, 0.0, 5.0).as_slice()]);
        assert!(unrolled_closed(&strip, Some(TAU)).is_some());
    }
}
