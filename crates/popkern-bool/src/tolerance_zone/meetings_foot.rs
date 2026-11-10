//! The foot of a point on a record curve within a parameter window,
//! warm-started from the previous foot.
//!
//! The meeting search refines each sampled minimum by golden section
//! and measures its zone by doubling and bisection: some 260 feet per
//! minimum, each a step from the last. A cold windowed projection
//! samples the window 65 times and polishes five seeds, about a
//! hundred evaluations a foot: a box united with a box spent 4.8 of
//! its 6 ms finding meetings of six straight records. Newton from the
//! previous foot reaches the same foot in a few evaluations; it is
//! kept only when no station of the window is nearer, and otherwise
//! the cold projection decides, so the foot is the one it would give.

use std::cell::Cell;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, closest_point_in_window};

/// Window stations the warm foot is checked against.
const CHECKS: usize = 4;
/// Newton steps from the warm seed.
const NEWTON_STEPS: usize = 12;

/// The parameter of `c` in `[w0, w1]` nearest `p`, starting from and
/// updating `seed`.
pub(super) fn warm_foot(
    c: &dyn CurveEval,
    p: &Point3,
    (w0, w1): (f64, f64),
    seed: &Cell<f64>,
) -> f64 {
    let (lo, hi) = (w0.min(w1), w0.max(w1));
    let dist = |t: f64| c.point_at(t).distance_to(p);
    let mut t = seed.get().clamp(lo, hi);
    for _ in 0..NEWTON_STEPS {
        let diff = c.point_at(t) - *p;
        let d1 = c.derivative_at(t, 1);
        let fp = d1.dot(&d1) + diff.dot(&c.derivative_at(t, 2));
        if fp.abs() < 1e-30 {
            break;
        }
        let next = (t - diff.dot(&d1) / fp).clamp(lo, hi);
        if dist(next) > dist(t) {
            break;
        }
        let done = (next - t).abs() < 1e-12 * (hi - lo).max(1e-300);
        t = next;
        if done {
            break;
        }
    }
    let best = dist(t);
    let beaten = (0..=CHECKS).any(|k| {
        let s = lo + (hi - lo) * k as f64 / CHECKS as f64;
        dist(s) < best - 1e-15
    });
    if beaten {
        t = closest_point_in_window(c, p, lo, hi);
    }
    seed.set(t);
    t
}
