//! True nearest-point measurement against a pcurve's surface trace.
//!
//! The trace of a pcurve is its image under the owning face's surface
//! map. Deviation measurements (stale-pcurve detection, export gap
//! reporting) need the distance from a 3D point to that trace. A
//! nearest-of-N-samples scan is floored at half the sample spacing —
//! on a closed rim circle that floor sits above the stale tolerance,
//! so geometrically coincident pcurves measured permanently "stale"
//! (PK-287). This module refines the coarse scan to the true local
//! minimum so the measurement can certify coincidence.

use crate::geometry::Point3;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceEval;

/// Nearest point on a pcurve trace to a query point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TraceNearest {
    /// Distance from the query point to the trace at `t`.
    pub distance: f64,
    /// Pcurve parameter of the nearest trace point found.
    pub t: f64,
}

/// Golden-section iterations for the local refinement. The bracket is
/// two coarse steps wide; 48 iterations shrink it by ~0.618^48 ≈ 9e-11,
/// far below any geometric tolerance in play.
const REFINE_ITERATIONS: usize = 48;

/// Distance from `point` to the surface trace of `pcurve`, refined to
/// the true local minimum.
///
/// Scans `coarse_samples + 1` uniformly spaced parameters to bracket
/// the global minimum, then runs golden-section refinement on the
/// bracketing interval. When the coarse best sits at a parameter-range
/// endpoint, the opposite end's bracket is refined as well: on a
/// closed trace both ends map to the same 3D point, and the true
/// minimum may sit just inside either side of the seam.
///
/// Returns `TraceNearest { distance: f64::INFINITY, t: t0 }` when the
/// parameter range is not finite.
pub fn nearest_on_trace(
    pcurve: &Pcurve,
    surface: &impl SurfaceEval,
    point: &Point3,
    coarse_samples: usize,
) -> TraceNearest {
    let samples = coarse_samples.max(2);
    let (t_start, t_end) = pcurve.parameter_range();
    if !t_start.is_finite() || !t_end.is_finite() {
        return TraceNearest {
            distance: f64::INFINITY,
            t: t_start,
        };
    }
    // Descending parameter ranges are legitimate (twin-side and
    // STEP descending-domain pcurves, PK-215). Sampling and
    // bracketing run over the ascending hull; returned `t` values
    // stay in the pcurve's own domain either way.
    let (t0, t1) = if t_start <= t_end {
        (t_start, t_end)
    } else {
        (t_end, t_start)
    };
    let distance_at = |t: f64| -> f64 {
        let uv = pcurve.point_at(t);
        surface.point_at(uv.x, uv.y).distance_to(point)
    };

    let mut best_index = 0usize;
    let mut best_distance = f64::INFINITY;
    for i in 0..=samples {
        let t = t0 + (t1 - t0) * (i as f64 / samples as f64);
        let d = distance_at(t);
        if d < best_distance {
            best_distance = d;
            best_index = i;
        }
    }
    let coarse_best = TraceNearest {
        distance: best_distance,
        t: t0 + (t1 - t0) * (best_index as f64 / samples as f64),
    };

    let step = (t1 - t0) / samples as f64;
    let bracket = |index: usize| -> (f64, f64) {
        let lo = t0 + step * index.saturating_sub(1) as f64;
        let hi = t0 + step * (index + 1).min(samples) as f64;
        (lo, hi)
    };

    // The coarse best is always a candidate: refinement may only
    // improve on it, never replace it with something worse.
    let mut best = coarse_best;
    let refined = refine_bracket(bracket(best_index), &distance_at);
    if refined.distance < best.distance {
        best = refined;
    }
    // Seam handling for closed traces: both ends map to the same 3D
    // point, so a coarse best at one end may hide a better minimum
    // just inside the other end.
    if best_index == 0 {
        let other = refine_bracket(bracket(samples), &distance_at);
        if other.distance < best.distance {
            best = other;
        }
    } else if best_index == samples {
        let other = refine_bracket(bracket(0), &distance_at);
        if other.distance < best.distance {
            best = other;
        }
    }
    best
}

/// Golden-section minimization of `distance_at` over `[lo, hi]`.
fn refine_bracket((lo, hi): (f64, f64), distance_at: &impl Fn(f64) -> f64) -> TraceNearest {
    const INV_PHI: f64 = 0.618_033_988_749_894_9;
    let (mut a, mut b) = (lo, hi);
    let mut c = b - (b - a) * INV_PHI;
    let mut d = a + (b - a) * INV_PHI;
    let mut fc = distance_at(c);
    let mut fd = distance_at(d);
    for _ in 0..REFINE_ITERATIONS {
        if fc < fd {
            b = d;
            d = c;
            fd = fc;
            c = b - (b - a) * INV_PHI;
            fc = distance_at(c);
        } else {
            a = c;
            c = d;
            fc = fd;
            d = a + (b - a) * INV_PHI;
            fd = distance_at(d);
        }
    }
    // Endpoints of the original bracket are candidates too: golden
    // section converges to an interior minimum, but the true minimum
    // of a clamped bracket can sit on its boundary.
    let mut best = TraceNearest {
        distance: distance_at(lo),
        t: lo,
    };
    for (t, dist) in [(hi, distance_at(hi)), (c, fc), (d, fd)] {
        if dist < best.distance {
            best = TraceNearest { distance: dist, t };
        }
    }
    best
}

#[cfg(test)]
#[path = "pcurve_trace_tests.rs"]
mod pcurve_trace_tests;
