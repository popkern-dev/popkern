//! The bidirectional march wrapper and its termination trace — split
//! from `marching.rs` (500-line cap rule).

use super::marching_helpers;
use super::{MarchSample, march_one_direction, principal_branch_seed};
use crate::geometry::intersection::degeneracy::DegeneracyHandler;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

/// March along a surface-surface intersection curve from a starting point.
///
/// Marches bidirectionally (+T and -T), merges results.
/// Returns `(samples, is_closed)`.
#[allow(dead_code, clippy::too_many_arguments)]
pub(crate) fn march(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    start: &MarchSample,
    tolerance: &Tolerance,
    handler: &dyn DegeneracyHandler,
) -> (Vec<MarchSample>, bool) {
    // The SEED must sit on the principal period branch every marched
    // step is wrapped onto — mechanism and measured rationale in
    // [`principal_branch_seed`] (PK-454).
    let start = principal_branch_seed(
        start, surface1, u1_range, v1_range, surface2, u2_range, v2_range,
    );
    let start = &start;

    // March in +T direction.
    let forward = march_one_direction(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, start, 1.0, tolerance, handler,
    );

    // A forward trace that CLOSED already covers the whole loop —
    // the backward march would re-trace it lap-for-lap (measured on
    // the coil-grind SSI: every closed marched loop carried 2x
    // samples as two duplicate laps that seam-splitting turned into
    // duplicate records and consolidation had to dedupe — PK-450).
    let forward_closed = marching_helpers::forward_trace_closed(&forward, tolerance);

    // March in -T direction (only when the loop is still open).
    let mut backward = if forward_closed {
        Vec::new()
    } else {
        march_one_direction(
            surface1, u1_range, v1_range, surface2, u2_range, v2_range, start, -1.0, tolerance,
            handler,
        )
    };

    // Reverse backward so it goes from far-negative to start.
    backward.reverse();

    // Merge: backward (without duplicate start) + forward.
    // backward ends at start, forward starts at start, so drop last of backward.
    if !backward.is_empty() {
        backward.pop(); // remove duplicate start point
    }
    backward.extend(forward);
    let samples = backward;

    // Check if closed: last sample within a step-sized tolerance of first.
    // The march steps in discrete increments, so exact closure is impossible.
    // Use the average step length as the closure tolerance.
    let avg_step = if samples.len() > 1 {
        let total_arc: f64 = samples
            .windows(2)
            .map(|w| w[0].point.distance_to(&w[1].point))
            .sum();
        total_arc / (samples.len() - 1) as f64
    } else {
        tolerance.linear()
    };
    let closure_tol = (avg_step * 2.0).max(tolerance.linear());
    let is_closed = if samples.len() > 10 {
        let first = &samples[0];
        let last = &samples[samples.len() - 1];
        first.point.distance_to(&last.point) < closure_tol
    } else {
        false
    };

    (samples, is_closed)
}

/// One `march` trace line per directional trace naming WHY it ended —
/// "step-budget" means the trace was TRUNCATED mid-curve (every
/// legitimate ending names itself); downstream honest-decline gates
/// catch the missing span, but the cause is only visible here.
pub(super) fn trace_march_termination(termination: &str, samples: &[MarchSample]) {
    if crate::trace::trace_enabled("march") {
        let last = &samples[samples.len() - 1].point;
        // Arc alongside the count: a budget-saturated march with arc
        // ≈ count × min-step is a tolerance-floor crawl (no new
        // geometry per step); one with real arc is a legitimate long
        // curve. The distinction picks the containment (PK-460).
        let arc: f64 = samples
            .windows(2)
            .map(|w| w[0].point.distance_to(&w[1].point))
            .sum();
        // Closest return to the start after the closure-check
        // grace window: a budget-saturated march whose return
        // distance stays large is an open drift (spiral/multi-sheet);
        // one that repeatedly passes near the start is a closure
        // detection miss (PK-460).
        let start = &samples[0].point;
        let closest_return = samples
            .iter()
            .skip(21)
            .map(|s| s.point.distance_to(start))
            .fold(f64::INFINITY, f64::min);
        eprintln!(
            "  [march] direction ended: {termination} after {} pts arc {arc:.4} return {closest_return:.2e} at ({:.4},{:.4},{:.4})",
            samples.len(),
            last.x,
            last.y,
            last.z
        );
    }
}
