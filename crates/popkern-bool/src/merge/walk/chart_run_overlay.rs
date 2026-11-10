//! Whether two chart runs leaving one node are one segment twice.
//!
//! The duplicate-departure check (I5) once read two runs as duplicates
//! when they left a node in the same chart direction and ended at the
//! same node. At a sphere's pole every run leaves along the meridian,
//! so every departure has the same chart direction whatever its
//! azimuth: Viviani's curve with the ball's pole at its double point
//! leaves the pole on four branches, the two halves of each lobe
//! ending at one node, and the face's walk fell back to a 3D outline
//! that kept the lobes in a subtraction (PK-1139). A duplicate overlays
//! its partner along its whole length; runs whose middles lie apart by
//! more than a hundredth of their extent are distinct.

use popkern_core::geometry::PortableMath;

use super::chart_embedding::ChartEmbedding;

/// The fraction of the runs' extent their middles may lie apart and
/// still be one segment.
const OVERLAY_FRACTION: f64 = 1e-2;

/// The point halfway along a chart polyline, and the polyline's length.
fn halfway(uv: &[(f64, f64)]) -> Option<((f64, f64), f64)> {
    let lengths: Vec<f64> = uv
        .windows(2)
        .map(|w| (w[1].0 - w[0].0).phypot(w[1].1 - w[0].1))
        .collect();
    let total: f64 = lengths.iter().sum();
    let mut left = 0.5 * total;
    for (w, &len) in uv.windows(2).zip(&lengths) {
        if left <= len {
            let s = if len > 0.0 { left / len } else { 0.0 };
            let p = (
                w[0].0 + s * (w[1].0 - w[0].0),
                w[0].1 + s * (w[1].1 - w[0].1),
            );
            return Some((p, total));
        }
        left -= len;
    }
    uv.last().map(|&p| (p, total))
}

/// Whether runs `a` and `b` (leaving one node toward one node) overlay:
/// their halfway points coincide within [`OVERLAY_FRACTION`] of the
/// longer run. Undecidable runs (empty polylines) count as overlaying,
/// as the direction test alone decided before.
pub(crate) fn chart_runs_overlay(a: &ChartEmbedding, b: &ChartEmbedding) -> bool {
    let (Some((pa, la)), Some((pb, lb))) = (halfway(&a.uv), halfway(&b.uv)) else {
        return true;
    };
    (pa.0 - pb.0).phypot(pa.1 - pb.1) <= OVERLAY_FRACTION * la.max(lb)
}

#[cfg(test)]
#[path = "chart_run_overlay_tests.rs"]
mod tests;
