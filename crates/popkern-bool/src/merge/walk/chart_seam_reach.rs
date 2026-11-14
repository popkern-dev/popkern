//! The periodic coordinate's span of the linear reach at each run's
//! start (PK-1185).
//!
//! A seam's border copies are placed by a window in the periodic
//! coordinate the face's other runs leave (`chart_branches`), whose
//! slack was a fixed fraction of the period. Beside a sphere's pole the
//! chart's longitude speed falls to the distance from the pole: a ring
//! passing 2e-6 round the pole crosses the seam 2.7e-6 from it, and the
//! ring's end there reads a longitude 4e-4 rad from the seam's own (1e-9
//! apart in space). Its two halves span the period exactly, the window
//! was 1.3e-5 rad wide, the seam misfit it, and the face fell to the
//! 3D walk, which kept the polar cap for the ball. The slack is the
//! reach measured in space: the reach over the chart's speed in the
//! periodic coordinate, under a quarter period. A point within the
//! reach of the pole itself has no longitude to misread and no slack:
//! its run's longitude is the chart's convention there, which the
//! placement reads as it always has (a quarter period of slack at the
//! poles a Viviani curve passes through moved its seam a period). The window's
//! own bounds, the extremes of the other runs, carry the same slack
//! where a run ends beside the pole: a ring passing 3e-6 beside it
//! crosses the seam there, its ends read 1e-3 rad past the seam, and
//! with the ring spanning the period the window was empty.

use std::collections::BTreeMap;

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::chart_subdivision::FaceChartSubdivision;

/// For each embedded half-edge of `sub`, the span in the periodic
/// coordinate (period `period`) of `reach` at the run's start and at
/// its end, 0 where that point lies within the reach of a pole.
pub(super) fn period_reach(
    sub: &FaceChartSubdivision,
    surface: &SurfaceGeometry,
    reach: f64,
    period: f64,
) -> BTreeMap<usize, (f64, f64)> {
    let cap = 0.25 * period;
    let span_at = |&(a, b): &(f64, f64)| {
        // A transposed chart walks the surface's v as its u.
        let speed = if sub.transposed {
            surface.partial_v(b, a).length()
        } else {
            surface.partial_u(a, b).length()
        };
        if speed * cap > reach {
            reach / speed
        } else {
            0.0
        }
    };
    sub.embeddings
        .iter()
        .filter_map(|(&he, emb)| Some((he, (span_at(emb.uv.first()?), span_at(emb.uv.last()?)))))
        .collect()
}
