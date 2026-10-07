//! The chart walk's own view of a fragment (PK-685), split from
//! `walk/mod.rs` under its line ratchet.

use std::collections::HashMap;

use super::chart_subdivision::FaceChartSubdivision;
use crate::periodic_chart_parity_torus::ChartLoop;

/// Orbit polylines of one face's chart walk, keyed by each orbit's
/// smallest half-edge (loop identity survives rotation).
pub(super) type OrbitPolylines = HashMap<usize, ChartLoop>;

/// A chart-walked fragment's loops as the walk embedded them: every
/// loop (outer first, then the inner loops in order) is a
/// material-on-LEFT orbit of the σ-oriented chart, σ self-calibrated
/// from the face's stored trim ring (`chart_subdivision`). Consumers
/// that need the region's side read it here rather than re-deriving
/// run directions: a closed ring's direction is the embedding's, and
/// a re-sampled run can disagree with it.
#[derive(Debug, Clone)]
pub(crate) struct FragmentChart {
    /// +1 when chart CCW is CCW about the outward normal, −1 otherwise.
    pub sigma: f64,
    /// Orbit loops (one run per half-edge embedding), outer loop first.
    pub loops: Vec<ChartLoop>,
}

/// The chart view of the fragment with loops `outer` and `inners`,
/// from the walk's orbit polylines keyed by each orbit's smallest
/// half-edge. `None` unless the chart walk produced every loop.
fn fragment_chart_of(
    (sigma, polylines): &(f64, OrbitPolylines),
    outer: &[usize],
    inners: &[Vec<usize>],
) -> Option<FragmentChart> {
    let polyline = |lp: &[usize]| polylines.get(lp.iter().min()?).cloned();
    let loops = std::iter::once(outer)
        .chain(inners.iter().map(|l| l.as_slice()))
        .map(polyline)
        .collect::<Option<Vec<_>>>()?;
    Some(FragmentChart {
        sigma: *sigma,
        loops,
    })
}

/// The chart walk's handedness and orbit polylines for one face.
pub(super) fn walk_view(
    sub: &FaceChartSubdivision,
    orbits: &[Vec<usize>],
) -> (f64, OrbitPolylines) {
    let polylines = orbits
        .iter()
        .filter_map(|o| {
            let runs = o
                .iter()
                .filter_map(|h| sub.embeddings.get(h))
                .map(|e| e.uv.as_slice());
            Some((*o.iter().min()?, ChartLoop::from_runs(runs)))
        })
        .collect();
    (sub.sigma, polylines)
}

/// Attach each fragment's chart view, for fragments of faces the chart
/// walk produced (`views`, keyed by face).
pub(super) fn attach(
    fragments: &mut [super::Fragment],
    views: &HashMap<crate::merge::graph::FaceRef, (f64, OrbitPolylines)>,
) {
    for fragment in fragments {
        if let Some(view) = views.get(&fragment.face) {
            fragment.chart = fragment_chart_of(view, &fragment.outer_loop, &fragment.inner_loops);
        }
    }
}
