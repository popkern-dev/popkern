//! Exact interior sample of a fragment on a doubly periodic chart
//! (PK-685), split from `fragment_interior` under its line cap.
//!
//! The decision is parity from a reference point known to lie in the
//! region (`crate::periodic_chart_parity_torus`). The reference comes
//! from, in order:
//!
//! 1. A SLIT of the fragment's own loops: a curve walked in both
//!    directions (a seam edge, a cut walked in and back out). It lies
//!    inside the region, so a step off its middle, to either side, is
//!    material. No walk handedness enters: the 3D-Newell walk claims
//!    a self-loop ring's direction by a fixed hand that need not match
//!    the material side, and a torus face's fragments are typically
//!    rings joined by seam slits.
//! 2. The chart walk's handedness: a step to the σ-left of the walk's
//!    own orbit polylines (material-on-left by construction).
//!
//! Anything else declines to the caller's nearest-boundary vote, with
//! the reason recorded.

use crate::merge::graph::IntersectionGraph;
use crate::merge::walk::Fragment;
use crate::periodic_chart_parity_torus::{
    ChartLoop, Periods, bounds_a_region, interior_reference, same_side, slit_reference,
    widest_material_arc,
};
use popkern_core::geometry::surface::SurfaceGeometry;

/// Why the exact torus decision declined (the caller falls back to the
/// nearest-boundary vote and records the reason).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TorusSampleDecline {
    /// No slit in the loops and no chart handedness (a 3D-Newell walk).
    NoReferenceSource,
    /// A loop could not be charted.
    Uncharted,
    /// The loops' windings do not cancel mod 2: they bound no region
    /// on the torus (an incomplete loop set), so parity is
    /// path-dependent.
    UnboundedLoops,
    /// Every candidate reference step is crowded below resolution.
    NoReference,
    /// Neither the centroid row nor the reference's row crosses
    /// material.
    NoMaterialArc,
}

impl std::fmt::Display for TorusSampleDecline {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::NoReferenceSource => "no slit and no chart handedness",
            Self::Uncharted => "a loop could not be charted",
            Self::UnboundedLoops => "loops bound no region (windings odd)",
            Self::NoReference => "no clear interior reference step",
            Self::NoMaterialArc => "no material arc on the sampled rows",
        })
    }
}

/// A fragment's region on a doubly periodic chart, decidable exactly:
/// its loops and a reference point inside it.
pub(crate) struct TorusRegion {
    loops: Vec<ChartLoop>,
    reference: (f64, f64),
    source: &'static str,
    periods: Periods,
}

impl TorusRegion {
    /// The fragment's region, from a slit reference else the chart
    /// walk's handedness.
    pub(crate) fn of(
        graph: &IntersectionGraph,
        fragment: &Fragment,
        surface: &SurfaceGeometry,
        periods: Periods,
    ) -> Result<Self, TorusSampleDecline> {
        let (loops, reference, source) = region_view(graph, fragment, surface, periods)?;
        Ok(Self {
            loops,
            reference,
            source,
            periods,
        })
    }

    /// Whether the chart point `x` lies in the region.
    pub(crate) fn contains(&self, x: (f64, f64)) -> bool {
        same_side(&self.loops, self.reference, x, self.periods)
    }
}

/// The widest material arc midpoint `(u, v)` of the fragment's region
/// on the row `cv`, else on the reference's own row (which meets
/// material by construction).
pub(crate) fn exact_torus_sample(
    graph: &IntersectionGraph,
    fragment: &Fragment,
    surface: &SurfaceGeometry,
    cv: f64,
    periods: Periods,
) -> Result<(f64, f64), TorusSampleDecline> {
    let trace = popkern_core::trace::trace_enabled("face-selection");
    let TorusRegion {
        loops,
        reference,
        source,
        ..
    } = TorusRegion::of(graph, fragment, surface, periods)?;
    for v in [cv, reference.1] {
        if let Some(u) = widest_material_arc(&loops, reference, v, periods) {
            if trace {
                eprintln!(
                    "      [selpoly-torus] {source} reference=({:.4},{:.4}) row v={v:.4} -> u={u:.4}",
                    reference.0, reference.1
                );
            }
            return Ok((u, v));
        }
    }
    Err(TorusSampleDecline::NoMaterialArc)
}

/// The region's loops, reference and reference source.
#[allow(clippy::type_complexity)]
fn region_view(
    graph: &IntersectionGraph,
    fragment: &Fragment,
    surface: &SurfaceGeometry,
    periods: Periods,
) -> Result<(Vec<ChartLoop>, (f64, f64), &'static str), TorusSampleDecline> {
    Ok(match slit_view(graph, fragment, surface, periods)? {
        Some((loops, reference)) => (loops, reference, "slit"),
        None => {
            let chart = fragment
                .chart
                .as_ref()
                .ok_or(TorusSampleDecline::NoReferenceSource)?;
            let loops = chart.loops.clone();
            if !bounds_a_region(&loops, periods) {
                return Err(TorusSampleDecline::UnboundedLoops);
            }
            let reference = interior_reference(&loops, chart.sigma, periods)
                .ok_or(TorusSampleDecline::NoReference)?;
            (loops, reference, "chart handedness")
        }
    })
}

/// The fragment's charted loops and a slit reference, when its loops
/// walk some edge in both directions; `Ok(None)` when they walk none.
#[allow(clippy::type_complexity)]
fn slit_view(
    graph: &IntersectionGraph,
    fragment: &Fragment,
    surface: &SurfaceGeometry,
    periods: Periods,
) -> Result<Option<(Vec<ChartLoop>, (f64, f64))>, TorusSampleDecline> {
    let walked: std::collections::HashSet<usize> = std::iter::once(&fragment.outer_loop)
        .chain(fragment.inner_loops.iter())
        .flatten()
        .copied()
        .collect();
    let is_slit = |h: usize| {
        let twin = graph.half_edges[h].twin;
        twin != h && walked.contains(&twin)
    };
    if !walked.iter().any(|&h| is_slit(h)) {
        return Ok(None);
    }
    let mut runs = Vec::new();
    for lp in std::iter::once(&fragment.outer_loop).chain(fragment.inner_loops.iter()) {
        runs.push(
            super::fragment_sampling::loop_uv_runs(graph, &fragment.face, lp, surface)
                .ok_or(TorusSampleDecline::Uncharted)?,
        );
    }
    let loops: Vec<ChartLoop> = runs
        .iter()
        .map(|r| ChartLoop::from_runs(r.iter().map(|(_, run)| run.as_slice())))
        .collect();
    if !bounds_a_region(&loops, periods) {
        return Err(TorusSampleDecline::UnboundedLoops);
    }
    // Slits in loop order (deterministic); the first with a clear step.
    let reference = runs
        .iter()
        .flatten()
        .filter(|(h, run)| is_slit(*h) && run.len() >= 2)
        .find_map(|(_, run)| slit_reference(&loops, run, periods))
        .ok_or(TorusSampleDecline::NoReference)?;
    Ok(Some((loops, reference)))
}
