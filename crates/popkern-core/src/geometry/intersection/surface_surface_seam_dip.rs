//! Marched traces that cross a chart seam without winding around the
//! periodic axis (PK-1030).
//!
//! The march normalizes every sample's angles into the window's
//! period, so a trace that dips across a chart seam and comes back
//! carries two opposite period jumps with no net winding. Fitted as
//! one record (the mixed-pair wrapped-record path), the uv pcurve
//! swept the whole chart at each jump. Measured: a rim fillet's torus
//! strip, whose v seam IS its rim, cut by walls y = ±0.9: the face
//! clip cut the oval at fake rim crossings into ~1e-3 pieces, 24 open
//! record ends per boolean, and the marching arm's walk left a loop
//! open by 0.87. Such a trace is split at its refined seam crossings
//! like any periodic pair's (`surface_surface_seam_split`): the seam
//! is where the face's chart ends, a rim or a seam edge.
//!
//! Unwinding the dip onto one branch instead serves the torus strip
//! but not a full-revolution face: on a hex head's chamfer cone the
//! unwound pcurve leaves the face's chart past its seam edge and the
//! merge strands the record (22 fastener and thread tests).
//!
//! A trace that winds around the axis (net displacement of half a
//! period or more) keeps the wrapped-record path it is tuned on.

use crate::geometry::intersection::marching::MarchSample;

/// Whether `samples` cross a periodic axis's seam and come back: they
/// jump by more than half a period, their minimal-image unwrapped run
/// returns within half a period of its start, and that run leaves the
/// range of the raw samples. The raw samples lie in the window's one
/// period, so only a run that crosses the seam leaves it; a sample
/// registered on another period image and straight back (a grazing
/// strip × cap contact: three samples, a run 3e-3 wide) jumps without
/// crossing and stays one record. `periods` are the (u1, v1, u2, v2)
/// periods; `None` or a non-finite period never dips.
pub(crate) fn dips_across_seam(samples: &[MarchSample], periods: [Option<f64>; 4]) -> bool {
    let coord = |s: &MarchSample, axis: usize| match axis {
        0 => s.u1,
        1 => s.v1,
        2 => s.u2,
        _ => s.v2,
    };
    periods.into_iter().enumerate().any(|(axis, period)| {
        let Some(p) = period.filter(|p| p.is_finite() && *p > 0.0) else {
            return false;
        };
        let Some(first) = samples.first().map(|s| coord(s, axis)) else {
            return false;
        };
        let (mut jumped, mut at) = (false, first);
        let (mut raw_lo, mut raw_hi) = (first, first);
        let (mut run_lo, mut run_hi) = (first, first);
        for w in samples.windows(2) {
            let (a, b) = (coord(&w[0], axis), coord(&w[1], axis));
            let d = b - a;
            jumped |= d.abs() > 0.5 * p;
            at += d - p * (d / p).round();
            (raw_lo, raw_hi) = (raw_lo.min(b), raw_hi.max(b));
            (run_lo, run_hi) = (run_lo.min(at), run_hi.max(at));
        }
        jumped && (at - first).abs() < 0.5 * p && (run_lo < raw_lo || run_hi > raw_hi)
    })
}

#[cfg(test)]
#[path = "surface_surface_seam_dip_tests.rs"]
mod tests;
