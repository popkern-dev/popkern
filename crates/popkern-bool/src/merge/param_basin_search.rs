//! Winding-aware closest-parameter search (PK-670).
//!
//! A curve that WINDS — a thread's root or crest edge, a coil seam —
//! passes near any query point once per turn, one distance basin per
//! turn. A fixed-count scan that keeps only its globally nearest
//! SAMPLE aliases once the stations fall below about two per turn: an
//! adjacent turn's sample can sit nearer than every sample of the
//! true turn, and the bracket refined around it converges one pitch
//! off (PK-669 measured the class in `split_pcurve_at_image`: 64
//! samples 78° apart on a 13.9-turn root edge). The searches here
//! scale the scan with the curve's winding — two stations per NURBS
//! knot span (a rational helix carries one span per quarter turn) or
//! per chart span its image crosses — and refine EVERY local minimum
//! of the scan, keeping the nearest refined point. The global-best
//! candidate is refined first, exactly as the fixed-count searches
//! did, so a search whose nearest sample already sat in the true
//! basin returns the same parameter.

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Scan station cap: a 1000-turn curve at four spans per turn.
const MAX_STATIONS: usize = 8192;
/// Refinement candidates per search (nearest sampled first).
const MAX_BASINS: usize = 64;

/// Distinct knot values strictly inside `(lo, hi)`, plus one: the
/// number of knot spans the window covers.
fn spans_in(knots: &[f64], lo: f64, hi: f64) -> usize {
    let mut count = 0usize;
    let mut last = f64::NEG_INFINITY;
    for &k in knots {
        if k > lo && k < hi && k > last {
            count += 1;
            last = k;
        }
    }
    count + 1
}

/// Scan stations for a search over the window `[lo, hi]` of a 3D
/// curve: `base`, raised to two per NURBS knot span in the window.
pub(crate) fn curve_scan_count(curve: &CurveGeometry, lo: f64, hi: f64, base: usize) -> usize {
    let (lo, hi) = (lo.min(hi), lo.max(hi));
    match curve {
        CurveGeometry::Nurbs(n) => (2 * spans_in(&n.knots, lo, hi)).clamp(base, MAX_STATIONS),
        _ => base,
    }
}

/// Scan stations for a search along a pcurve's IMAGE on `surface`:
/// `base`, raised to two per chart span the pcurve crosses — knot
/// spans of a NURBS host, eight per period of a periodic analytic
/// host — and two per knot span of a NURBS pcurve itself.
pub(crate) fn pcurve_image_scan_count(
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    base: usize,
) -> usize {
    let (t0, t1) = (pcurve.t_start, pcurve.t_end);
    let (mut u_lo, mut u_hi, mut v_lo, mut v_hi) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for i in 0..=base {
        let uv = pcurve.point_at(t0 + (t1 - t0) * i as f64 / base as f64);
        u_lo = u_lo.min(uv.x);
        u_hi = u_hi.max(uv.x);
        v_lo = v_lo.min(uv.y);
        v_hi = v_hi.max(uv.y);
    }
    let chart = match surface {
        SurfaceGeometry::Nurbs(s) => {
            2 * (spans_in(&s.knots_u, u_lo, u_hi) + spans_in(&s.knots_v, v_lo, v_hi))
        }
        _ => {
            let (pu, pv) = surface.parameter_periodicity();
            let wraps = |period: Option<f64>, extent: f64| match period {
                Some(p) if p > 0.0 && extent.is_finite() => (extent / p).ceil() as usize,
                _ => 0,
            };
            8 * (wraps(pu, u_hi - u_lo) + wraps(pv, v_hi - v_lo))
        }
    };
    let own = match &pcurve.geometry {
        PcurveGeometry::Nurbs(n) => 2 * spans_in(&n.knots, t0.min(t1), t0.max(t1)),
        _ => 0,
    };
    chart.max(own).clamp(base, MAX_STATIONS)
}

/// The parameter in `[lo, hi]` minimizing `dist`, with that distance.
///
/// The historical fixed-count search runs first: `legacy + 1` evenly
/// spaced stations, the FIRST strict minimum refined by the caller's
/// own `refine(center, half_width)`. When that lands within `exact`
/// (a point ON the curve) it is the answer, bit-identical to the
/// fixed-count search. Otherwise every local minimum of a
/// `dense + 1`-station scan is refined too (nearest sampled first, at
/// most [`MAX_BASINS`]) and the nearest refined point wins — a basin
/// replaces the historical answer only when strictly nearer, so the
/// adjacent-turn alias is corrected and nothing else moves.
pub(crate) fn nearest_param_all_basins(
    dist: &dyn Fn(f64) -> f64,
    (lo, hi): (f64, f64),
    (legacy, dense): (usize, usize),
    exact: f64,
    refine: &dyn Fn(f64, f64) -> f64,
) -> (f64, f64) {
    let scan = |count: usize| -> Vec<f64> {
        (0..=count)
            .map(|i| dist(lo + (hi - lo) * i as f64 / count as f64))
            .collect()
    };
    let legacy = legacy.max(1);
    let coarse = scan(legacy);
    let first_best =
        (0..coarse.len()).fold(0usize, |b, i| if coarse[i] < coarse[b] { i } else { b });
    let t = refine(
        lo + (hi - lo) * first_best as f64 / legacy as f64,
        (hi - lo) / legacy as f64,
    );
    let mut best = (t, dist(t));
    if best.1 <= exact {
        return best;
    }
    let dense = dense.max(legacy);
    let fine = if dense == legacy { coarse } else { scan(dense) };
    let mut minima: Vec<usize> = (0..fine.len())
        .filter(|&i| {
            let left = i == 0 || fine[i] <= fine[i - 1];
            let right = i + 1 == fine.len() || fine[i] <= fine[i + 1];
            left && right
        })
        .collect();
    minima.sort_by(|&a, &b| fine[a].total_cmp(&fine[b]).then(a.cmp(&b)));
    for index in minima.into_iter().take(MAX_BASINS) {
        let t = refine(
            lo + (hi - lo) * index as f64 / dense as f64,
            (hi - lo) / dense as f64,
        );
        let d = dist(t);
        if d < best.1 {
            best = (t, d);
        }
        if best.1 <= exact {
            break;
        }
    }
    best
}

#[cfg(test)]
#[path = "param_basin_search_tests.rs"]
mod tests;
