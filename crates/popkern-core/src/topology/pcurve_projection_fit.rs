//! Measured interpolation of projected pcurve samples (PK-555).
//!
//! A dense recompute projects samples of the edge curve into the face's
//! uv chart and interpolates them by curve FRACTION, so the pcurve runs
//! co-parameterized with its edge. Two fixed choices made that
//! interpolant unmeasured:
//! - a fixed 23 samples: a tilted drill's cylinder trim stood 7.7e-6
//!   off its edge (7.7x the model tolerance);
//! - one smooth cubic across the edge curve's parametric C0 knots: a
//!   rim written as rational quadratic arcs joined by double knots has
//!   a uv track with a derivative jump at each join; a single
//!   interpolant converges only linearly there and overshoots the
//!   trace (PK-558: 3.4e-4 off a tilted drill's exit rim, which no
//!   export rebuild within the declared accuracy could repair).
//!
//! The interpolant is therefore fitted PIECEWISE between the edge
//! curve's parametric breakpoints (C0 joins at those fractions), and
//! the sample density is doubled until the pcurve's surface image
//! follows the edge curve, as a point set, within a fraction of the
//! model's linear tolerance, the refinement stopping at the first
//! doubling that fails to halve the gap (a floor no sample count
//! moves: a vertex pinned off the curve's end).
//!
//! Uniform doubling cannot follow a track that turns fast in a short
//! stretch: a ring passing a sphere's pole within 1e-3 sweeps half a
//! turn of longitude within its distance from the pole, and 736
//! uniform samples left a rod's ring on a turned ball 1e-2 off its edge
//! (the union declined). Past the doublings, samples go where the image
//! strays (`locally_refined`): intervals straying past the target are
//! halved, the rest graded 2:1, judged co-parameterized at quarter
//! points, the better of the two fits kept (PK-1178).

use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{NurbsCurve, Point3};

/// Samples of the historical fixed recompute (the OCCT default); the
/// initial sampling spreads this many over the segments.
pub(super) const PROJECTION_SAMPLES: usize = 23;
/// Doublings of each segment's span count while the gap converges:
/// 23 samples reach 736. Three (177) left a full-turn rim on a
/// radius-10 sphere 1e-3 off its edge — cubic interpolation error
/// falls 16× per doubling, from 1.5e-2 at 23 (PK-646).
const MAX_DOUBLINGS: usize = 5;
/// Fewest spans per segment (cubic interpolation needs four samples).
const MIN_SEGMENT_SPANS: usize = 3;
/// Initial gaps beyond this multiple of the target AND beyond
/// `MISFIT_EXTENT_FRACTION` of the edge's own extent may be a branch or
/// topology misfit on intermediate boolean bodies (measured up to 15)
/// rather than interpolation error: they are refined only while each
/// doubling at least halves them, as interpolation error does (cubic:
/// 16x). A gap that size can be interpolation error: a loop passing
/// 0.05 from a sphere's pole sweeps half a turn of longitude there,
/// and 23 samples left its pcurve 4.4e-2 off the edge (PK-1055).
const REFINABLE_GAP_FACTOR: f64 = 1.0e4;
/// A gap this fraction of the edge's chord extent or below is
/// interpolation error at any model scale: the absolute floor alone
/// shipped a radius-10 sphere's full-turn rim pcurve 1.5e-2 off its
/// edge unrefined (PK-646), recorded as an import tolerance.
const MISFIT_EXTENT_FRACTION: f64 = 1.0e-2;
/// Newton steps locating the edge-curve point nearest an image sample.
const NEAREST_NEWTON_STEPS: usize = 4;
/// Fraction of the model's linear tolerance the pcurve's image may
/// stray from the edge curve.
const PROJECTION_FIT_FRACTION: f64 = 0.1;
/// uv distance within which two pieces meet at a breakpoint: both
/// interpolate the same sample, so they share it to round-off.
const JOIN_TOLERANCE: f64 = 1e-9;

/// Curve fractions of the edge curve's parametric breakpoints inside
/// `(t_lo, t_hi)`, ascending: knots of a NURBS edge curve whose
/// multiplicity reaches its degree (C0 in the parameter).
pub(super) fn breakpoint_fractions(
    curve_on_edge: &CurveOnEdge,
    (t_lo, t_hi): (f64, f64),
) -> Vec<f64> {
    let CurveGeometry::Nurbs(nurbs) = &curve_on_edge.geometry else {
        return Vec::new();
    };
    let span = t_hi - t_lo;
    if span == 0.0 || span.is_nan() {
        return Vec::new();
    }
    let mut fractions = Vec::new();
    let mut k = 0;
    while k < nurbs.knots.len() {
        let knot = nurbs.knots[k];
        let multiplicity = nurbs.knots[k..].iter().take_while(|&&x| x == knot).count();
        let f = (knot - t_lo) / span;
        if multiplicity >= nurbs.degree && f > 1e-9 && f < 1.0 - 1e-9 {
            fractions.push(f);
        }
        k += multiplicity;
    }
    fractions.sort_by(f64::total_cmp);
    fractions
}

/// Sample fractions: `spans[i]` uniform steps in segment `i` of
/// `[0, breaks..., 1]`, segment bounds shared once.
pub(super) fn segment_fractions(breaks: &[f64], spans: &[usize]) -> Vec<f64> {
    let bounds: Vec<f64> = std::iter::once(0.0)
        .chain(breaks.iter().copied())
        .chain(std::iter::once(1.0))
        .collect();
    let mut fractions = vec![0.0];
    for (pair, &n) in bounds.windows(2).zip(spans) {
        fractions.extend((1..=n).map(|i| pair[0] + (pair[1] - pair[0]) * i as f64 / n as f64));
    }
    fractions
}

/// Initial span count per segment: the historical sample budget spread
/// by segment length, at least `MIN_SEGMENT_SPANS` each.
pub(super) fn initial_spans(breaks: &[f64]) -> Vec<usize> {
    let bounds: Vec<f64> = std::iter::once(0.0)
        .chain(breaks.iter().copied())
        .chain(std::iter::once(1.0))
        .collect();
    bounds
        .windows(2)
        .map(|w| {
            let share = (PROJECTION_SAMPLES - 1) as f64 * (w[1] - w[0]);
            (share.round() as usize).max(MIN_SEGMENT_SPANS)
        })
        .collect()
}

/// Interpolate `points` (one per fraction of `segment_fractions(breaks,
/// spans)`) piecewise between the breakpoints, joined C0.
fn piecewise_interpolant(points: &[Point3], breaks: &[f64], spans: &[usize]) -> Option<NurbsCurve> {
    let fractions = segment_fractions(breaks, spans);
    if points.len() != fractions.len() {
        return None;
    }
    let mut start = 0;
    let mut joined: Option<NurbsCurve> = None;
    for &n in spans {
        let end = start + n;
        let (pts, params) = (points.get(start..=end)?, fractions.get(start..=end)?);
        let degree = 3.min(pts.len() - 1);
        let piece = crate::geometry::fitting::fit_curve(pts, params, degree, pts.len()).ok()?;
        joined = Some(match joined {
            None => piece,
            Some(left) => left.concatenate_clamped(&piece, JOIN_TOLERANCE).ok()?,
        });
        start = end;
    }
    joined
}

/// The recompute pcurve through projected samples (see the module
/// docs). `sample` projects the samples at the given curve fractions,
/// returning one uv point per fraction or fewer when a projection
/// failed. `None` when the initial sampling cannot be interpolated.
pub(super) fn measured_projection_interpolant(
    surface: &SurfaceGeometry,
    curve_on_edge: &CurveOnEdge,
    (t_lo, t_hi): (f64, f64),
    linear_tolerance: f64,
    sample: &dyn Fn(&[f64]) -> Vec<Point3>,
) -> Option<NurbsCurve> {
    let breaks = breakpoint_fractions(curve_on_edge, (t_lo, t_hi));
    let mut spans = initial_spans(&breaks);
    let points = sample(&segment_fractions(&breaks, &spans));
    let Some(mut best) = piecewise_interpolant(&points, &breaks, &spans) else {
        // A walk that dropped samples: the historical interpolant over
        // the surviving samples, unrefined.
        return uniform_interpolant(&points);
    };
    let target = PROJECTION_FIT_FRACTION * linear_tolerance;
    let mut best_gap = image_gap(&best, surface, curve_on_edge, (t_lo, t_hi), &breaks, &spans);
    let extent = edge_extent(curve_on_edge, (t_lo, t_hi));
    let misfit_sized =
        best_gap > REFINABLE_GAP_FACTOR * target && best_gap > MISFIT_EXTENT_FRACTION * extent;
    let mut best_spans = spans.clone();
    // Whether any doubling improved the fit.
    let mut improved = false;
    for _ in 0..MAX_DOUBLINGS {
        if best_gap <= target {
            break;
        }
        spans.iter_mut().for_each(|n| *n *= 2);
        let points = sample(&segment_fractions(&breaks, &spans));
        let Some(fit) = piecewise_interpolant(&points, &breaks, &spans) else {
            break;
        };
        let gap = image_gap(&fit, surface, curve_on_edge, (t_lo, t_hi), &breaks, &spans);
        // Stop only when a doubling no longer improves the fit: the
        // former "must halve" rule quit a near-tangent tee's
        // intersection pcurves at 44 spans when the steep end spans
        // improved 1.05e-3 -> 8.5e-4, serving them 1.03e-4 off their
        // edges (PK-885); two more doublings reach 2e-7.
        if gap >= best_gap || (misfit_sized && gap > 0.5 * best_gap) {
            break;
        }
        (best, best_gap, best_spans, improved) = (fit, gap, spans.clone(), true);
    }
    // Where uniform doubling cannot follow the track (a ring passing a
    // sphere's pole sweeps half a turn of longitude within its distance
    // from the pole, PK-1178), samples go where the image strays.
    // Only where doubling followed the track some way, or the chart has
    // a singular row a sweep can turn round: a gap no added samples move
    // at all on a regular chart is a misfit (a cube fillet's rails on an
    // intermediate body, 2.8e-2 off a 0.1 edge at every density, its
    // local rounds only spent five times the allocations), while beside
    // a pole doubling moves nothing until its spans reach the sweep's
    // length (a ring 1e-4 round the pole, 0.1 off at 22 spans and not
    // halved at 44, PK-1185).
    let singular =
        crate::geometry::chart::ChartDescriptor::has_singular_row(surface, linear_tolerance);
    if best_gap > target && (improved || singular) {
        let start = segment_fractions(&breaks, &best_spans);
        let range = (t_lo, t_hi);
        if let Some((fit, fractions)) = locally_refined(
            surface,
            curve_on_edge,
            range,
            &breaks,
            start,
            target,
            sample,
        ) {
            // Judged with the doubling's fit on one dense set: midpoints
            // alone let a graded cubic stray between them.
            let judge = dense_fractions(&fractions);
            let gap =
                |f: &NurbsCurve| max_gap(&interval_gaps(f, surface, curve_on_edge, range, &judge));
            if gap(&fit) < gap(&best) {
                best = fit;
            }
        }
    }
    Some(best)
}

/// Rounds of local refinement: enough to halve an interval of a
/// uniform sampling down to a sweep a millionth of the edge long and
/// then converge there.
const LOCAL_ROUNDS: usize = 40;
/// The most samples local refinement takes.
const MAX_LOCAL_SAMPLES: usize = 4096;

/// Refine `fractions` where the interpolant's image strays past
/// `target`: each round fits the samples piecewise between `breaks`,
/// judges each interval at its quarter points, halves those straying
/// within a tenth of the worst and grades the rest (no interval over
/// twice a neighbour), stopping when the rounds stall; the best fit and
/// its fractions.
fn locally_refined(
    surface: &SurfaceGeometry,
    curve_on_edge: &CurveOnEdge,
    range: (f64, f64),
    breaks: &[f64],
    mut fractions: Vec<f64>,
    target: f64,
    sample: &dyn Fn(&[f64]) -> Vec<Point3>,
) -> Option<(NurbsCurve, Vec<f64>)> {
    let periods = crate::geometry::chart::ChartDescriptor::evaluation_periods_of(surface);
    let mut best: Option<(NurbsCurve, f64, Vec<f64>)> = None;
    let mut history: Vec<f64> = Vec::new();
    for _ in 0..LOCAL_ROUNDS {
        let points = sample(&fractions);
        let fit = interpolant_at(&points, &fractions, breaks)?;
        let judge = dense_fractions(&fractions);
        let gaps = interval_gaps(&fit, surface, curve_on_edge, range, &judge);
        let interval_gap: Vec<f64> = gaps.chunks(4).map(max_gap).collect();
        let worst = max_gap(&interval_gap);
        // An interval past the target still turning the chart past a
        // quarter period is a sweep round a singular row its samples have
        // not reached: the image cannot settle there, or ringing beside
        // it, until the interval is as short as the sweep (a ring 1e-4
        // round a pole on an edge sampled at 0.14, ten halvings), so the
        // stall rule waits for it.
        let sweeping = interval_gap.iter().zip(points.windows(2)).any(|(&g, w)| {
            g > target && turns_past_quarter(periods, (w[0].x, w[0].y), (w[1].x, w[1].y))
        });
        let stalled = best.as_ref().is_some_and(|b| worst >= b.1)
            || history.len() >= 2 && worst > 0.5 * history[history.len() - 2];
        if stalled && !sweeping {
            break;
        }
        history.push(worst);
        // The intervals straying within a tenth of the worst: splitting
        // every one past the target at first doubled the whole edge,
        // which leaves a pole sweep where it was (a ring 1e-4 round the
        // pole: 6.8e-2 at 23 samples, 6.1e-2 at 45, PK-1185).
        let cut = target.max(0.1 * worst);
        if best.as_ref().is_none_or(|b| worst < b.1) {
            best = Some((fit, worst, fractions.clone()));
        }
        if worst <= target || fractions.len() >= MAX_LOCAL_SAMPLES {
            break;
        }
        let mut next = Vec::with_capacity(2 * fractions.len());
        for (w, &g) in fractions.windows(2).zip(&interval_gap) {
            next.push(w[0]);
            if g > cut {
                next.push(0.5 * (w[0] + w[1]));
            }
        }
        next.extend(fractions.last().copied());
        fractions = graded(next);
    }
    best.map(|(fit, _, fractions)| (fit, fractions))
}

/// Whether the chart step from `a` to `b` turns a periodic coordinate
/// (periods `periods`) past a quarter period.
fn turns_past_quarter(periods: (Option<f64>, Option<f64>), a: (f64, f64), b: (f64, f64)) -> bool {
    [(b.0 - a.0, periods.0), (b.1 - a.1, periods.1)]
        .iter()
        .any(|&(d, p)| p.is_some_and(|p| p > 0.0 && d.abs() > 0.25 * p))
}

/// `fractions` with every interval more than twice a neighbour's
/// halved until none is: a cubic through samples graded more steeply
/// rings between them.
fn graded(mut fractions: Vec<f64>) -> Vec<f64> {
    loop {
        let len: Vec<f64> = fractions.windows(2).map(|w| w[1] - w[0]).collect();
        let wide: Vec<usize> = (0..len.len())
            .filter(|&i| {
                let left = i.checked_sub(1).map(|j| len[j]);
                let right = len.get(i + 1).copied();
                [left, right]
                    .into_iter()
                    .flatten()
                    .any(|n| len[i] > 2.0 * n)
            })
            .collect();
        if wide.is_empty() || fractions.len() >= MAX_LOCAL_SAMPLES {
            return fractions;
        }
        for &i in wide.iter().rev() {
            let mid = 0.5 * (fractions[i] + fractions[i + 1]);
            fractions.insert(i + 1, mid);
        }
    }
}

/// Each interval of `fractions` at its start and three quarter points,
/// then the last fraction: gaps over the result's intervals come four
/// to an interval of `fractions`.
fn dense_fractions(fractions: &[f64]) -> Vec<f64> {
    let mut out = Vec::with_capacity(4 * fractions.len());
    for w in fractions.windows(2) {
        for k in 0..4 {
            out.push(w[0] + (w[1] - w[0]) * f64::from(k) / 4.0);
        }
    }
    out.extend(fractions.last().copied());
    out
}

/// The largest of `gaps` (nought for none).
fn max_gap(gaps: &[f64]) -> f64 {
    gaps.iter().copied().fold(0.0_f64, f64::max)
}

/// Interpolate `points` at `fractions` (ascending, holding every
/// breakpoint) piecewise between `breaks`, joined C0.
fn interpolant_at(points: &[Point3], fractions: &[f64], breaks: &[f64]) -> Option<NurbsCurve> {
    if points.len() != fractions.len() {
        return None;
    }
    let mut cuts: Vec<usize> = vec![0];
    for &b in breaks {
        let k = fractions.iter().position(|&f| (f - b).abs() <= 1e-12)?;
        cuts.push(k);
    }
    cuts.push(fractions.len() - 1);
    let mut joined: Option<NurbsCurve> = None;
    for w in cuts.windows(2) {
        let (pts, params) = (points.get(w[0]..=w[1])?, fractions.get(w[0]..=w[1])?);
        if pts.len() < 2 {
            return None;
        }
        let degree = 3.min(pts.len() - 1);
        let piece = crate::geometry::fitting::fit_curve(pts, params, degree, pts.len()).ok()?;
        joined = Some(match joined {
            None => piece,
            Some(left) => left.concatenate_clamped(&piece, JOIN_TOLERANCE).ok()?,
        });
    }
    joined
}

/// Per interval of `fractions`, the distance from the image of `fit`
/// at its midpoint to the edge curve at the same fraction (the
/// agreeing direction of the two): co-parameterized, as the mesher
/// pairs a canonical station's point with the pcurve at its fraction
/// (a fit on the curve as a point set but drifting in parameter failed
/// that round trip, its loop fell to private samples and the measured
/// volume lost 2e-5).
fn interval_gaps(
    fit: &NurbsCurve,
    surface: &SurfaceGeometry,
    curve_on_edge: &CurveOnEdge,
    (t_lo, t_hi): (f64, f64),
    fractions: &[f64],
) -> Vec<f64> {
    let curve = &curve_on_edge.geometry;
    let along = |reversed: bool| -> Vec<f64> {
        fractions
            .windows(2)
            .map(|w| {
                let f = 0.5 * (w[0] + w[1]);
                let uv = fit.point_at(f);
                let image = surface.point_at(uv.x, uv.y);
                let g = if reversed { 1.0 - f } else { f };
                curve.point_at(t_lo + (t_hi - t_lo) * g).distance_to(&image)
            })
            .collect()
    };
    let (a, b) = (along(false), along(true));
    if max_gap(&a) <= max_gap(&b) { a } else { b }
}

/// The edge curve's polygonal length over 32 chords.
fn edge_extent(curve_on_edge: &CurveOnEdge, (t_lo, t_hi): (f64, f64)) -> f64 {
    let at = |k: usize| {
        curve_on_edge
            .geometry
            .point_at(t_lo + (t_hi - t_lo) * k as f64 / 32.0)
    };
    (0..32).map(|k| at(k).distance_to(&at(k + 1))).sum()
}

/// The historical interpolant: uniform parameters over the samples.
fn uniform_interpolant(points: &[Point3]) -> Option<NurbsCurve> {
    if points.len() < 4 {
        return None;
    }
    let params: Vec<f64> = (0..points.len())
        .map(|i| i as f64 / (points.len() - 1) as f64)
        .collect();
    crate::geometry::fitting::fit_curve(points, &params, 3, points.len()).ok()
}

/// Largest distance from the pcurve's surface image to the edge curve
/// as a point set, at every mid-sample fraction: the nearest curve
/// point is found by Newton steps seeded at the same fraction (the
/// interpolant's parameterization may drift from the curve's near a
/// pinned vertex without moving its image off the curve).
fn image_gap(
    fit: &NurbsCurve,
    surface: &SurfaceGeometry,
    curve_on_edge: &CurveOnEdge,
    (t_lo, t_hi): (f64, f64),
    breaks: &[f64],
    spans: &[usize],
) -> f64 {
    let curve = &curve_on_edge.geometry;
    let (lo, hi) = (t_lo.min(t_hi), t_lo.max(t_hi));
    let fractions = segment_fractions(breaks, spans);
    // The walk's samples may run either way along the curve's
    // parameter: measure both correspondences, keep the agreeing one.
    let gap_along = |reversed: bool| {
        fractions
            .windows(2)
            .map(|w| {
                let f = 0.5 * (w[0] + w[1]);
                let uv = fit.point_at(f);
                let image = surface.point_at(uv.x, uv.y);
                let g = if reversed { 1.0 - f } else { f };
                let mut t = t_lo + (t_hi - t_lo) * g;
                for _ in 0..NEAREST_NEWTON_STEPS {
                    let d = curve.derivative_at(t, 1);
                    let speed = d.dot(&d);
                    if speed.is_nan() || speed <= 0.0 {
                        break;
                    }
                    t = (t - (curve.point_at(t) - image).dot(&d) / speed).clamp(lo, hi);
                }
                curve.point_at(t).distance_to(&image)
            })
            .fold(0.0_f64, f64::max)
    };
    gap_along(false).min(gap_along(true))
}

#[cfg(test)]
#[path = "pcurve_projection_fit_tests.rs"]
mod tests;
