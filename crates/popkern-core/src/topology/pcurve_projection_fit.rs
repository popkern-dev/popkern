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

use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{NurbsCurve, Point3};

/// Samples of the historical fixed recompute (the OCCT default); the
/// initial sampling spreads this many over the segments.
pub(super) const PROJECTION_SAMPLES: usize = 23;
/// Doublings of each segment's span count: 23 samples reach 177.
const MAX_DOUBLINGS: usize = 3;
/// Fewest spans per segment (cubic interpolation needs four samples).
const MIN_SEGMENT_SPANS: usize = 3;
/// Initial gaps beyond this multiple of the target are not
/// interpolation error (a branch or topology misfit on intermediate
/// boolean bodies, measured up to 15) and are left unrefined.
const REFINABLE_GAP_FACTOR: f64 = 1.0e4;
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
    if best_gap > REFINABLE_GAP_FACTOR * target {
        return Some(best);
    }
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
        let converging = gap < 0.5 * best_gap;
        if gap < best_gap {
            (best, best_gap) = (fit, gap);
        }
        if !converging {
            break;
        }
    }
    Some(best)
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
