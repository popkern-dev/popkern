//! Root isolation inside one (curve segment, surface cell) pair, and the
//! angular-window contract for periodic conics (PK-781).
//!
//! **Isolation.** The general solver seeds Newton once per overlapping
//! pair. A pair can hold two crossings (a circle through a rod crosses
//! it twice inside one quarter-turn Bézier segment), and the tuned seed
//! finds one. Measured on a 2460-case sweep of circles and arcs through
//! rods of radius 1, 0.2 and 0.02 through tangency: with the period
//! lifts below and no isolation, 16 circle cases still returned one
//! root of two. A pair stays AMBIGUOUS while a piece of the curve
//! segment that no root explains still meets the cell's box: after a
//! root at t* the two side pieces (t* excluded by a gap) are reseeded
//! from their samples nearest the cell, recursively. A piece whose box
//! clears the cell, or whose polish diverges or falls back out of it
//! (onto the root its sibling excluded), is resolved. A pair whose
//! tuned polish found nothing is not isolated: no sweep case needed it.
//! (An earlier variant that also halved failed pairs and fallen-out
//! pieces cost 11-16% on the release M6 thread-cut serve; this one
//! measures within 2% of the solver without isolation.)
//! Depth and the polishes per pair are bounded (`MAX_DEPTH`,
//! `POLISH_BUDGET`).
//!
//! **Periodic windows.** Circles, arcs and ellipses are 2π-periodic in
//! their angle parameter, and an edge may query its curve on a window a
//! whole turn from the stored angles (a lifted seam window, a reversed
//! arc stored above 2π). The Bézier decomposition follows the stored
//! angles, so a root came back one turn from the window and was
//! filtered out: 614 of 615 shifted-arc sweep cases found nothing.
//! Of the 112 circle misses on the sweep, 96 were the surface side: a
//! cylinder crossing just below its seam polishes to u < 0. Roots now
//! lift by whole periods into the requested windows (t, u and v).

use std::f64::consts::TAU;

use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::intersection::surface_curve_newton::newton_surface_curve;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

/// A polished root: curve parameter, surface parameters, point.
pub(super) type Root = (f64, f64, f64, Point3);

/// Subdivision depth bound per pair.
const MAX_DEPTH: u32 = 4;

/// Samples per piece for its box and its nearest-sample seed.
const PIECE_SAMPLES: usize = 8;

/// The parameter period of a periodic conic (2π in its angle).
pub(super) fn curve_period(curve: &CurveGeometry) -> Option<f64> {
    matches!(
        curve,
        CurveGeometry::Circle { .. } | CurveGeometry::Arc { .. } | CurveGeometry::Ellipse { .. }
    )
    .then_some(TAU)
}

/// `x` moved into `window` (with `slack`), by whole periods when the
/// parameter is periodic; `None` when no representative lands inside.
/// A polished root's parameters are wherever Newton walked: a crossing
/// just below a cylinder's seam polishes to u = −0.1, outside [0, 2π),
/// and was dropped (measured on the PK-781 sweep).
pub(super) fn lift_into_window(
    x: f64,
    window: (f64, f64),
    slack: f64,
    period: Option<f64>,
) -> Option<f64> {
    let (w0, w1) = sorted(window);
    let inside = |x: f64| x >= w0 - slack && x <= w1 + slack;
    if inside(x) {
        return Some(x);
    }
    let period = period?;
    let lifted = x + ((w0 - slack - x) / period).ceil() * period;
    inside(lifted).then_some(lifted)
}

fn sorted((a, b): (f64, f64)) -> (f64, f64) {
    if a <= b { (a, b) } else { (b, a) }
}

/// One overlapping pair: the curve, the surface, and the cell.
pub(super) struct PairCell<'a> {
    pub curve: &'a CurveGeometry,
    pub surface: &'a SurfaceGeometry,
    pub lo: Point3,
    pub hi: Point3,
    pub u_mid: f64,
    pub v_mid: f64,
    pub margin: f64,
    pub tolerance: &'a Tolerance,
}

impl PairCell<'_> {
    fn centre(&self) -> Point3 {
        Point3::new(
            0.5 * (self.lo.x + self.hi.x),
            0.5 * (self.lo.y + self.hi.y),
            0.5 * (self.lo.z + self.hi.z),
        )
    }

    /// Whether the piece's sampled box (grown by the sample spacing's
    /// sagitta allowance) meets the cell.
    fn piece_meets(&self, (a, b): (f64, f64)) -> bool {
        let pts: Vec<Point3> = (0..=PIECE_SAMPLES)
            .map(|k| {
                self.curve
                    .point_at(a + (b - a) * k as f64 / PIECE_SAMPLES as f64)
            })
            .collect();
        let spacing = pts
            .windows(2)
            .map(|w| w[0].distance_to(&w[1]))
            .fold(0.0_f64, f64::max);
        // A sample chord of length h on a curve whose turn across it is
        // at most a quarter turn strays at most h/4 from the chord.
        let pad = self.margin + 0.25 * spacing;
        let (mut lo, mut hi) = (pts[0], pts[0]);
        for p in &pts[1..] {
            lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        lo.x - pad <= self.hi.x
            && hi.x + pad >= self.lo.x
            && lo.y - pad <= self.hi.y
            && hi.y + pad >= self.lo.y
            && lo.z - pad <= self.hi.z
            && hi.z + pad >= self.lo.z
    }

    fn nearest_sample(&self, (a, b): (f64, f64)) -> f64 {
        let c = self.centre();
        (0..=PIECE_SAMPLES)
            .map(|k| a + (b - a) * k as f64 / PIECE_SAMPLES as f64)
            .min_by(|&x, &y| {
                let dx = self.curve.point_at(x).distance_to(&c);
                let dy = self.curve.point_at(y).distance_to(&c);
                dx.total_cmp(&dy)
            })
            .unwrap_or(0.5 * (a + b))
    }
}

/// Roots of the pair beyond `first`, the pair's tuned first polish
/// (which keeps its basin: see the solver's seed notes). `piece` is the
/// curve segment's parameter window; new roots append to `found` (which
/// also holds every root already known, so a re-found root is not new).
pub(super) fn isolate_further_roots(
    pair: &PairCell<'_>,
    piece: (f64, f64),
    first: Option<f64>,
    dedup: f64,
    found: &mut Vec<Root>,
) {
    let mut budget = POLISH_BUDGET;
    split_and_seed(pair, sorted(piece), first, dedup, found, 0, &mut budget);
}

/// Extra polishes per pair, whatever the subdivision tree's shape.
const POLISH_BUDGET: u32 = 12;

fn split_and_seed(
    pair: &PairCell<'_>,
    (a, b): (f64, f64),
    root: Option<f64>,
    dedup: f64,
    found: &mut Vec<Root>,
    depth: u32,
    budget: &mut u32,
) {
    if depth >= MAX_DEPTH || b - a <= 0.0 {
        return;
    }
    // The gap excluding a found root from its side pieces: wide enough
    // that the root's own tolerance band does not keep them ambiguous.
    let gap = (1e-6 * (b - a)).max(1e-12);
    let pieces: Vec<(f64, f64)> = match root {
        Some(t) if t > a + gap && t < b - gap => vec![(a, t - gap), (t + gap, b)],
        // A root on the piece's edge (the first polish landing on a
        // segment end) explains nothing inside it: halve the piece.
        _ => {
            let m = 0.5 * (a + b);
            vec![(a, m), (m, b)]
        }
    };
    for piece in pieces {
        if *budget == 0 {
            return;
        }
        if piece.1 - piece.0 <= gap || !pair.piece_meets(piece) {
            continue;
        }
        *budget -= 1;
        let seed = pair.nearest_sample(piece);
        let solved = newton_surface_curve(
            pair.curve,
            seed,
            pair.surface,
            pair.u_mid,
            pair.v_mid,
            pair.tolerance,
        );
        match solved {
            Some(r) if r.0 > piece.0 && r.0 < piece.1 => {
                let known = found
                    .iter()
                    .any(|k| (k.0 - r.0).abs() < dedup || k.3.distance_to(&r.3) < dedup);
                if !known {
                    found.push(r);
                }
                split_and_seed(pair, piece, Some(r.0), dedup, found, depth + 1, budget);
            }
            // The polish fell out of the piece (back onto the root its
            // sibling excluded) or diverged: the nearest-sample basin
            // holds no further root.
            _ => {}
        }
    }
}

#[cfg(test)]
#[path = "surface_curve_pair_roots_tests.rs"]
mod tests;
