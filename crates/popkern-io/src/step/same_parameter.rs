//! SameParameter rebuild of foreign pcurves on import (PK-552).
//!
//! popkern pcurves run co-parameterized with their edge curve: at a
//! given fraction of the edge, the pcurve's surface image and the edge
//! curve coincide. Foreign writers do not all keep that. OpenCASCADE
//! writes a plane circle's pcurve as the whole rational circle: longer
//! than the edge, and not parameterized by angle. A pcurve whose
//! surface image lies on the edge curve but fails the same-fraction
//! check is rebuilt on the edge's own parameterization:
//! - each edge sample is projected onto the pcurve's surface image,
//!   continuing from the previous sample (wrapping across a closed
//!   pcurve's seam, unwrapping a periodic surface's u or v);
//! - the projected uv track is interpolated against the edge
//!   parameter, doubling the span count until the rebuilt pcurve's
//!   image is within tolerance of the edge curve at every probe.
//!
//! This also trims a longer or closed pcurve to the edge's vertex range.
//! A pcurve whose image leaves the edge curve by more than the
//! tolerance is not repaired; the caller keeps its original decision.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::fitting::interpolate_curve;
use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceOnFace};

use crate::step::diagnostics::{DiagnosticLog, Phase};

/// Samples of the pcurve image searched for the first edge point.
const IMAGE_SAMPLES: usize = 129;
/// Samples in each continuation window.
const WINDOW_SAMPLES: usize = 33;
/// Continuation window half-width as a fraction of the pcurve domain.
const WINDOW_FRACTION: f64 = 0.125;
/// Widening factor and limit (in domain widths) of the search along an
/// unbounded linear pcurve, whose file domain need not cover the edge.
const LINE_WIDENING: f64 = 4.0;
const LINE_MAX_WIDTHS: f64 = 1.0e9;
/// Golden-section refinement steps (bracket shrinks by 0.618 each).
const GOLDEN_STEPS: usize = 64;
/// Initial and maximum edge spans of the rebuilt pcurve.
const INITIAL_SPANS: usize = 16;
const MAX_SPANS: usize = 256;
/// Degree of the rebuilt pcurve.
const FIT_DEGREE: usize = 3;
/// Fraction of the tolerance the rebuilt pcurve aims for, leaving
/// the rest of the tolerance to the file's own geometric gap.
const FIT_TARGET_FRACTION: f64 = 0.1;

/// A rebuilt pcurve and its largest distance from the edge curve at the
/// same parameter (nodes and mid-span probes).
pub(crate) struct SameParameterRepair {
    pub pcurve: Pcurve,
    pub max_gap: f64,
}

/// Rebuild `candidate` on the edge's parameterization in half-edge
/// direction, or `None` when its surface image does not follow the
/// edge curve within `tolerance`.
pub(crate) fn same_parameter_repair(
    candidate: &Pcurve,
    curve: &CurveOnEdge,
    surface: &SurfaceOnFace,
    edge_forward: bool,
    tolerance: f64,
) -> Option<SameParameterRepair> {
    let track = PcurveTrack::new(candidate, surface, tolerance)?;
    let (e0, e1) = if edge_forward {
        (curve.t_start, curve.t_end)
    } else {
        (curve.t_end, curve.t_start)
    };
    let edge_point = |t: f64| curve.geometry.point_at(t);
    let mut spans = INITIAL_SPANS;
    let mut best: Option<SameParameterRepair> = None;
    loop {
        let nodes = track.follow(e0, e1, spans, &edge_point)?;
        let repair = fit_nodes(&nodes, e0, e1, surface, &edge_point)?;
        let done = repair.max_gap <= FIT_TARGET_FRACTION * tolerance || spans >= MAX_SPANS;
        if best.as_ref().is_none_or(|b| repair.max_gap < b.max_gap) {
            best = Some(repair);
        }
        if done {
            return best.filter(|b| b.max_gap <= tolerance);
        }
        spans *= 2;
    }
}

/// Counts the repairs made while importing one body, for one summary
/// diagnostic instead of one per half-edge.
#[derive(Debug, Default)]
pub(crate) struct SameParameterTally {
    repaired: usize,
    max_gap: f64,
}

impl SameParameterTally {
    pub(crate) fn record(&mut self, max_gap: f64) {
        self.repaired += 1;
        self.max_gap = self.max_gap.max(max_gap);
    }

    pub(crate) fn report(&self, log: &mut DiagnosticLog) {
        if self.repaired > 0 {
            log.warn(
                Phase::MapTopology,
                None,
                format!(
                    "rebuilt {} STEP pcurve(s) on their edge curve's parameterization \
                     (trimmed / same-parameter; max residual {:.3e})",
                    self.repaired, self.max_gap
                ),
            );
        }
    }
}

/// A pcurve's surface image, searchable by projection.
struct PcurveTrack<'a> {
    pcurve: &'a Pcurve,
    surface: &'a SurfaceOnFace,
    s0: f64,
    s1: f64,
    closed: bool,
    /// A linear pcurve: an infinite line in uv (a STEP LINE resolves to
    /// a unit segment of it), evaluated without clamping.
    unbounded: bool,
    tolerance: f64,
}

/// One node of the followed track: edge parameter and surface uv.
struct TrackNode {
    t: f64,
    uv: (f64, f64),
}

impl<'a> PcurveTrack<'a> {
    fn new(pcurve: &'a Pcurve, surface: &'a SurfaceOnFace, tolerance: f64) -> Option<Self> {
        let (s0, s1) = if pcurve.t_start <= pcurve.t_end {
            (pcurve.t_start, pcurve.t_end)
        } else {
            (pcurve.t_end, pcurve.t_start)
        };
        let span = s1 - s0;
        if span.is_nan() || span <= 0.0 {
            return None;
        }
        let mut track = Self {
            pcurve,
            surface,
            s0,
            s1,
            closed: false,
            unbounded: matches!(pcurve.geometry, PcurveGeometry::Linear { .. }),
            tolerance,
        };
        track.closed = track.image(s0).distance_to(&track.image(s1)) <= tolerance;
        Some(track)
    }

    /// uv at domain parameter `s` (wrapped into the domain when closed).
    fn uv(&self, s: f64) -> Point3 {
        let span = self.s1 - self.s0;
        let s = if self.unbounded {
            s
        } else if self.closed {
            self.s0 + (s - self.s0).rem_euclid(span)
        } else {
            s.clamp(self.s0, self.s1)
        };
        self.pcurve.point_at(s)
    }

    fn image(&self, s: f64) -> Point3 {
        let uv = self.uv(s);
        self.surface.geometry.point_at(uv.x, uv.y)
    }

    /// The parameter in `[lo, hi]` whose image is nearest `target`.
    fn nearest(&self, target: Point3, lo: f64, hi: f64, samples: usize) -> (f64, f64) {
        let (lo, hi) = if self.closed || self.unbounded {
            (lo, hi)
        } else {
            (lo.max(self.s0), hi.min(self.s1))
        };
        let step = (hi - lo) / (samples - 1) as f64;
        let dist = |s: f64| self.image(s).distance_to(&target);
        let best = (0..samples)
            .map(|k| lo + step * k as f64)
            .min_by(|a, b| dist(*a).total_cmp(&dist(*b)))
            .unwrap_or(lo);
        let (mut a, mut b) = ((best - step).max(lo), (best + step).min(hi));
        let ratio = 0.5 * (5.0_f64.sqrt() - 1.0);
        for _ in 0..GOLDEN_STEPS {
            let (c, d) = (b - ratio * (b - a), a + ratio * (b - a));
            if dist(c) < dist(d) {
                b = d;
            } else {
                a = c;
            }
        }
        let s = 0.5 * (a + b);
        (s, dist(s))
    }

    /// Speed of the image at `s`: `|d image / ds|` (central difference).
    fn speed(&self, s: f64) -> f64 {
        let h = 1e-6 * (self.s1 - self.s0).max(1e-12);
        self.image(s + h).distance_to(&self.image(s - h)) / (2.0 * h)
    }

    /// The nearest image point over the whole domain; along an
    /// unbounded line the search widens until the target is on the
    /// image or the widening limit is reached.
    fn nearest_anywhere(&self, target: Point3) -> (f64, f64) {
        let width = self.s1 - self.s0;
        let (mut lo, mut hi) = (self.s0, self.s1);
        loop {
            let found = self.nearest(target, lo, hi, IMAGE_SAMPLES);
            if !self.unbounded || found.1 <= self.tolerance || hi - lo > LINE_MAX_WIDTHS * width {
                return found;
            }
            let grow = 0.5 * (LINE_WIDENING - 1.0) * (hi - lo);
            (lo, hi) = (lo - grow, hi + grow);
        }
    }

    /// Project `spans + 1` edge samples onto the image in order.
    /// Declines when a sample is off the image or the projection
    /// reverses direction.
    fn follow(
        &self,
        e0: f64,
        e1: f64,
        spans: usize,
        edge_point: &dyn Fn(f64) -> Point3,
    ) -> Option<Vec<TrackNode>> {
        let fixed_window = WINDOW_FRACTION * (self.s1 - self.s0);
        let mut nodes: Vec<TrackNode> = Vec::with_capacity(spans + 1);
        let mut previous: Option<(f64, Point3)> = None;
        let mut direction = 0.0_f64;
        for i in 0..=spans {
            let t = e0 + (e1 - e0) * i as f64 / spans as f64;
            let target = edge_point(t);
            let (s, gap) = match previous {
                Some((p, last_target)) => {
                    // Continue from the previous node: the window covers
                    // twice the parameter step the image speed implies,
                    // so a periodic image (a line wrapping a cylinder)
                    // cannot jump to another copy.
                    let step_window = 2.0 * target.distance_to(&last_target) / self.speed(p);
                    let window = if self.unbounded || !step_window.is_finite() {
                        step_window
                    } else {
                        step_window.max(fixed_window)
                    };
                    if window.is_finite() && window > 0.0 {
                        self.nearest(target, p - window, p + window, WINDOW_SAMPLES)
                    } else {
                        self.nearest_anywhere(target)
                    }
                }
                None => self.nearest_anywhere(target),
            };
            if gap > self.tolerance {
                return None;
            }
            if let Some((p, _)) = previous {
                let step = s - p;
                if direction * step < 0.0 {
                    return None;
                }
                if step.abs() > 0.0 {
                    direction = step.signum();
                }
            }
            let uv = self.uv(s);
            let uv = match nodes.last() {
                Some(last) => unwrap_periodic(self.surface, last.uv, (uv.x, uv.y)),
                None => (uv.x, uv.y),
            };
            nodes.push(TrackNode { t, uv });
            previous = Some((s, target));
        }
        Some(nodes)
    }
}

/// `uv` shifted by whole surface periods to lie nearest `reference`.
fn unwrap_periodic(surface: &SurfaceOnFace, reference: (f64, f64), uv: (f64, f64)) -> (f64, f64) {
    let (pu, pv) = surface.geometry.parameter_periodicity();
    let near = |r: f64, x: f64, period: Option<f64>| match period {
        Some(p) if p > 0.0 => x + p * ((r - x) / p).round(),
        _ => x,
    };
    (near(reference.0, uv.0, pu), near(reference.1, uv.1, pv))
}

/// Interpolate the nodes against the edge parameter and measure the
/// rebuilt pcurve's image against the edge curve at nodes and mid-spans.
fn fit_nodes(
    nodes: &[TrackNode],
    e0: f64,
    e1: f64,
    surface: &SurfaceOnFace,
    edge_point: &dyn Fn(f64) -> Point3,
) -> Option<SameParameterRepair> {
    let mut ordered: Vec<&TrackNode> = nodes.iter().collect();
    if e0 > e1 {
        ordered.reverse();
    }
    let params: Vec<f64> = ordered.iter().map(|n| n.t).collect();
    let points: Vec<Point3> = ordered
        .iter()
        .map(|n| Point3::new(n.uv.0, n.uv.1, 0.0))
        .collect();
    let fit = interpolate_curve(&points, &params, FIT_DEGREE).ok()?;
    let gap_at = |t: f64| {
        let uv = fit.point_at(t);
        surface
            .geometry
            .point_at(uv.x, uv.y)
            .distance_to(&edge_point(t))
    };
    let mut max_gap = 0.0_f64;
    for pair in params.windows(2) {
        max_gap = max_gap
            .max(gap_at(pair[0]))
            .max(gap_at(0.5 * (pair[0] + pair[1])));
    }
    if let Some(&last) = params.last() {
        max_gap = max_gap.max(gap_at(last));
    }
    Some(SameParameterRepair {
        pcurve: Pcurve {
            geometry: PcurveGeometry::Nurbs(fit),
            t_start: e0,
            t_end: e1,
        },
        max_gap,
    })
}

#[cfg(test)]
#[path = "same_parameter_tests.rs"]
mod tests;
