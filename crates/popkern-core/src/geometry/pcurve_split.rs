//! Image-honest pcurve splitting.
//!
//! Splitting a pcurve by a parameter FRACTION trusts that its
//! parameter range agrees with the edge curve's — a lie whenever an
//! earlier surgery re-ranged one of them (PK-391). The image is the
//! truth: this module locates the pcurve parameter whose image on the
//! host surface is nearest the 3D split point and cuts there.
//!
//! Unlike the topology-internal split helper (which image-searches
//! only NURBS pcurves and falls back to the fraction), this public
//! entry point image-searches EVERY variant and reports failure
//! instead of guessing: a seam iso-line's split must land exactly at
//! the split node or the chart-branch bookkeeping downstream reads
//! the wrong seam side (PK-441: the boolean graph split dropped seam
//! pcurves entirely, the chart embedding re-projected them onto an
//! arbitrary branch of the closed chart, and the walk substrate
//! declined on a degree imbalance).

use crate::geometry::Point3;
use crate::geometry::pcurve::{Pcurve, PcurveGeometry};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

impl Pcurve {
    /// Split at normalised fraction `t ∈ (0, 1)` of the parameter
    /// range, returning `(first, second)` in parameter order. Every
    /// analytic variant splits closed-form; a NURBS pcurve splits by
    /// knot insertion. Returns `None` when the NURBS split fails or
    /// the parameter span is degenerate.
    pub fn split_at_fraction(&self, t: f64) -> Option<(Pcurve, Pcurve)> {
        let (p_start, p_end) = (self.t_start, self.t_end);
        if (p_end - p_start).abs() <= 0.0 {
            return None;
        }
        let split_param = p_start + t * (p_end - p_start);
        match &self.geometry {
            PcurveGeometry::Linear { start, end } => {
                let mid = (
                    start.0 + t * (end.0 - start.0),
                    start.1 + t * (end.1 - start.1),
                );
                Some((
                    Pcurve::linear(*start, mid, p_start, split_param),
                    Pcurve::linear(mid, *end, split_param, p_end),
                ))
            }
            PcurveGeometry::Elliptical {
                center,
                u_radius,
                v_radius,
                rotation,
                start_angle,
                end_angle,
            } => {
                // Angle interpolates linearly with normalised parameter.
                let mid_angle = start_angle + t * (end_angle - start_angle);
                Some((
                    Pcurve::elliptical(
                        *center,
                        *u_radius,
                        *v_radius,
                        *rotation,
                        *start_angle,
                        mid_angle,
                        p_start,
                        split_param,
                    ),
                    Pcurve::elliptical(
                        *center,
                        *u_radius,
                        *v_radius,
                        *rotation,
                        mid_angle,
                        *end_angle,
                        split_param,
                        p_end,
                    ),
                ))
            }
            PcurveGeometry::Sinusoidal {
                u_offset,
                v_const,
                v_cos,
                v_sin,
            } => Some((
                Pcurve::sinusoidal(*u_offset, *v_const, *v_cos, *v_sin, p_start, split_param),
                Pcurve::sinusoidal(*u_offset, *v_const, *v_cos, *v_sin, split_param, p_end),
            )),
            PcurveGeometry::Nurbs(curve) => match curve.split_at(split_param) {
                Ok((left, right)) => Some((Pcurve::from_nurbs(left), Pcurve::from_nurbs(right))),
                Err(_) => None,
            },
            PcurveGeometry::CylinderArcViaRationalQuad {
                theta_a,
                theta_b,
                v_const,
            } => {
                // The angular formula already handles arbitrary
                // parameter sub-ranges; only the windows change.
                let piece = |a: f64, b: f64| Pcurve {
                    geometry: PcurveGeometry::CylinderArcViaRationalQuad {
                        theta_a: *theta_a,
                        theta_b: *theta_b,
                        v_const: *v_const,
                    },
                    t_start: a,
                    t_end: b,
                };
                Some((piece(p_start, split_param), piece(split_param, p_end)))
            }
        }
    }
}

/// Split `pcurve` at the parameter whose image on `surface` is
/// nearest `split_pos`.
///
/// Returns the `(first, second)` pieces in the pcurve's own parameter
/// order, or `None` when the image search does not reach
/// `split_pos` within `10 × tolerance` or the cut would be an
/// endpoint (callers keep the unsplit pcurve, or drop it, on their
/// own policy).
///
/// The image search is a coarse scan refined by golden section. A
/// pcurve whose image WINDS many times (a thread tool's root edge on
/// its swept chart: 13.9 turns of helix) has one distance basin per
/// winding, so the coarse grid must resolve every basin and the
/// refinement must not trust only the globally-nearest SAMPLE: on the
/// M2.5 × 6 rod (PK-669) 64 samples spaced 78° along a 0.98-radius
/// helix put an adjacent turn's sample (0.45 = one pitch away) nearer
/// than the true turn's (up to 0.65 away), the refinement converged
/// one pitch off, the split declined, both pieces lost their pcurves,
/// and the 16-sample projection fallback could not follow 13 turns —
/// the root face never split and the boolean declined. The scan is
/// therefore densified per chart span ([`image_scan_count`]) and
/// every local minimum of it is a refinement candidate, nearest
/// first. The global-best candidate is still tried first, so a split
/// whose nearest sample already sat in the right basin resolves in
/// that same basin.
pub fn split_pcurve_at_image(
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    split_pos: &Point3,
    tolerance: f64,
) -> Option<(Pcurve, Pcurve)> {
    let (p_start, p_end) = (pcurve.t_start, pcurve.t_end);
    let span = p_end - p_start;
    if span.abs() <= 0.0 {
        return None;
    }
    let image_distance = |f: f64| {
        let uv = pcurve.point_at(p_start + f * span);
        surface.point_at(uv.x, uv.y).distance_to(split_pos)
    };
    let count = image_scan_count(pcurve, surface);
    let scan: Vec<f64> = (0..=count)
        .map(|i| image_distance(i as f64 / count as f64))
        .collect();
    let window = 1.0 / count as f64;
    let refine = |center: f64| {
        let (mut lo, mut hi) = (
            (center - window).clamp(0.0, 1.0),
            (center + window).clamp(0.0, 1.0),
        );
        for _ in 0..48 {
            let m1 = lo + (hi - lo) * 0.382;
            let m2 = lo + (hi - lo) * 0.618;
            if image_distance(m1) < image_distance(m2) {
                hi = m2;
            } else {
                lo = m1;
            }
        }
        0.5 * (lo + hi)
    };
    // Bounded: the candidates are ordered nearest-first, so a decline
    // costs at most MAX_BASINS refinements however noisy the scan.
    const MAX_BASINS: usize = 64;
    for index in scan_basins(&scan).into_iter().take(MAX_BASINS) {
        let refined = refine(index as f64 / count as f64);
        if image_distance(refined) > tolerance * 10.0 {
            continue;
        }
        if refined <= 1e-9 || refined >= 1.0 - 1e-9 {
            return None;
        }
        return pcurve.split_at_fraction(refined);
    }
    None
}

/// Coarse scan sample count for an image search: 64, raised to two
/// samples per chart span the pcurve's image crosses — knot spans of
/// a NURBS host (a rational helix sweep carries one per quarter
/// turn), eight per period of a periodic analytic host — so every
/// winding's distance basin is resolved. Capped at 4096.
fn image_scan_count(pcurve: &Pcurve, surface: &SurfaceGeometry) -> usize {
    const BASE: usize = 64;
    let (mut u_lo, mut u_hi, mut v_lo, mut v_hi) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for i in 0..=BASE {
        let uv = pcurve
            .point_at(pcurve.t_start + (pcurve.t_end - pcurve.t_start) * i as f64 / BASE as f64);
        u_lo = u_lo.min(uv.x);
        u_hi = u_hi.max(uv.x);
        v_lo = v_lo.min(uv.y);
        v_hi = v_hi.max(uv.y);
    }
    let spans_in = |knots: &[f64], lo: f64, hi: f64| {
        let mut count = 0usize;
        let mut last = f64::NEG_INFINITY;
        for &k in knots {
            if k > lo && k < hi && k > last {
                count += 1;
                last = k;
            }
        }
        count + 1
    };
    let wanted = match surface {
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
    wanted.clamp(BASE, 4096)
}

/// Refinement candidates of a coarse distance scan: the global
/// nearest sample first (the historical single candidate), then every
/// other local minimum in increasing sampled distance.
fn scan_basins(scan: &[f64]) -> Vec<usize> {
    let n = scan.len();
    let mut minima: Vec<usize> = (0..n)
        .filter(|&i| {
            let left = i == 0 || scan[i] <= scan[i - 1];
            let right = i + 1 == n || scan[i] <= scan[i + 1];
            left && right
        })
        .collect();
    minima.sort_by(|&a, &b| scan[a].total_cmp(&scan[b]).then(a.cmp(&b)));
    minima
}

#[cfg(test)]
#[path = "pcurve_split_tests.rs"]
mod pcurve_split_tests;
