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
    const COARSE: usize = 64;
    let mut best = (0.0, image_distance(0.0));
    for i in 1..=COARSE {
        let f = i as f64 / COARSE as f64;
        let d = image_distance(f);
        if d < best.1 {
            best = (f, d);
        }
    }
    let window = 1.0 / COARSE as f64;
    let (mut lo, mut hi) = (
        (best.0 - window).clamp(0.0, 1.0),
        (best.0 + window).clamp(0.0, 1.0),
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
    let refined = 0.5 * (lo + hi);
    if image_distance(refined) > tolerance * 10.0 {
        return None;
    }
    if refined <= 1e-9 || refined >= 1.0 - 1e-9 {
        return None;
    }
    pcurve.split_at_fraction(refined)
}

#[cfg(test)]
#[path = "pcurve_split_tests.rs"]
mod pcurve_split_tests;
