//! Local helpers for split.rs — kept in a sibling file to fit the
//! 500-line per-file structural cap. Used only by `split.rs` and
//! `split_boundary.rs` siblings; not part of the public API.

use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::pcurve::{Pcurve, PcurveGeometry};

/// Split a pcurve at parameter `t` (in [0,1] of the edge).
///
/// Returns `(left, right)` where `left` covers the first portion and
/// `right` covers the second. If the pcurve is `None` or splitting
/// fails, returns `(None, None)`. Handles both `Linear` and `Nurbs`
/// pcurve variants.
/// Split a pcurve at the parameter whose IMAGE on `surface` is nearest
/// `split_pos`, falling back to the parameter fraction `t` when the
/// pcurve cannot be evaluated against a surface or its image never
/// comes near the point.
///
/// A fraction split trusts that the pcurve's parameter range agrees
/// with the edge curve's. When an earlier surgery shortened the edge
/// without trimming its pcurve (measured: a strip rail's side-plane
/// pcurve still spanned the full plate height while the edge ended at
/// the rim fillet's rail), the fraction lands the cut 25 % past the
/// split point and the left piece's image overshoots the edge — every
/// distance-to-curve consumer then reads a gap of up to 0.7 r (PK-391).
/// The image is the truth; the parameter is derived from it.
pub(crate) fn split_pcurve_at_point(
    pcurve: &Option<Pcurve>,
    t: f64,
    surface: Option<&crate::geometry::surface::SurfaceGeometry>,
    split_pos: &crate::geometry::Point3,
    tolerance: f64,
) -> (Option<Pcurve>, Option<Pcurve>) {
    use crate::geometry::surface::SurfaceEval;
    let (Some(pc), Some(surface)) = (pcurve, surface) else {
        return split_pcurve_at(pcurve, t);
    };
    if !matches!(pc.geometry, PcurveGeometry::Nurbs(_)) {
        return split_pcurve_at(pcurve, t);
    }
    let (p_start, p_end) = (pc.t_start, pc.t_end);
    let span = p_end - p_start;
    if span.abs() <= 0.0 {
        return split_pcurve_at(pcurve, t);
    }
    let image_distance = |f: f64| {
        let uv = pc.point_at(p_start + f * span);
        surface.point_at(uv.x, uv.y).distance_to(split_pos)
    };
    const COARSE: usize = 64;
    let mut best = (t, image_distance(t));
    for i in 0..=COARSE {
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
    for _ in 0..40 {
        let m1 = lo + (hi - lo) * 0.382;
        let m2 = lo + (hi - lo) * 0.618;
        if image_distance(m1) < image_distance(m2) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let refined = 0.5 * (lo + hi);
    let refined_distance = image_distance(refined);
    // Take the image-honest parameter only when the image actually
    // reaches the split point and it is an interior cut; otherwise the
    // pcurve is not a faithful image of this edge at all and the
    // fraction is the least surprising choice.
    if refined_distance <= tolerance * 10.0 && refined > 1e-9 && refined < 1.0 - 1e-9 {
        split_pcurve_at(pcurve, refined)
    } else {
        split_pcurve_at(pcurve, t)
    }
}

pub(crate) fn split_pcurve_at(pcurve: &Option<Pcurve>, t: f64) -> (Option<Pcurve>, Option<Pcurve>) {
    let Some(pc) = pcurve else {
        return (None, None);
    };
    let (p_start, p_end) = (pc.t_start, pc.t_end);
    let split_param = p_start + t * (p_end - p_start);
    match &pc.geometry {
        PcurveGeometry::Linear { start, end } => {
            let span = p_end - p_start;
            let s = if span.abs() > 0.0 { t } else { 0.0 };
            let mid = (
                start.0 + s * (end.0 - start.0),
                start.1 + s * (end.1 - start.1),
            );
            (
                Some(Pcurve::linear(*start, mid, p_start, split_param)),
                Some(Pcurve::linear(mid, *end, split_param, p_end)),
            )
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
            let left = Pcurve::elliptical(
                *center,
                *u_radius,
                *v_radius,
                *rotation,
                *start_angle,
                mid_angle,
                p_start,
                split_param,
            );
            let right = Pcurve::elliptical(
                *center,
                *u_radius,
                *v_radius,
                *rotation,
                mid_angle,
                *end_angle,
                split_param,
                p_end,
            );
            (Some(left), Some(right))
        }
        PcurveGeometry::Sinusoidal {
            u_offset,
            v_const,
            v_cos,
            v_sin,
        } => {
            // Sinusoidal form is closed-form on any sub-range —
            // restrict the parameter window without changing the
            // analytic formula.
            let left =
                Pcurve::sinusoidal(*u_offset, *v_const, *v_cos, *v_sin, p_start, split_param);
            let right = Pcurve::sinusoidal(*u_offset, *v_const, *v_cos, *v_sin, split_param, p_end);
            (Some(left), Some(right))
        }
        PcurveGeometry::Nurbs(curve) => match curve.split_at(split_param) {
            Ok((left, right)) => (
                Some(Pcurve::from_nurbs(left)),
                Some(Pcurve::from_nurbs(right)),
            ),
            Err(_) => (None, None),
        },
        PcurveGeometry::CylinderArcViaRationalQuad {
            theta_a,
            theta_b,
            v_const,
        } => {
            // Splitting the rational quadratic at any sub-range
            // produces two pieces of the SAME analytic curve with
            // unchanged angular span — only the parameter sub-
            // ranges change.  The angular formula already handles
            // arbitrary sub-ranges via `Pcurve.t_start, t_end`.
            let left = Pcurve {
                geometry: PcurveGeometry::CylinderArcViaRationalQuad {
                    theta_a: *theta_a,
                    theta_b: *theta_b,
                    v_const: *v_const,
                },
                t_start: p_start,
                t_end: split_param,
            };
            let right = Pcurve {
                geometry: PcurveGeometry::CylinderArcViaRationalQuad {
                    theta_a: *theta_a,
                    theta_b: *theta_b,
                    v_const: *v_const,
                },
                t_start: split_param,
                t_end: p_end,
            };
            (Some(left), Some(right))
        }
    }
}

/// Reverse a NURBS curve (swap direction of traversal).
///
/// The reversed curve evaluates the same geometry but in the opposite
/// direction: `reversed.point_at(t_start)` equals
/// `original.point_at(t_end)`.
pub(crate) fn reverse_nurbs_curve(curve: &NurbsCurve) -> NurbsCurve {
    use crate::geometry::curve::CurveEval;
    let (t_start, t_end) = curve.parameter_range();
    let mut reversed_knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_start + t_end - k)
        .collect();
    reversed_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mut reversed_cp = curve.control_points.clone();
    reversed_cp.reverse();

    let mut reversed_weights = curve.weights.clone();
    reversed_weights.reverse();

    NurbsCurve {
        degree: curve.degree,
        knots: reversed_knots,
        control_points: reversed_cp,
        weights: reversed_weights,
    }
}
