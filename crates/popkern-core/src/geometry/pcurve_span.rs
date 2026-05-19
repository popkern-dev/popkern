//! Pcurve span orientation and sub-range helpers.

use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::pcurve::{Pcurve, PcurveGeometry};

impl Pcurve {
    /// Return the same UV locus traversed in the opposite direction over this
    /// pcurve's parameter range.
    pub fn reversed(&self) -> Self {
        let geometry = match &self.geometry {
            PcurveGeometry::Linear { start, end } => PcurveGeometry::Linear {
                start: *end,
                end: *start,
            },
            PcurveGeometry::Elliptical {
                center,
                u_radius,
                v_radius,
                rotation,
                start_angle,
                end_angle,
            } => PcurveGeometry::Elliptical {
                center: *center,
                u_radius: *u_radius,
                v_radius: *v_radius,
                rotation: *rotation,
                start_angle: *end_angle,
                end_angle: *start_angle,
            },
            PcurveGeometry::CylinderArcViaRationalQuad {
                theta_a,
                theta_b,
                v_const,
            } => PcurveGeometry::CylinderArcViaRationalQuad {
                theta_a: *theta_b,
                theta_b: *theta_a,
                v_const: *v_const,
            },
            PcurveGeometry::Nurbs(curve) => PcurveGeometry::Nurbs(reverse_nurbs_curve(curve)),
            PcurveGeometry::Sinusoidal { .. } => {
                PcurveGeometry::Nurbs(reverse_nurbs_curve(&self.to_nurbs()))
            }
        };
        Self {
            geometry,
            t_start: self.t_start,
            t_end: self.t_end,
        }
    }

    /// Restrict this pcurve to a sub-range while preserving its UV locus.
    pub fn slice_between(&self, t_start: f64, t_end: f64) -> Option<Self> {
        if !t_start.is_finite() || !t_end.is_finite() || (t_end - t_start).abs() <= 1e-14 {
            return None;
        }
        let (lo, hi, reversed) = if t_start <= t_end {
            (t_start, t_end, false)
        } else {
            (t_end, t_start, true)
        };
        let geometry = match &self.geometry {
            PcurveGeometry::Linear { .. } => {
                let start = self.point_at(lo);
                let end = self.point_at(hi);
                PcurveGeometry::Linear {
                    start: (start.x, start.y),
                    end: (end.x, end.y),
                }
            }
            PcurveGeometry::Elliptical {
                center,
                u_radius,
                v_radius,
                rotation,
                start_angle,
                end_angle,
            } => {
                let s0 = normalised_parameter(self, lo);
                let s1 = normalised_parameter(self, hi);
                PcurveGeometry::Elliptical {
                    center: *center,
                    u_radius: *u_radius,
                    v_radius: *v_radius,
                    rotation: *rotation,
                    start_angle: start_angle + s0 * (end_angle - start_angle),
                    end_angle: start_angle + s1 * (end_angle - start_angle),
                }
            }
            PcurveGeometry::Sinusoidal {
                u_offset,
                v_const,
                v_cos,
                v_sin,
            } => PcurveGeometry::Sinusoidal {
                u_offset: *u_offset,
                v_const: *v_const,
                v_cos: *v_cos,
                v_sin: *v_sin,
            },
            PcurveGeometry::CylinderArcViaRationalQuad { v_const, .. } => {
                let start = self.point_at(lo);
                let end = self.point_at(hi);
                PcurveGeometry::CylinderArcViaRationalQuad {
                    theta_a: start.x,
                    theta_b: end.x,
                    v_const: *v_const,
                }
            }
            PcurveGeometry::Nurbs(curve) => {
                PcurveGeometry::Nurbs(slice_nurbs_curve(curve, lo, hi)?)
            }
        };
        let sliced = Self {
            geometry,
            t_start: lo,
            t_end: hi,
        };
        Some(if reversed { sliced.reversed() } else { sliced })
    }
}

fn normalised_parameter(pcurve: &Pcurve, t: f64) -> f64 {
    let span = pcurve.t_end - pcurve.t_start;
    if span.abs() > 0.0 {
        (t - pcurve.t_start) / span
    } else {
        0.0
    }
}

fn reverse_nurbs_curve(curve: &NurbsCurve) -> NurbsCurve {
    use crate::geometry::curve::CurveEval;
    let (t_min, t_max) = curve.parameter_range();
    let mut knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_min + t_max - k)
        .collect();
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut control_points = curve.control_points.clone();
    control_points.reverse();
    let mut weights = curve.weights.clone();
    weights.reverse();
    NurbsCurve {
        degree: curve.degree,
        knots,
        control_points,
        weights,
    }
}

fn slice_nurbs_curve(curve: &NurbsCurve, lo: f64, hi: f64) -> Option<NurbsCurve> {
    use crate::geometry::curve::CurveEval;
    let (full_min, full_max) = curve.parameter_range();
    let eps = 1.0e-12 * (full_max - full_min).abs().max(1.0);
    if lo < full_min - eps || hi > full_max + eps || hi <= lo + eps {
        return None;
    }
    let lo = lo.max(full_min);
    let hi = hi.min(full_max);
    let at_start = (lo - full_min).abs() <= eps;
    let at_end = (hi - full_max).abs() <= eps;
    let mut current = curve.clone();
    if !at_start {
        current = current.split_at(lo).ok()?.1;
    }
    if !at_end {
        let (cur_min, _) = current.parameter_range();
        let cut = if at_start { hi } else { cur_min + (hi - lo) };
        current = current.split_at(cut).ok()?.0;
    }
    Some(current)
}
