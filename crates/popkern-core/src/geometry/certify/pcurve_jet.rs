//! Pcurve jet enclosures: `(u, v)(t)` and its first two derivatives
//! over a parameter interval, per pcurve kind, replicating the
//! kernel's evaluator (`Pcurve::point_at`).

use super::CertifyDecline;
use super::bernstein::{self, Coef};
use super::interval::Iv;
use super::interval_vec::IvVec2;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::pcurve::{Pcurve, PcurveGeometry};
use std::ops::Div;

/// Enclosures of a pcurve's value and first two derivatives with
/// respect to its own parameter over an interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PcurveJet {
    /// `(u, v)`.
    pub value: IvVec2,
    /// `d(u, v)/dt`.
    pub d1: IvVec2,
    /// `d²(u, v)/dt²`.
    pub d2: IvVec2,
}

impl PcurveJet {
    fn is_finite(&self) -> bool {
        self.value.is_finite() && self.d1.is_finite() && self.d2.is_finite()
    }
}

/// A pcurve prepared for rigorous jet enclosure.
#[derive(Debug, Clone)]
pub struct PcurveJetSource {
    kind: Kind,
    t_range: (f64, f64),
}

#[derive(Debug, Clone)]
enum Kind {
    Linear {
        start: (f64, f64),
        end: (f64, f64),
    },
    Elliptical {
        center: (f64, f64),
        radii: (f64, f64),
        rotation: f64,
        angles: (f64, f64),
    },
    Sinusoidal {
        u_offset: f64,
        v_const: f64,
        v_cos: f64,
        v_sin: f64,
    },
    ArcViaRationalQuad {
        theta: (f64, f64),
        v_const: f64,
    },
    Nurbs {
        degree: usize,
        knots: Vec<f64>,
        homogeneous: Vec<Coef<3>>,
    },
}

fn pt(x: f64) -> Iv {
    Iv::point(x)
}

impl PcurveJetSource {
    /// Prepare `pcurve`, or say why it cannot be certified.
    pub fn new(pcurve: &Pcurve) -> Result<Self, CertifyDecline> {
        let t_range = (pcurve.t_start, pcurve.t_end);
        let normalised = matches!(
            pcurve.geometry,
            PcurveGeometry::Linear { .. }
                | PcurveGeometry::Elliptical { .. }
                | PcurveGeometry::CylinderArcViaRationalQuad { .. }
        );
        if normalised && t_range.0 == t_range.1 {
            return Err(CertifyDecline::DegeneratePcurve);
        }
        let kind = match &pcurve.geometry {
            PcurveGeometry::Linear { start, end } => Kind::Linear {
                start: *start,
                end: *end,
            },
            PcurveGeometry::Elliptical {
                center,
                u_radius,
                v_radius,
                rotation,
                start_angle,
                end_angle,
            } => Kind::Elliptical {
                center: *center,
                radii: (*u_radius, *v_radius),
                rotation: *rotation,
                angles: (*start_angle, *end_angle),
            },
            PcurveGeometry::Sinusoidal {
                u_offset,
                v_const,
                v_cos,
                v_sin,
            } => Kind::Sinusoidal {
                u_offset: *u_offset,
                v_const: *v_const,
                v_cos: *v_cos,
                v_sin: *v_sin,
            },
            PcurveGeometry::CylinderArcViaRationalQuad {
                theta_a,
                theta_b,
                v_const,
            } => {
                // atan2 reduces to atan only while cos α > 0.
                if (theta_b - theta_a).abs() >= std::f64::consts::PI {
                    return Err(CertifyDecline::UnsupportedPcurve(
                        "rational-quad arc with |θb − θa| ≥ π",
                    ));
                }
                Kind::ArcViaRationalQuad {
                    theta: (*theta_a, *theta_b),
                    v_const: *v_const,
                }
            }
            PcurveGeometry::Nurbs(curve) => nurbs_kind(curve)?,
        };
        Ok(Self { kind, t_range })
    }

    /// Interior knots of a NURBS pcurve (empty for the other kinds).
    pub fn breaks(&self) -> Vec<f64> {
        match &self.kind {
            Kind::Nurbs { degree, knots, .. } => {
                let (lo, hi) = (knots[*degree], knots[knots.len() - degree - 1]);
                let mut out: Vec<f64> = knots
                    .iter()
                    .copied()
                    .filter(|&k| k > lo && k < hi)
                    .collect();
                out.dedup();
                out
            }
            _ => Vec::new(),
        }
    }

    /// `s = (t − t_start)/(t_end − t_start)` and `ds/dt`.
    fn normalised(&self, t: Iv) -> (Iv, Iv) {
        let span = pt(self.t_range.1) - pt(self.t_range.0);
        ((t - pt(self.t_range.0)).div(span), pt(1.0).div(span))
    }

    /// Enclosure of the jet over the parameter interval `t`.
    pub fn jet(&self, t: Iv) -> Result<PcurveJet, CertifyDecline> {
        let zero = IvVec2::new(pt(0.0), pt(0.0));
        let jet = match &self.kind {
            Kind::Linear { start, end } => {
                let (s, ds) = self.normalised(t);
                let delta = IvVec2::new(pt(end.0) - pt(start.0), pt(end.1) - pt(start.1));
                PcurveJet {
                    value: IvVec2::new(pt(start.0), pt(start.1)) + delta.scale(s),
                    d1: delta.scale(ds),
                    d2: zero,
                }
            }
            Kind::Elliptical {
                center,
                radii,
                rotation,
                angles,
            } => {
                let (s, ds) = self.normalised(t);
                let sweep = pt(angles.1) - pt(angles.0);
                let alpha = pt(angles.0) + s * sweep;
                let da = sweep * ds;
                let (ca, sa) = (alpha.cos(), alpha.sin());
                let (cr, sr) = (pt(*rotation).cos(), pt(*rotation).sin());
                let (ru, rv) = (pt(radii.0), pt(radii.1));
                // Local (ru cos α, rv sin α) and its α-derivatives.
                let local = [
                    (ru * ca, rv * sa),
                    (-(ru * sa), rv * ca),
                    (-(ru * ca), -(rv * sa)),
                ];
                let rotate = |(lu, lv): (Iv, Iv)| IvVec2::new(lu * cr - lv * sr, lu * sr + lv * cr);
                PcurveJet {
                    value: IvVec2::new(pt(center.0), pt(center.1)) + rotate(local[0]),
                    d1: rotate(local[1]).scale(da),
                    d2: rotate(local[2]).scale(da.sqr()),
                }
            }
            Kind::Sinusoidal {
                u_offset,
                v_const,
                v_cos,
                v_sin,
            } => {
                let (c, s) = (t.cos(), t.sin());
                let (kc, ks) = (pt(*v_cos), pt(*v_sin));
                PcurveJet {
                    value: IvVec2::new(pt(*u_offset) + t, pt(*v_const) + kc * c + ks * s),
                    d1: IvVec2::new(pt(1.0), ks * c - kc * s),
                    d2: IvVec2::new(pt(0.0), -(kc * c) - ks * s),
                }
            }
            Kind::ArcViaRationalQuad { theta, v_const } => {
                let (s, ds) = self.normalised(t);
                arc_via_rational_quad(*theta, *v_const, s, ds)
            }
            Kind::Nurbs {
                degree,
                knots,
                homogeneous,
            } => nurbs_jet(*degree, knots, homogeneous, t)?,
        };
        if jet.is_finite() {
            Ok(jet)
        } else {
            Err(CertifyDecline::NonFinite)
        }
    }
}

/// `θ(s) = θm + atan(y/x)` with `y = sin α (2s − 1)`,
/// `x = cos α (1 − 2s(1−s)) + 2s(1−s)` (positive for |α| < π/2, so
/// atan2 is atan); derivatives by the quotient rule on `N/D`,
/// `N = x y′ − y x′`, `D = x² + y²`.
fn arc_via_rational_quad(theta: (f64, f64), v_const: f64, s: Iv, ds: Iv) -> PcurveJet {
    let alpha = (pt(theta.1) - pt(theta.0)).scale(0.5);
    let theta_m = (pt(theta.0) + pt(theta.1)).scale(0.5);
    let (sa, ca) = (alpha.sin(), alpha.cos());
    let d = (s * (pt(1.0) - s)).scale(2.0);
    let dd = (pt(1.0) - s.scale(2.0)).scale(2.0);
    let ddd = pt(-4.0);
    let one_minus_c = pt(1.0) - ca;
    let y = sa * (s.scale(2.0) - pt(1.0));
    let y1 = sa.scale(2.0);
    let x = ca + one_minus_c * d;
    let x1 = one_minus_c * dd;
    let x2 = one_minus_c * ddd;
    let n = x * y1 - y * x1;
    let n1 = -(y * x2);
    let den = x.sqr() + y.sqr();
    let den1 = (x * x1 + y * y1).scale(2.0);
    let dtheta = n.div(den);
    let d2theta = (n1 * den - n * den1).div(den.sqr());
    let value = theta_m + y.div(x).atan();
    PcurveJet {
        value: IvVec2::new(value, pt(v_const)),
        d1: IvVec2::new(dtheta * ds, pt(0.0)),
        d2: IvVec2::new(d2theta * ds.sqr(), pt(0.0)),
    }
}

fn nurbs_kind(curve: &NurbsCurve) -> Result<Kind, CertifyDecline> {
    let n = curve.control_points.len();
    if curve.weights.len() != n || curve.knots.len() != n + curve.degree + 1 || n == 0 {
        return Err(CertifyDecline::MalformedNurbs);
    }
    let mut homogeneous = Vec::with_capacity(n);
    for (p, &w) in curve.control_points.iter().zip(curve.weights.iter()) {
        if w <= 0.0 || !w.is_finite() {
            return Err(CertifyDecline::NonPositiveWeight);
        }
        let wi = pt(w);
        homogeneous.push([pt(p.x) * wi, pt(p.y) * wi, wi]);
    }
    Ok(Kind::Nurbs {
        degree: curve.degree,
        knots: curve.knots.clone(),
        homogeneous,
    })
}

fn nurbs_jet(
    degree: usize,
    knots: &[f64],
    homogeneous: &[Coef<3>],
    t: Iv,
) -> Result<PcurveJet, CertifyDecline> {
    let n = homogeneous.len();
    let spans: Vec<usize> = (degree..n).filter(|&k| knots[k] < knots[k + 1]).collect();
    let last = *spans.last().ok_or(CertifyDecline::MalformedNurbs)?;
    let locate = |x: f64| {
        spans
            .iter()
            .copied()
            .find(|&k| x < knots[k + 1])
            .unwrap_or(last)
    };
    let (a, b) = (locate(t.lo), locate(t.hi));
    let span = if a == b || t.hi == knots[a + 1] {
        a
    } else {
        return Err(CertifyDecline::SpanStraddle);
    };
    let ctrl = &homogeneous[span - degree..=span];
    let base = bernstein::bspline_span_to_bernstein(degree, knots, span, ctrl);
    let width = pt(knots[span + 1]) - pt(knots[span]);
    let first = bernstein::derivative_net(&base, width);
    let second = bernstein::derivative_net(&first, width);
    let (lo, hi) = (
        bernstein::local_param(t.lo, knots[span], knots[span + 1]),
        bernstein::local_param(t.hi, knots[span], knots[span + 1]),
    );
    let enclose = |net: &[Coef<3>]| -> Coef<3> {
        if net.is_empty() {
            [pt(0.0); 3]
        } else {
            bernstein::hull(&bernstein::restrict(net, lo, hi))
        }
    };
    let (h0, h1, h2) = (enclose(&base), enclose(&first), enclose(&second));
    let w = h0[2];
    if w.contains_zero() {
        return Err(CertifyDecline::DenominatorContainsZero);
    }
    let c0 = [h0[0].div(w), h0[1].div(w)];
    let c1 = [0, 1].map(|k| (h1[k] - h1[2] * c0[k]).div(w));
    let c2 = [0, 1].map(|k| (h2[k] - (h1[2] * c1[k]).scale(2.0) - h2[2] * c0[k]).div(w));
    Ok(PcurveJet {
        value: IvVec2::new(c0[0], c0[1]),
        d1: IvVec2::new(c1[0], c1[1]),
        d2: IvVec2::new(c2[0], c2[1]),
    })
}
