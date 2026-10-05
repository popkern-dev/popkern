//! Pcurve jet enclosures: `(u, v)(t)` and its derivatives up to
//! order [`PCURVE_JET_TERMS`]` − 1` over a parameter interval, per
//! pcurve kind, replicating the kernel's evaluator
//! (`Pcurve::point_at`).

use super::bernstein::{self, Coef};
use super::interval::Iv;
use super::interval_vec::IvVec2;
use super::pcurve_affine::AffinePolyline;
use super::taylor::{SEAM_TAYLOR_ORDER, Taylor};
use super::{CertifyDecline, UnsupportedGeometry};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::pcurve::{Pcurve, PcurveGeometry};
use std::ops::Div;

/// Derivative count a pcurve jet carries: the seam position series
/// runs one order past the Taylor order (its derivative is the seam
/// tangent series).
pub const PCURVE_JET_TERMS: usize = SEAM_TAYLOR_ORDER + 2;

/// Enclosures of a pcurve's value and derivatives `d[k] = dᵏ(u, v)/dtᵏ`
/// with respect to its own parameter over an interval.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PcurveJet {
    /// Derivatives, `d[0]` the value.
    pub d: [IvVec2; PCURVE_JET_TERMS],
}

impl PcurveJet {
    fn is_finite(&self) -> bool {
        self.d.iter().all(|x| x.is_finite())
    }

    /// Component-wise hull.
    pub fn hull(&self, o: &PcurveJet) -> PcurveJet {
        PcurveJet {
            d: std::array::from_fn(|k| self.d[k].hull(o.d[k])),
        }
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
    /// A degree-1 pcurve that is one affine map up to roundoff.
    Affine(AffinePolyline),
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
                        UnsupportedGeometry::RationalQuadArcHalfTurn,
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

    /// Enclosure over any parameter interval, hulled across the NURBS
    /// knots it contains; the flag says whether it straddled one.
    pub fn jet_piecewise(&self, t: Iv) -> Result<(PcurveJet, bool), CertifyDecline> {
        let mut cuts = vec![t.lo];
        cuts.extend(self.breaks().into_iter().filter(|&k| k > t.lo && k < t.hi));
        cuts.push(t.hi);
        let mut out: Option<PcurveJet> = None;
        for w in cuts.windows(2) {
            let jet = self.jet(Iv::new(w[0], w[1]))?;
            out = Some(match out {
                Some(acc) => acc.hull(&jet),
                None => jet,
            });
        }
        out.map(|jet| (jet, cuts.len() > 2))
            .ok_or(CertifyDecline::MalformedNurbs)
    }

    /// NURBS knot spans meeting `t` (`None`: a single smooth branch,
    /// for the analytic kinds).
    pub fn spans_touching(&self, t: Iv) -> Vec<Option<usize>> {
        match &self.kind {
            Kind::Nurbs {
                degree,
                knots,
                homogeneous,
            } => nurbs_spans_touching(*degree, knots, homogeneous.len(), t)
                .into_iter()
                .map(Some)
                .collect(),
            _ => vec![None],
        }
    }

    /// Jet over `t` of one NURBS span's (extended) polynomial, or of
    /// the analytic pcurve.
    pub fn jet_in(&self, t: Iv, span: Option<usize>) -> Result<PcurveJet, CertifyDecline> {
        match (&self.kind, span) {
            (
                Kind::Nurbs {
                    degree,
                    knots,
                    homogeneous,
                },
                Some(span),
            ) => {
                let jet = PcurveJet {
                    d: nurbs_jet_in_span(*degree, knots, homogeneous, t, span)?,
                };
                if jet.is_finite() {
                    Ok(jet)
                } else {
                    Err(CertifyDecline::NonFinite)
                }
            }
            _ => self.jet(t),
        }
    }

    /// Enclosure of the jet over the parameter interval `t`.
    pub fn jet(&self, t: Iv) -> Result<PcurveJet, CertifyDecline> {
        let zero = IvVec2::new(pt(0.0), pt(0.0));
        let mut d = [zero; PCURVE_JET_TERMS];
        match &self.kind {
            Kind::Linear { start, end } => {
                let (s, ds) = self.normalised(t);
                let delta = IvVec2::new(pt(end.0) - pt(start.0), pt(end.1) - pt(start.1));
                d[0] = IvVec2::new(pt(start.0), pt(start.1)) + delta.scale(s);
                d[1] = delta.scale(ds);
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
                let (ca, sa) = (Iv::cos(alpha), Iv::sin(alpha));
                let (cr, sr) = (Iv::cos(pt(*rotation)), Iv::sin(pt(*rotation)));
                let (ru, rv) = (pt(radii.0), pt(radii.1));
                let mut rate = pt(1.0);
                for (k, slot) in d.iter_mut().enumerate() {
                    // dᵏ/dαᵏ of (ru cos α, rv sin α), rotated.
                    let (lu, lv) = (ru * cos_deriv(k, ca, sa), rv * sin_deriv(k, ca, sa));
                    let rotated = IvVec2::new(lu * cr - lv * sr, lu * sr + lv * cr);
                    *slot = if k == 0 {
                        IvVec2::new(pt(center.0), pt(center.1)) + rotated
                    } else {
                        rotated.scale(rate)
                    };
                    rate = rate * da;
                }
            }
            Kind::Sinusoidal {
                u_offset,
                v_const,
                v_cos,
                v_sin,
            } => {
                let (c, s) = (Iv::cos(t), Iv::sin(t));
                let (kc, ks) = (pt(*v_cos), pt(*v_sin));
                for (k, slot) in d.iter_mut().enumerate() {
                    let v = kc * cos_deriv(k, c, s) + ks * sin_deriv(k, c, s);
                    let u = match k {
                        0 => pt(*u_offset) + t,
                        1 => pt(1.0),
                        _ => pt(0.0),
                    };
                    *slot = if k == 0 {
                        IvVec2::new(u, pt(*v_const) + v)
                    } else {
                        IvVec2::new(u, v)
                    };
                }
            }
            Kind::ArcViaRationalQuad { theta, v_const } => {
                let (s, ds) = self.normalised(t);
                d = arc_via_rational_quad(*theta, *v_const, s, ds);
            }
            Kind::Nurbs {
                degree,
                knots,
                homogeneous,
            } => d = nurbs_jet(*degree, knots, homogeneous, t)?,
            Kind::Affine(affine) => d = affine.jet(t),
        }
        let jet = PcurveJet { d };
        if jet.is_finite() {
            Ok(jet)
        } else {
            Err(CertifyDecline::NonFinite)
        }
    }
}

/// `dᵏ/dxᵏ cos x` from enclosures of `cos x`, `sin x`.
fn cos_deriv(k: usize, c: Iv, s: Iv) -> Iv {
    match k % 4 {
        0 => c,
        1 => -s,
        2 => -c,
        _ => s,
    }
}

/// `dᵏ/dxᵏ sin x`.
fn sin_deriv(k: usize, c: Iv, s: Iv) -> Iv {
    match k % 4 {
        0 => s,
        1 => c,
        2 => -s,
        _ => -c,
    }
}

/// `θ(s) = θm + atan(y/x)` with `y = sin α (2s − 1)`,
/// `x = cos α (1 − 2s(1−s)) + 2s(1−s)` (positive for |α| < π/2, so
/// atan2 is atan), expanded as a Taylor series about `s` (exact
/// polynomials `x`, `y`, then Taylor division and atan).
fn arc_via_rational_quad(
    theta: (f64, f64),
    v_const: f64,
    s: Iv,
    ds: Iv,
) -> [IvVec2; PCURVE_JET_TERMS] {
    let alpha = (pt(theta.1) - pt(theta.0)).scale(0.5);
    let theta_m = (pt(theta.0) + pt(theta.1)).scale(0.5);
    let (sa, ca) = (Iv::sin(alpha), Iv::cos(alpha));
    let one_minus_c = pt(1.0) - ca;
    // x(s+δ) = x + x′δ + ½x″δ², y(s+δ) = y + y′δ.
    let d = (s * (pt(1.0) - s)).scale(2.0);
    let mut x = [pt(0.0); PCURVE_JET_TERMS];
    x[0] = ca + one_minus_c * d;
    x[1] = one_minus_c * (pt(1.0) - s.scale(2.0)).scale(2.0);
    x[2] = one_minus_c.scale(-2.0);
    let mut y = [pt(0.0); PCURVE_JET_TERMS];
    y[0] = sa * (s.scale(2.0) - pt(1.0));
    y[1] = sa.scale(2.0);
    let angle = Taylor::atan(Taylor::new(y) / Taylor::new(x));
    let mut out = [IvVec2::new(pt(0.0), pt(0.0)); PCURVE_JET_TERMS];
    let (mut factorial, mut rate) = (1.0, pt(1.0));
    for (k, slot) in out.iter_mut().enumerate() {
        if k > 0 {
            factorial *= k as f64;
        }
        let dk = angle.c[k].scale(factorial) * rate;
        *slot = if k == 0 {
            IvVec2::new(theta_m + angle.c[0], pt(v_const))
        } else {
            IvVec2::new(dk, pt(0.0))
        };
        rate = rate * ds;
    }
    out
}

fn nurbs_kind(curve: &NurbsCurve) -> Result<Kind, CertifyDecline> {
    if let Some(affine) = AffinePolyline::of(curve) {
        return Ok(Kind::Affine(affine));
    }
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
) -> Result<[IvVec2; PCURVE_JET_TERMS], CertifyDecline> {
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
    nurbs_jet_in_span(degree, knots, homogeneous, t, span)
}

/// The NURBS pcurve's knot spans whose closed intervals meet `t`.
fn nurbs_spans_touching(degree: usize, knots: &[f64], n: usize, t: Iv) -> Vec<usize> {
    let spans: Vec<usize> = (degree..n).filter(|&k| knots[k] < knots[k + 1]).collect();
    let (first, last) = (spans.first().copied(), spans.last().copied());
    spans
        .iter()
        .copied()
        .filter(|&k| {
            (t.hi >= knots[k] || Some(k) == first) && (t.lo <= knots[k + 1] || Some(k) == last)
        })
        .collect()
}

/// Jet of span `span`'s polynomial over `t`, extended beyond the span
/// where `t` leaves it.
fn nurbs_jet_in_span(
    degree: usize,
    knots: &[f64],
    homogeneous: &[Coef<3>],
    t: Iv,
    span: usize,
) -> Result<[IvVec2; PCURVE_JET_TERMS], CertifyDecline> {
    let ctrl = &homogeneous[span - degree..=span];
    let width = pt(knots[span + 1]) - pt(knots[span]);
    let (lo, hi) = (
        bernstein::local_param(t.lo, knots[span], knots[span + 1]),
        bernstein::local_param(t.hi, knots[span], knots[span + 1]),
    );
    // Homogeneous derivative hulls h[k] = dᵏ(wu, wv, w)/dtᵏ.
    let mut net = bernstein::bspline_span_to_bernstein(degree, knots, span, ctrl);
    let mut h = [[pt(0.0); 3]; PCURVE_JET_TERMS];
    for slot in h.iter_mut() {
        if net.is_empty() {
            break;
        }
        *slot = bernstein::hull(&bernstein::restrict(&net, lo, hi));
        net = bernstein::derivative_net(&net, width);
    }
    let w = h[0][2];
    if w.contains_zero() {
        return Err(CertifyDecline::DenominatorContainsZero);
    }
    // Leibniz: C⁽ᵏ⁾ = (A⁽ᵏ⁾ − Σ_{j≥1} C(k,j) w⁽ʲ⁾ C⁽ᵏ⁻ʲ⁾) / w.
    let binom = |n: usize, k: usize| -> f64 {
        (0..k).fold(1.0, |acc, m| acc * (n - m) as f64 / (m + 1) as f64)
    };
    let mut c = [[pt(0.0); 2]; PCURVE_JET_TERMS];
    for k in 0..PCURVE_JET_TERMS {
        for axis in 0..2 {
            let mut acc = h[k][axis];
            for j in 1..=k {
                acc = acc - (h[j][2] * c[k - j][axis]).scale(binom(k, j));
            }
            c[k][axis] = acc.div(w);
        }
    }
    Ok(c.map(|[u, v]| IvVec2::new(u, v)))
}
