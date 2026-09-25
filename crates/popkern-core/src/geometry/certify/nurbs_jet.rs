//! Rigorous jet enclosures of rational tensor-product B-spline
//! surfaces over parameter boxes.
//!
//! Per knot span the homogeneous surface is converted to a Bernstein
//! net (interval blossoms), differentiated as Bernstein nets, and
//! restricted to the query box by interval de Casteljau; the hulls
//! enclose every homogeneous partial over the box. Euclidean partials
//! follow by the Leibniz quotient rule in interval arithmetic.

use super::CertifyDecline;
use super::bernstein::{self, Coef};
use super::interval::Iv;
use super::interval_vec::IvVec3;
use super::surface_jet::{JET_ORDERS, SurfaceJet, jet_slot};
use crate::geometry::nurbs_surface::NurbsSurface;
use std::collections::HashMap;
use std::ops::Div;

type Net = Vec<Vec<Coef<4>>>;

/// Most knot-span pieces one straddling box may be cut into.
const MAX_KNOT_PIECES: usize = 64;

/// Bernstein nets of every homogeneous partial up to order 6 on one
/// knot span.
#[derive(Debug, Clone)]
struct SpanNets {
    u: (f64, f64),
    v: (f64, f64),
    /// Indexed like [`JET_ORDERS`].
    nets: Vec<Net>,
}

/// A NURBS surface prepared for jet enclosure. Span nets are built on
/// first use and cached.
#[derive(Debug, Clone)]
pub struct NurbsJetSource {
    degree: (usize, usize),
    knots_u: Vec<f64>,
    knots_v: Vec<f64>,
    num: (usize, usize),
    /// Homogeneous control points `(w·x, w·y, w·z, w)`, index
    /// `i_u · num_v + i_v`.
    homogeneous: Vec<Coef<4>>,
    cache: HashMap<(usize, usize), SpanNets>,
}

impl NurbsJetSource {
    /// Prepare `surface`. Weights must be positive, except entries
    /// listed in `zero_weight` (ε-regularized Gregory corners whose
    /// exact rational form has weight and numerator zero), which are
    /// set to an exact zero coefficient.
    pub fn new(surface: &NurbsSurface, zero_weight: &[usize]) -> Result<Self, CertifyDecline> {
        let expected = surface.num_u * surface.num_v;
        if surface.control_points.len() != expected
            || surface.weights.len() != expected
            || surface.knots_u.len() != surface.num_u + surface.degree_u + 1
            || surface.knots_v.len() != surface.num_v + surface.degree_v + 1
        {
            return Err(CertifyDecline::MalformedNurbs);
        }
        let mut homogeneous = Vec::with_capacity(expected);
        for (k, (p, &w)) in surface
            .control_points
            .iter()
            .zip(surface.weights.iter())
            .enumerate()
        {
            if zero_weight.contains(&k) {
                homogeneous.push([Iv::point(0.0); 4]);
                continue;
            }
            if w <= 0.0 || !w.is_finite() {
                return Err(CertifyDecline::NonPositiveWeight);
            }
            let w_iv = Iv::point(w);
            homogeneous.push([
                Iv::point(p.x) * w_iv,
                Iv::point(p.y) * w_iv,
                Iv::point(p.z) * w_iv,
                w_iv,
            ]);
        }
        Ok(Self {
            degree: (surface.degree_u, surface.degree_v),
            knots_u: surface.knots_u.clone(),
            knots_v: surface.knots_v.clone(),
            num: (surface.num_u, surface.num_v),
            homogeneous,
            cache: HashMap::new(),
        })
    }

    /// Distinct interior knot values in u and v: seam intervals are
    /// split wherever a pcurve's image crosses one.
    pub fn breaks(&self) -> (Vec<f64>, Vec<f64>) {
        (
            interior_breaks(&self.knots_u, self.degree.0, self.num.0),
            interior_breaks(&self.knots_v, self.degree.1, self.num.1),
        )
    }

    /// Enclosure of every partial up to order 6 over a box that may
    /// straddle knots: the box is cut at every interior knot it
    /// contains and the per-span enclosures are hulled. Returns
    /// whether a cut happened — across a knot the partials above the
    /// continuity order are two-valued, so a consumer must not apply
    /// a mean-value form over a straddled box.
    pub fn jet_piecewise(&mut self, u: Iv, v: Iv) -> Result<(SurfaceJet, bool), CertifyDecline> {
        let (bu, bv) = self.breaks();
        let pieces = |x: Iv, breaks: &[f64]| -> Vec<Iv> {
            let mut cuts: Vec<f64> = vec![x.lo];
            cuts.extend(breaks.iter().copied().filter(|&k| k > x.lo && k < x.hi));
            cuts.push(x.hi);
            cuts.windows(2).map(|w| Iv::new(w[0], w[1])).collect()
        };
        let (pu, pv) = (pieces(u, &bu), pieces(v, &bv));
        if pu.len() * pv.len() > MAX_KNOT_PIECES {
            return Err(CertifyDecline::SpanStraddle);
        }
        let mut out: Option<SurfaceJet> = None;
        for &x in &pu {
            for &y in &pv {
                let jet = self.jet(x, y)?;
                out = Some(match out {
                    Some(acc) => acc.hull(&jet),
                    None => jet,
                });
            }
        }
        let straddled = pu.len() * pv.len() > 1;
        out.map(|jet| (jet, straddled))
            .ok_or(CertifyDecline::MalformedNurbs)
    }

    /// Enclosure of every partial up to order 6 over the box.
    pub fn jet(&mut self, u: Iv, v: Iv) -> Result<SurfaceJet, CertifyDecline> {
        let su = span_of(&self.knots_u, self.degree.0, self.num.0, u)?;
        let sv = span_of(&self.knots_v, self.degree.1, self.num.1, v)?;
        self.jet_in_span(u, v, (su, sv))
    }

    /// Knot spans `(u, v)` whose closed parameter intervals meet the
    /// box (the first/last span also covers extrapolation beyond the
    /// domain).
    pub fn spans_touching(&self, u: Iv, v: Iv) -> Result<Vec<(usize, usize)>, CertifyDecline> {
        if !u.is_finite() || !v.is_finite() {
            return Err(CertifyDecline::NonFinite);
        }
        let touching = |knots: &[f64], degree: usize, num: usize, x: Iv| -> Vec<usize> {
            let spans: Vec<usize> = (degree..num).filter(|&k| knots[k] < knots[k + 1]).collect();
            let (first, last) = (spans.first().copied(), spans.last().copied());
            spans
                .iter()
                .copied()
                .filter(|&k| {
                    let lo_ok = x.hi >= knots[k] || Some(k) == first;
                    let hi_ok = x.lo <= knots[k + 1] || Some(k) == last;
                    lo_ok && hi_ok
                })
                .collect()
        };
        let su = touching(&self.knots_u, self.degree.0, self.num.0, u);
        let sv = touching(&self.knots_v, self.degree.1, self.num.1, v);
        Ok(su
            .iter()
            .flat_map(|&a| sv.iter().map(move |&b| (a, b)))
            .collect())
    }

    /// Enclosure over the box of the polynomial of span `(su, sv)`,
    /// extended beyond that span where the box leaves it: a smooth
    /// function that equals the surface wherever the span applies.
    pub fn jet_in_span(
        &mut self,
        u: Iv,
        v: Iv,
        (su, sv): (usize, usize),
    ) -> Result<SurfaceJet, CertifyDecline> {
        if !self.cache.contains_key(&(su, sv)) {
            let nets = self.span_nets(su, sv);
            self.cache.insert((su, sv), nets);
        }
        let span = self
            .cache
            .get(&(su, sv))
            .ok_or(CertifyDecline::MalformedNurbs)?;
        let (au, bu) = (
            bernstein::local_param(u.lo, span.u.0, span.u.1),
            bernstein::local_param(u.hi, span.u.0, span.u.1),
        );
        let (av, bv) = (
            bernstein::local_param(v.lo, span.v.0, span.v.1),
            bernstein::local_param(v.hi, span.v.0, span.v.1),
        );
        let homogeneous: Vec<Coef<4>> = span
            .nets
            .iter()
            .map(|net| {
                if net.is_empty() || net[0].is_empty() {
                    return [Iv::point(0.0); 4];
                }
                let across_u: Vec<Vec<Coef<4>>> = transpose(net)
                    .iter()
                    .map(|column| bernstein::restrict(column, au, bu))
                    .collect();
                let rows = transpose(&across_u);
                let restricted: Vec<Coef<4>> = rows
                    .iter()
                    .flat_map(|row| bernstein::restrict(row, av, bv))
                    .collect();
                bernstein::hull(&restricted)
            })
            .collect();
        rational_quotient(&homogeneous)
    }

    fn span_nets(&self, su: usize, sv: usize) -> SpanNets {
        let (p, q) = self.degree;
        let num_v = self.num.1;
        // v-direction first: each u-row of the window → Bernstein in v.
        let rows: Net = (su - p..=su)
            .map(|i| {
                let ctrl: Vec<Coef<4>> = (sv - q..=sv)
                    .map(|j| self.homogeneous[i * num_v + j])
                    .collect();
                bernstein::bspline_span_to_bernstein(q, &self.knots_v, sv, &ctrl)
            })
            .collect();
        let columns = transpose(&rows);
        let base: Net = transpose(
            &columns
                .iter()
                .map(|col| bernstein::bspline_span_to_bernstein(p, &self.knots_u, su, col))
                .collect::<Vec<_>>(),
        );
        let u = (self.knots_u[su], self.knots_u[su + 1]);
        let v = (self.knots_v[sv], self.knots_v[sv + 1]);
        let (wu, wv) = (
            Iv::point(u.1) - Iv::point(u.0),
            Iv::point(v.1) - Iv::point(v.0),
        );
        let nets = JET_ORDERS
            .iter()
            .map(|&(a, b)| {
                let mut net = base.clone();
                for _ in 0..a {
                    net = transpose(
                        &transpose(&net)
                            .iter()
                            .map(|col| bernstein::derivative_net(col, wu))
                            .collect::<Vec<_>>(),
                    );
                }
                for _ in 0..b {
                    net = net
                        .iter()
                        .map(|row| bernstein::derivative_net(row, wv))
                        .collect();
                }
                net
            })
            .collect();
        SpanNets { u, v, nets }
    }
}

/// `net[i][j]` → `out[j][i]`.
fn transpose(net: &[Vec<Coef<4>>]) -> Net {
    let cols = net.first().map_or(0, Vec::len);
    (0..cols)
        .map(|j| net.iter().map(|row| row[j]).collect())
        .collect()
}

fn interior_breaks(knots: &[f64], degree: usize, num: usize) -> Vec<f64> {
    let (lo, hi) = (knots[degree], knots[num]);
    let mut out: Vec<f64> = knots
        .iter()
        .copied()
        .filter(|&k| k > lo && k < hi)
        .collect();
    out.dedup();
    out
}

/// The non-empty knot span containing `[x.lo, x.hi]`. Boxes beyond
/// the domain extend the first or last span's polynomial (as the
/// kernel's evaluator does); a box straddling an interior knot is
/// declined, since the caller splits there.
fn span_of(knots: &[f64], degree: usize, num: usize, x: Iv) -> Result<usize, CertifyDecline> {
    if !x.is_finite() {
        return Err(CertifyDecline::NonFinite);
    }
    let spans: Vec<usize> = (degree..num).filter(|&k| knots[k] < knots[k + 1]).collect();
    let Some(&last) = spans.last() else {
        return Err(CertifyDecline::MalformedNurbs);
    };
    let locate = |t: f64| {
        spans
            .iter()
            .copied()
            .find(|&k| t < knots[k + 1])
            .unwrap_or(last)
    };
    let (a, b) = (locate(x.lo), locate(x.hi));
    if a == b {
        return Ok(a);
    }
    // x.hi sitting exactly on the knot that closes x.lo's span belongs
    // to x.lo's span (closed span).
    if x.hi == knots[a + 1] {
        return Ok(a);
    }
    Err(CertifyDecline::SpanStraddle)
}

/// Euclidean partials from homogeneous ones by the Leibniz quotient:
/// `S_kl = (A_kl − Σ_{(i,j)≠0} C(k,i) C(l,j) w_ij S_{k−i,l−j}) / w`.
pub(super) fn rational_quotient(homogeneous: &[Coef<4>]) -> Result<SurfaceJet, CertifyDecline> {
    let w = homogeneous[jet_slot(0, 0)][3];
    if w.contains_zero() || !w.is_finite() {
        return Err(CertifyDecline::DenominatorContainsZero);
    }
    let binom = |n: usize, k: usize| -> f64 {
        (0..k).fold(1.0, |acc, m| acc * (n - m) as f64 / (m + 1) as f64)
    };
    let mut jet = SurfaceJet::entire();
    for &(k, l) in JET_ORDERS.iter() {
        let a = homogeneous[jet_slot(k, l)];
        let mut acc = IvVec3::new(a[0], a[1], a[2]);
        for i in 0..=k {
            for j in 0..=l {
                if i == 0 && j == 0 {
                    continue;
                }
                let c = binom(k, i) * binom(l, j);
                let wij = homogeneous[jet_slot(i, j)][3];
                acc = acc - jet.partial(k - i, l - j).scale(wij * Iv::point(c));
            }
        }
        jet.set(k, l, IvVec3::new(acc.x.div(w), acc.y.div(w), acc.z.div(w)));
    }
    Ok(jet)
}
