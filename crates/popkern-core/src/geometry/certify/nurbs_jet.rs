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

#[path = "span_box.rs"]
mod span_box;

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
    /// Box enclosure buffers and iso-line caches (see [`span_box`]).
    boxes: span_box::SpanBoxEncloser,
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
            boxes: span_box::SpanBoxEncloser::default(),
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
        let homogeneous = self.boxes.enclose((su, sv), &span.nets, (au, bu), (av, bv));
        rational_quotient(&homogeneous)
    }

    /// An enclosure of the surface's position alone over a box inside
    /// the domain: the hull of the positive-weight entries of each
    /// span's homogeneous net restricted to the box. With non-negative
    /// weights the rational patch is a convex combination of those
    /// entries, and an entry of weight exactly zero has a zero
    /// numerator (an ε-regularized Gregory corner), so it contributes
    /// nothing. Unlike [`NurbsJetSource::jet`] nothing is divided by
    /// the weight function, so a box touching a corner where it
    /// vanishes still encloses (PK-904).
    ///
    /// # Errors
    ///
    /// [`CertifyDecline::DenominatorContainsZero`] when an entry not
    /// known to be zero has a weight enclosure reaching zero;
    /// [`CertifyDecline::SpanStraddle`] for a box leaving the domain.
    pub fn position_hull(&mut self, u: Iv, v: Iv) -> Result<IvVec3, CertifyDecline> {
        let (p, q) = self.degree;
        let domain = |knots: &[f64], degree: usize, num: usize| (knots[degree], knots[num]);
        let (du, dv) = (
            domain(&self.knots_u, p, self.num.0),
            domain(&self.knots_v, q, self.num.1),
        );
        if u.lo < du.0 || u.hi > du.1 || v.lo < dv.0 || v.hi > dv.1 {
            return Err(CertifyDecline::SpanStraddle);
        }
        let mut out: Option<IvVec3> = None;
        for (su, sv) in self.spans_touching(u, v)? {
            if !self.cache.contains_key(&(su, sv)) {
                let nets = self.span_nets(su, sv);
                self.cache.insert((su, sv), nets);
            }
            let span = self
                .cache
                .get(&(su, sv))
                .ok_or(CertifyDecline::MalformedNurbs)?;
            // The box is clipped to the span, so its local parameters lie
            // in [0, 1]: outward rounding past either end would mix a
            // zero-weight corner entry into a weight straddling zero.
            // A box end on a span end is exactly 0 or 1 (the zero
            // tracking below needs them exact).
            let unit = |x: f64, s: (f64, f64)| {
                if x == s.0 {
                    Iv::point(0.0)
                } else if x == s.1 {
                    Iv::point(1.0)
                } else {
                    let t = bernstein::local_param(x, s.0, s.1);
                    Iv::new(t.lo.clamp(0.0, 1.0), t.hi.clamp(0.0, 1.0))
                }
            };
            let local = |x: Iv, s: (f64, f64)| {
                let (lo, hi) = (x.lo.max(s.0), x.hi.min(s.1));
                (unit(lo, s), unit(hi.max(lo), s))
            };
            let ((au, bu), (av, bv)) = (local(u, span.u), local(v, span.v));
            // On a Bézier span the control window is the Bernstein net
            // exactly; the converted net carries rounding on its zeros.
            let zero = self.zero_entries((su, sv));
            let window: Option<Net> = zero.as_ref().map(|_| {
                (su - p..=su)
                    .map(|i| {
                        (sv - q..=sv)
                            .map(|j| self.homogeneous[i * self.num.1 + j])
                            .collect()
                    })
                    .collect()
            });
            let base = window.as_ref().unwrap_or(&span.nets[jet_slot(0, 0)]);
            let rows: Net = base
                .iter()
                .map(|row| bernstein::restrict(row, av, bv))
                .collect();
            let restricted: Net = transpose(&rows)
                .iter()
                .map(|col| bernstein::restrict(col, au, bu))
                .collect();
            // The original entries restricted entry `k` (of `m`)
            // combines with positive coefficients: its blossom's
            // arguments are `a` m−1−k times and `b` k times, or the one
            // point `a ∪ b` n times on a point-width box.
            let reach = |k: usize, m: usize, n: usize, a: Iv, b: Iv| {
                if m == 1 {
                    return if a.hi <= 0.0 && b.hi <= 0.0 {
                        0..=0
                    } else if a.lo >= 1.0 && b.lo >= 1.0 {
                        n..=n
                    } else {
                        0..=n
                    };
                }
                let lo = if b.lo >= 1.0 { k } else { 0 };
                let hi = if a.hi <= 0.0 { k } else { n };
                lo..=hi
            };
            // `restricted` is indexed `[i_v][i_u]` (rows were restricted
            // in v, then transposed and restricted in u).
            for (j, row) in restricted.iter().enumerate() {
                for (i, entry) in row.iter().enumerate() {
                    let w = entry[3];
                    let null = zero.as_ref().is_some_and(|z| {
                        reach(i, row.len(), p, au, bu)
                            .all(|a| reach(j, restricted.len(), q, av, bv).all(|b| z[a][b]))
                    });
                    if null {
                        continue;
                    }
                    if w.lo <= 0.0 {
                        return Err(CertifyDecline::DenominatorContainsZero);
                    }
                    let at = IvVec3::new(entry[0].div(w), entry[1].div(w), entry[2].div(w));
                    out = Some(out.map_or(at, |acc| acc.hull(at)));
                }
            }
        }
        out.ok_or(CertifyDecline::DenominatorContainsZero)
    }

    /// Which Bernstein entries of span `(su, sv)` are exactly zero
    /// (ε-regularized corners), indexed `[i_u][i_v]`: on a Bézier span
    /// (every inner knot at the span's ends) the Bernstein net is the
    /// control window, so its exact-zero control points. `None` on
    /// other spans, where every entry mixes several control points.
    fn zero_entries(&self, (su, sv): (usize, usize)) -> Option<Vec<Vec<bool>>> {
        let (p, q) = self.degree;
        let bezier = |knots: &[f64], s: usize, d: usize| {
            knots[s + 1 - d..=s].iter().all(|&k| k == knots[s])
                && knots[s + 1..=s + d].iter().all(|&k| k == knots[s + 1])
        };
        if !bezier(&self.knots_u, su, p) || !bezier(&self.knots_v, sv, q) {
            return None;
        }
        let num_v = self.num.1;
        Some(
            (su - p..=su)
                .map(|i| {
                    (sv - q..=sv)
                        .map(|j| {
                            self.homogeneous[i * num_v + j]
                                .iter()
                                .all(|c| c.lo == 0.0 && c.hi == 0.0)
                        })
                        .collect()
                })
                .collect(),
        )
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
