//! The position-only enclosure of a NURBS jet source (PK-904): where the
//! weight function vanishes, the jets' quotient declines but the hull
//! of the positive-weight net entries still encloses.

use super::super::surface_jet::ParamMap;
use super::*;

impl NurbsJetSource {
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

    /// An upper bound of `|S(u, v) − L(u, v)|` over the box, `L` the
    /// affine or bilinear pairing `l`, from the nets alone (PK-909). On
    /// each span the restricted homogeneous net and the net of `W·L`
    /// (`W` the weight function) are taken to one degree higher in each
    /// direction, where both are exact; `S − L` is then the convex
    /// combination, with weights `b_k W_k / W`, of `D_k / W_k`, `D` the
    /// difference of their numerators. Entries of exactly zero weight
    /// (an activated corner's) have `D_k = 0` and contribute nothing, so
    /// the bound holds at the corner; it shrinks with the box like the
    /// Taylor bound, where the position hull stalls at first order.
    ///
    /// # Errors
    ///
    /// As [`NurbsJetSource::position_hull`].
    pub fn difference_hull(&mut self, u: Iv, v: Iv, l: &ParamMap) -> Result<f64, CertifyDecline> {
        let (p, q) = self.degree;
        let domain = |knots: &[f64], degree: usize, num: usize| (knots[degree], knots[num]);
        let (du, dv) = (
            domain(&self.knots_u, p, self.num.0),
            domain(&self.knots_v, q, self.num.1),
        );
        if u.lo < du.0 || u.hi > du.1 || v.lo < dv.0 || v.hi > dv.1 {
            return Err(CertifyDecline::SpanStraddle);
        }
        let mut worst = 0.0_f64;
        for (su, sv) in self.spans_touching(u, v)? {
            let (restricted, null, (ub, vb)) = self.restricted((su, sv), u, v)?;
            // L over the clipped box in its local (s, t) ∈ [0, 1]²:
            // with α = α0 + s·dα, β = β0 + t·dβ the map's own
            // coordinates, L = a + s·b + t·c + s·t·d.
            let (a0, da) = (
                (Iv::point(ub.0) - l.u0) / l.hu,
                (Iv::point(ub.1) - Iv::point(ub.0)) / l.hu,
            );
            let (b0, db) = (
                (Iv::point(vb.0) - l.v0) / l.hv,
                (Iv::point(vb.1) - Iv::point(vb.0)) / l.hv,
            );
            let a = l.p + l.q.scale(a0) + l.r.scale(b0) + l.t.scale(a0 * b0);
            let b = (l.q + l.t.scale(b0)).scale(da);
            let c = (l.r + l.t.scale(a0)).scale(db);
            let d = l.t.scale(da * db);
            let (nv, nu) = (restricted.len(), restricted.first().map_or(0, Vec::len));
            let weight = |i: usize, j: usize| restricted[j][i][3];
            // Elevated entry (i, j) of degree (nu, nv) from the
            // degree-(nu − 1, nv − 1) entries `f(i', j')`.
            let elevate = |f: &dyn Fn(usize, usize) -> Coef<4>, i: usize, j: usize| -> Coef<4> {
                let mut acc = [Iv::point(0.0); 4];
                for (ii, ai) in [(i.wrapping_sub(1), i), (i, nu - i)] {
                    if ii >= nu || ai == 0 {
                        continue;
                    }
                    for (jj, aj) in [(j.wrapping_sub(1), j), (j, nv - j)] {
                        if jj >= nv || aj == 0 {
                            continue;
                        }
                        let k = Iv::point((ai * aj) as f64 / (nu * nv) as f64);
                        let e = f(ii, jj);
                        for m in 0..4 {
                            acc[m] = acc[m] + e[m] * k;
                        }
                    }
                }
                acc
            };
            let entry = |i: usize, j: usize| restricted[j][i];
            let w_only = |i: usize, j: usize| [weight(i, j); 4];
            // s·W and t·W, each one degree up in its own direction and
            // elevated in the other.
            let s_w = |i: usize, j: usize| -> Iv {
                let mut acc = Iv::point(0.0);
                if i >= 1 {
                    for (jj, aj) in [(j.wrapping_sub(1), j), (j, nv - j)] {
                        if jj < nv && aj != 0 {
                            acc = acc
                                + weight(i - 1, jj) * Iv::point((i * aj) as f64 / (nu * nv) as f64);
                        }
                    }
                }
                acc
            };
            let t_w = |i: usize, j: usize| -> Iv {
                let mut acc = Iv::point(0.0);
                if j >= 1 {
                    for (ii, ai) in [(i.wrapping_sub(1), i), (i, nu - i)] {
                        if ii < nu && ai != 0 {
                            acc = acc
                                + weight(ii, j - 1) * Iv::point((ai * j) as f64 / (nu * nv) as f64);
                        }
                    }
                }
                acc
            };
            // s·t·W, one degree up in both.
            let st_w = |i: usize, j: usize| -> Iv {
                if i >= 1 && j >= 1 {
                    weight(i - 1, j - 1) * Iv::point((i * j) as f64 / (nu * nv) as f64)
                } else {
                    Iv::point(0.0)
                }
            };
            for i in 0..=nu {
                for j in 0..=nv {
                    // Exactly zero when every entry it combines is.
                    let zero = [
                        (i.wrapping_sub(1), j.wrapping_sub(1)),
                        (i.wrapping_sub(1), j),
                        (i, j.wrapping_sub(1)),
                        (i, j),
                    ]
                    .iter()
                    .filter(|&&(ii, jj)| ii < nu && jj < nv)
                    .all(|&(ii, jj)| null[jj][ii]);
                    if zero {
                        continue;
                    }
                    let r = elevate(&entry, i, j);
                    let w = elevate(&w_only, i, j)[3];
                    if w.lo <= 0.0 {
                        return Err(CertifyDecline::DenominatorContainsZero);
                    }
                    let wl =
                        a.scale(w) + b.scale(s_w(i, j)) + c.scale(t_w(i, j)) + d.scale(st_w(i, j));
                    let diff = IvVec3::new(r[0], r[1], r[2]) - wl;
                    worst = worst.max((diff.norm() / w).hi);
                }
            }
        }
        Ok(worst)
    }

    /// Span `(su, sv)`'s homogeneous net restricted to the box (indexed
    /// `[i_v][i_u]`), which of its entries are exactly zero, and the box
    /// clipped to the span.
    #[allow(clippy::type_complexity)]
    fn restricted(
        &mut self,
        (su, sv): (usize, usize),
        u: Iv,
        v: Iv,
    ) -> Result<(Net, Vec<Vec<bool>>, ((f64, f64), (f64, f64))), CertifyDecline> {
        let (p, q) = self.degree;
        if !self.cache.contains_key(&(su, sv)) {
            let nets = self.span_nets(su, sv);
            self.cache.insert((su, sv), nets);
        }
        let span = self
            .cache
            .get(&(su, sv))
            .ok_or(CertifyDecline::MalformedNurbs)?;
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
        let clip = |x: Iv, s: (f64, f64)| {
            let lo = x.lo.max(s.0);
            (lo, x.hi.min(s.1).max(lo))
        };
        let (cu, cv) = (clip(u, span.u), clip(v, span.v));
        let ((au, bu), (av, bv)) = (
            (unit(cu.0, span.u), unit(cu.1, span.u)),
            (unit(cv.0, span.v), unit(cv.1, span.v)),
        );
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
        let null = restricted
            .iter()
            .enumerate()
            .map(|(j, row)| {
                (0..row.len())
                    .map(|i| {
                        zero.as_ref().is_some_and(|z| {
                            reach(i, row.len(), p, au, bu)
                                .all(|a| reach(j, restricted.len(), q, av, bv).all(|b| z[a][b]))
                        })
                    })
                    .collect()
            })
            .collect();
        Ok((restricted, null, (cu, cv)))
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
}
