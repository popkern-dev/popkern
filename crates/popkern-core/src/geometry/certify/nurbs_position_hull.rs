//! The position-only enclosure of a NURBS jet source (PK-904): where the
//! weight function vanishes, the jets' quotient declines but the hull
//! of the positive-weight net entries still encloses.

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
