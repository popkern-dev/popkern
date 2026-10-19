//! Periodic NURBS rings (PK-890, `docs/architecture/closed-nurbs-rings.md`).
//!
//! A NURBS surface is periodic in a direction when it is stored in
//! periodic form there: an unclamped knot vector whose spacing repeats
//! one period on, and control rows whose last `p` repeat the first `p`
//! (`p` the degree). Such a ring is smooth across its seam to the
//! continuity of its knots, like a cylinder's lateral face. The test
//! reads the representation alone, so periodicity survives
//! serialisation and copies with no flag.
//!
//! A clamped surface whose end rows merely coincide is CLOSED but not
//! periodic: its seam is a chart boundary, as OpenCASCADE separates
//! `IsUClosed` from `IsUPeriodic`. Coil-spring tubes and revolved
//! NURBS walls are clamped rings and keep their open-chart handling.

use super::NurbsSurface;

/// Roundoff budget, in ulps of the net's coordinate scale and of the
/// period, for two control points or two knot gaps to count as equal.
///
/// Periodic form is a statement about the representation, not a
/// modelling tolerance: a producer writes the wrapped rows and knots
/// by copying them, so they agree to the last few ulps or not at all.
const PERIODIC_ULPS: f64 = 64.0;

impl NurbsSurface {
    /// The period in `u` when the surface is stored in periodic form in
    /// `u`, else `None`. See the module docs.
    pub fn u_period(&self) -> Option<f64> {
        let grid = |i_u: usize, i_v: usize| (i_u, i_v);
        self.period_along(self.degree_u, &self.knots_u, self.num_u, self.num_v, grid)
    }

    /// The period in `v` when the surface is stored in periodic form in
    /// `v`, else `None`.
    pub fn v_period(&self) -> Option<f64> {
        let grid = |i_v: usize, i_u: usize| (i_u, i_v);
        self.period_along(self.degree_v, &self.knots_v, self.num_v, self.num_u, grid)
    }

    /// The periodic-form test along one direction. `grid(i, j)` maps
    /// the index `i` along the tested direction and `j` across it to
    /// `(i_u, i_v)`.
    fn period_along(
        &self,
        degree: usize,
        knots: &[f64],
        count: usize,
        across: usize,
        grid: impl Fn(usize, usize) -> (usize, usize),
    ) -> Option<f64> {
        // `count - degree` distinct rows, at least two spans of them.
        let distinct = count.checked_sub(degree)?;
        if degree == 0 || distinct < 2 || knots.len() != count + degree + 1 {
            return None;
        }
        let period = knots[count] - knots[degree];
        if !(period > 0.0 && period.is_finite()) {
            return None;
        }
        // Knot spacing repeats one period on: k[i + distinct] = k[i] + T
        // for every knot that has a partner.
        let knot_gap = PERIODIC_ULPS * f64::EPSILON * period.max(knots[count].abs());
        let knots_wrap = (0..knots.len() - distinct)
            .all(|i| (knots[i + distinct] - knots[i] - period).abs() <= knot_gap);
        if !knots_wrap {
            return None;
        }
        let scale = self
            .control_points
            .iter()
            .map(|p| p.x.abs().max(p.y.abs()).max(p.z.abs()))
            .fold(1.0_f64, f64::max);
        let point_gap = PERIODIC_ULPS * f64::EPSILON * scale;
        let rows_wrap = (0..degree).all(|i| {
            (0..across).all(|j| {
                let (a, b) = (grid(i, j), grid(i + distinct, j));
                let (wa, wb) = (self.wt(a.0, a.1), self.wt(b.0, b.1));
                self.cp(a.0, a.1).distance_to(self.cp(b.0, b.1)) <= point_gap
                    && (wa - wb).abs() <= PERIODIC_ULPS * f64::EPSILON * wa.abs().max(wb.abs())
            })
        });
        rows_wrap.then_some(period)
    }

    /// `(u, v)` moved into the knot domain along each periodic
    /// direction.
    ///
    /// Parameters already inside the domain return unchanged, so
    /// in-domain evaluation is bit-identical to an open chart's and the
    /// periodic-form test runs only for out-of-domain queries. An open
    /// direction is left alone: evaluation there extrapolates its end
    /// span as before.
    pub(super) fn wrap_into_domain(&self, u: f64, v: f64) -> (f64, f64) {
        (
            wrap_axis(u, &self.knots_u, self.degree_u, self.num_u, || {
                self.u_period()
            }),
            wrap_axis(v, &self.knots_v, self.degree_v, self.num_v, || {
                self.v_period()
            }),
        )
    }
}

fn wrap_axis(
    t: f64,
    knots: &[f64],
    degree: usize,
    count: usize,
    period: impl FnOnce() -> Option<f64>,
) -> f64 {
    let (Some(&lo), Some(&hi)) = (knots.get(degree), knots.get(count)) else {
        return t;
    };
    if (lo..=hi).contains(&t) || !t.is_finite() {
        return t;
    }
    match period() {
        Some(period) => lo + (t - lo).rem_euclid(period),
        None => t,
    }
}

#[cfg(test)]
#[path = "periodic_form_tests.rs"]
mod tests;
