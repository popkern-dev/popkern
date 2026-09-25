//! Interval Bernstein and B-spline primitives.
//!
//! Coefficients are `[Iv; D]` homogeneous tuples (`D = 4` for
//! surfaces, `D = 3` for pcurves). Everything runs in outward-rounded
//! interval arithmetic, so each output encloses the exact coefficient
//! of the exact input.
//!
//! Enclosure principle: a polynomial restricted to a parameter
//! interval, written in the Bernstein basis over that interval, lies
//! in the convex hull of its coefficients, so the coefficient hull is
//! a rigorous enclosure of its range there.

use super::interval::Iv;
use std::ops::Div;

/// A homogeneous coefficient tuple.
pub type Coef<const D: usize> = [Iv; D];

fn lerp<const D: usize>(a: &Coef<D>, b: &Coef<D>, t: Iv) -> Coef<D> {
    let one_minus = Iv::point(1.0) - t;
    let mut out = *a;
    for k in 0..D {
        out[k] = one_minus * a[k] + t * b[k];
    }
    out
}

/// Blossom of the B-spline segment on knot span `span` (so
/// `knots[span] ≤ knots[span + 1]`, with `ctrl` the `degree + 1`
/// control points `span − degree ..= span`) at the arguments `args`
/// (length `degree`), by de Boor's algorithm with one argument per
/// level.
pub fn bspline_blossom<const D: usize>(
    degree: usize,
    knots: &[f64],
    span: usize,
    ctrl: &[Coef<D>],
    args: &[f64],
) -> Coef<D> {
    let mut work: Vec<Coef<D>> = ctrl.to_vec();
    for (r, &t) in args.iter().enumerate().take(degree) {
        let level = r + 1;
        for i in (level..=degree).rev() {
            let k = span - degree + i;
            let (left, right) = (knots[k], knots[k + degree + 1 - level]);
            let alpha = (Iv::point(t) - Iv::point(left)).div(Iv::point(right) - Iv::point(left));
            work[i] = lerp(&work[i - 1], &work[i], alpha);
        }
    }
    work[degree]
}

/// Bernstein coefficients of the B-spline segment on knot span `span`
/// over `[knots[span], knots[span + 1]]`: blossom values at
/// `(a^{degree − j}, b^{j})`.
pub fn bspline_span_to_bernstein<const D: usize>(
    degree: usize,
    knots: &[f64],
    span: usize,
    ctrl: &[Coef<D>],
) -> Vec<Coef<D>> {
    let (a, b) = (knots[span], knots[span + 1]);
    (0..=degree)
        .map(|j| {
            let args: Vec<f64> = (0..degree)
                .map(|r| if r < degree - j { a } else { b })
                .collect();
            bspline_blossom(degree, knots, span, ctrl, &args)
        })
        .collect()
}

/// Bernstein coefficients of the derivative of a degree-`n`
/// polynomial on an interval of width `span_width`:
/// `n · (c_{i+1} − c_i) / span_width`.
pub fn derivative_net<const D: usize>(net: &[Coef<D>], span_width: Iv) -> Vec<Coef<D>> {
    let n = net.len().saturating_sub(1);
    let factor = Iv::point(n as f64).div(span_width);
    net.windows(2)
        .map(|w| {
            let mut out = w[0];
            for k in 0..D {
                out[k] = (w[1][k] - w[0][k]) * factor;
            }
            out
        })
        .collect()
}

/// De Casteljau at a (possibly interval) local parameter `t`,
/// returning the value and the left and right sub-nets.
fn de_casteljau_split<const D: usize>(net: &[Coef<D>], t: Iv) -> (Vec<Coef<D>>, Vec<Coef<D>>) {
    let n = net.len();
    let mut work = net.to_vec();
    let mut left = Vec::with_capacity(n);
    let mut right = Vec::with_capacity(n);
    left.push(work[0]);
    right.push(work[n - 1]);
    for level in 1..n {
        for i in 0..n - level {
            work[i] = lerp(&work[i], &work[i + 1], t);
        }
        left.push(work[0]);
        right.push(work[n - 1 - level]);
    }
    right.reverse();
    (left, right)
}

/// Enclosure of the net's polynomial at local parameter `t`.
pub fn de_casteljau_point<const D: usize>(net: &[Coef<D>], t: Iv) -> Coef<D> {
    let mut work = net.to_vec();
    let n = work.len();
    for level in 1..n {
        for i in 0..n - level {
            work[i] = lerp(&work[i], &work[i + 1], t);
        }
    }
    work[0]
}

/// Bernstein net of the polynomial restricted to the local interval
/// `[alpha, beta]` (endpoints given as enclosures). When the two
/// enclosures are identical the restriction is a point, and a single
/// coefficient (the point value) is returned.
pub fn restrict<const D: usize>(net: &[Coef<D>], alpha: Iv, beta: Iv) -> Vec<Coef<D>> {
    if alpha == beta || net.len() == 1 {
        return vec![de_casteljau_point(net, alpha.hull(beta))];
    }
    // [0, 1] → [0, β] (left part), then [0, β] → [α/β, 1] (right).
    let (left, _) = de_casteljau_split(net, beta);
    let ratio = alpha.div(beta);
    if !ratio.is_finite() {
        // β can reach 0 only when α does too; restrict the other way.
        let (_, right) = de_casteljau_split(net, alpha);
        let ratio = (beta - alpha).div(Iv::point(1.0) - alpha);
        let (inner, _) = de_casteljau_split(&right, ratio);
        return inner;
    }
    let (_, inner) = de_casteljau_split(&left, ratio);
    inner
}

/// Component-wise hull of every coefficient.
pub fn hull<const D: usize>(net: &[Coef<D>]) -> Coef<D> {
    let mut out = net[0];
    for c in net.iter().skip(1) {
        for k in 0..D {
            out[k] = out[k].hull(c[k]);
        }
    }
    out
}

/// Local parameter enclosure `(x − a) / (b − a)` of a float `x` in a
/// knot span `[a, b]`.
pub fn local_param(x: f64, a: f64, b: f64) -> Iv {
    (Iv::point(x) - Iv::point(a)).div(Iv::point(b) - Iv::point(a))
}
