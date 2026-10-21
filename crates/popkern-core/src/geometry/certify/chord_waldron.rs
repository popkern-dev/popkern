//! Waldron's bound of a function's deviation from its linear
//! interpolant over a triangle, in the smallest of three metric forms
//! (see `chord_triangle`).

use super::super::interval::Iv;
use super::chord_triangle::{Uv, coord};

/// `(M/2) R²` for the smallest valid `M R²` given bounds `a ≥ |f_uu|`,
/// `b ≥ |f_uv|`, `c ≥ |f_vv|` over triangle `t`: the deviation of `f`
/// from its linear interpolant at `t`'s vertices.
pub(super) fn waldron(t: &[Uv; 3], [a, b, c]: [f64; 3]) -> f64 {
    // Bounds of |D²_d S| = |d1² S_uu + 2 d1 d2 S_uv + d2² S_vv|, all
    // valid, the smallest kept: isotropic, λ_max([[a, b], [b, c]])
    // times |d|²; and diagonal, (a + εb) d1² + (c + b/ε) d2² for any
    // ε > 0 (as 2|d1 d2| ≤ ε d1² + d2²/ε), the metric of the
    // parameters scaled by sqrt(a + εb) and sqrt(c + b/ε). ε = 1, and
    // ε = the triangle's v/u extent ratio, which balances the cross
    // term on long strips (a cylinder's or a far cone's: the isotropic
    // form charges the axial length with the ring curvature).
    let (half_sum, half_diff) = (Iv::point(a + c).scale(0.5), Iv::point(a - c).scale(0.5));
    let m = (half_sum + (half_diff.sqr() + Iv::point(b).sqr()).sqrt()).hi;
    let isotropic = (Iv::point(m) * Iv::point(enclosing_radius_sqr(t)))
        .scale(0.5)
        .hi;
    let extent = |k: usize| {
        let xs = t.iter().map(|q| coord(*q, k));
        xs.clone().fold(f64::NEG_INFINITY, f64::max) - xs.fold(f64::INFINITY, f64::min)
    };
    let ratio = (extent(1) / extent(0)).clamp(1e-6, 1e6);
    let diagonal = |eps: f64| {
        let (eps, b) = (Iv::point(eps), Iv::point(b));
        let su = (Iv::point(a) + eps * b).sqrt().hi;
        let sv = (Iv::point(c) + b / eps).sqrt().hi;
        let scaled = t.map(|q| (q.0 * su, q.1 * sv));
        Iv::point(enclosing_radius_sqr(&scaled)).scale(0.5).hi
    };
    let quad = if ratio.is_finite() {
        diagonal(ratio)
    } else {
        f64::INFINITY
    };
    isotropic.min(diagonal(1.0)).min(quad) * (1.0 + 1e-12)
}

/// Squared radius of the smallest circle holding the triangle.
fn enclosing_radius_sqr(t: &[Uv; 3]) -> f64 {
    let d2 = |a: Uv, b: Uv| (a.0 - b.0).powi(2) + (a.1 - b.1).powi(2);
    let sides = [d2(t[1], t[2]), d2(t[0], t[2]), d2(t[0], t[1])];
    let long = sides.iter().copied().fold(0.0, f64::max);
    let others: f64 = sides.iter().sum::<f64>() - long;
    if long >= others {
        // Right or obtuse: the longest side's midpoint.
        return long / 4.0;
    }
    // Acute: the circumradius, R² = a²b²c² / (16 · area²).
    let area2 = ((t[1].0 - t[0].0) * (t[2].1 - t[0].1) - (t[2].0 - t[0].0) * (t[1].1 - t[0].1))
        .powi(2)
        / 4.0;
    if area2 <= 0.0 {
        return long / 4.0;
    }
    (sides[0] * sides[1] * sides[2] / (16.0 * area2)).max(long / 4.0)
}
