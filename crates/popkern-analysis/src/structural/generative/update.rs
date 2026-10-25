//! The optimality criteria update (top88 eqn 3): each free design
//! variable moves to x·√(−∂c/∂x / (λ ∂V/∂x)), clamped to the move limit
//! and [0, 1], with the Lagrange multiplier λ found by bisection so the
//! filtered design meets the volume target.

use super::filter::DensityFilter;
use super::print::PrintFilter;
use super::projection::Projection;

/// Bisection steps on λ: halving 200 times takes any bracket to the
/// last bit, so the volume is met to round-off.
const BISECTIONS: usize = 200;

/// What one update needs besides the design.
pub struct UpdateInput<'a> {
    /// ∂c/∂x per entry (≤ 0 for a compliance).
    pub dc: &'a [f64],
    /// ∂V/∂x per entry (> 0).
    pub dv: &'a [f64],
    /// Each entry's density floor (kept material).
    pub floors: &'a [f64],
    /// Entry volumes, mm³.
    pub volumes: &'a [f64],
    /// The density filter.
    pub filter: &'a DensityFilter,
    /// The projection after the filter, when on.
    pub projection: Option<Projection>,
    /// The print filter after the projection, when on.
    pub print: Option<&'a PrintFilter>,
    /// Target material volume, mm³.
    pub target: f64,
    /// Largest change of a variable per update.
    pub move_limit: f64,
}

/// The material volume of design `x`: Σ v_e x̄_e.
pub fn volume(
    x: &[f64],
    floors: &[f64],
    volumes: &[f64],
    filter: &DensityFilter,
    projection: Option<Projection>,
    print: Option<&PrintFilter>,
) -> f64 {
    physical(x, floors, filter, projection, print)
        .iter()
        .zip(volumes)
        .map(|(p, v)| p * v)
        .sum()
}

/// The physical densities of design `x`: the filtered density,
/// projected when the projection is on, printed when the print filter
/// is on, at least the entry's floor.
pub fn physical(
    x: &[f64],
    floors: &[f64],
    filter: &DensityFilter,
    projection: Option<Projection>,
    print: Option<&PrintFilter>,
) -> Vec<f64> {
    let mut p = filter.apply(x);
    if let Some(pr) = projection {
        for p in p.iter_mut() {
            *p = pr.apply(*p);
        }
    }
    match print {
        Some(f) => f.apply(&p, floors).densities,
        None => {
            for (p, &f) in p.iter_mut().zip(floors) {
                *p = p.max(f);
            }
            p
        }
    }
}

/// The updated design and its multiplier.
pub fn optimality_criteria(x: &[f64], input: &UpdateInput<'_>) -> (Vec<f64>, f64) {
    let step = |lambda: f64| -> Vec<f64> {
        x.iter()
            .enumerate()
            .map(|(e, &xe)| {
                let lo = (xe - input.move_limit).max(input.floors[e]);
                let hi = (xe + input.move_limit).min(1.0).max(lo);
                let ratio = if input.dv[e] > 0.0 && lambda > 0.0 {
                    (-input.dc[e]).max(0.0) / (input.dv[e] * lambda)
                } else {
                    0.0
                };
                (xe * ratio.sqrt()).clamp(lo, hi)
            })
            .collect()
    };
    let over = |lambda: f64| {
        volume(
            &step(lambda),
            input.floors,
            input.volumes,
            input.filter,
            input.projection,
            input.print,
        ) > input.target
    };
    // Bracket λ: grow the upper end until the design fits.
    let (mut l1, mut l2) = (0.0, 1.0);
    while over(l2) && l2 < 1e300 {
        l1 = l2;
        l2 *= 2.0;
    }
    for _ in 0..BISECTIONS {
        let mid = 0.5 * (l1 + l2);
        if mid <= l1 || mid >= l2 {
            break;
        }
        if over(mid) {
            l1 = mid;
        } else {
            l2 = mid;
        }
    }
    (step(l2), l2)
}
