//! Wrapping values by a period: a curve's or a chart's angle brought into
//! a window. Named so callers state what they do (the chart service owns
//! a surface's periods and asks here for the arithmetic; a circle's
//! curve parameter wraps here too).

/// `x` brought into `[lo, lo + period)` by whole periods.
pub fn wrap_from(x: f64, lo: f64, period: f64) -> f64 {
    lo + (x - lo).rem_euclid(period)
}

/// The offset from `b` to `a` reduced into `[0, period)`.
pub fn offset_mod(a: f64, b: f64, period: f64) -> f64 {
    (a - b).rem_euclid(period)
}

/// `x` reduced into `[0, period)` (exactly `x.rem_euclid(period)`).
pub fn modulo(x: f64, period: f64) -> f64 {
    x.rem_euclid(period)
}
