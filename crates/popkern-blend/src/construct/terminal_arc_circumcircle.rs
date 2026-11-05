//! Circumcircle of three points in a plane chart, for fitting a terminal arc through sampled contacts.

/// Circumcircle through three 2-D points; returns (cu, cv, r).
pub(super) fn circumcircle_2d(
    a: (f64, f64),
    b: (f64, f64),
    c: (f64, f64),
) -> Option<(f64, f64, f64)> {
    let d = 2.0 * (a.0 * (b.1 - c.1) + b.0 * (c.1 - a.1) + c.0 * (a.1 - b.1));
    if d.abs() < 1e-30 {
        return None;
    }
    let a2 = a.0 * a.0 + a.1 * a.1;
    let b2 = b.0 * b.0 + b.1 * b.1;
    let c2 = c.0 * c.0 + c.1 * c.1;
    let cu = (a2 * (b.1 - c.1) + b2 * (c.1 - a.1) + c2 * (a.1 - b.1)) / d;
    let cv = (a2 * (c.0 - b.0) + b2 * (a.0 - c.0) + c2 * (b.0 - a.0)) / d;
    let r = ((a.0 - cu).powi(2) + (a.1 - cv).powi(2)).sqrt();
    Some((cu, cv, r))
}
