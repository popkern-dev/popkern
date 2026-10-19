//! The contact solver's iteration window along a walk (split from
//! `walker.rs` under the 500-line cap): the periodic re-centring every
//! sample after the first gets, and the first sample's exception when
//! its seed sits on the end of a full-period window (PK-766, PK-888).

/// Recenter the contact solver's UV range around `current_uv` for
/// periodic dimensions, leaving non-periodic dimensions untouched.
///
/// For a periodic dimension with period `T`, the new bounds are
/// `(current_u - 2·T, current_u + 2·T)` — four periods wide,
/// centered on the previous solution. Newton's natural step from a
/// good warm start is at most one period (typically a small fraction
/// of one), so the bounds are unreachable in normal use.
///
/// The 4×period factor (rather than the obvious 2×period centered
/// window) is deliberate: with a 2-period-wide window, Newton can
/// converge to a periodically-equivalent point that lands *exactly*
/// on one of the new boundaries (the same 3D point as `current_u`
/// reached by adding/subtracting one period). The post-convergence
/// `uv_on_boundary` check in the contact solver then fires and
/// reports a spurious `OutOfBounds`. Widening to 4×period puts the
/// boundaries 2 periods away from `current_u`, which is far beyond
/// any plausible step.
///
/// Non-periodic dimensions inherit `base_range` unchanged so the
/// walker still respects real face bounds (e.g. the cone's axial v
/// parameter).
pub(crate) fn recenter_range(
    base_range: ((f64, f64), (f64, f64)),
    current_uv: (f64, f64),
    periodicity: (Option<f64>, Option<f64>),
) -> ((f64, f64), (f64, f64)) {
    let (period_u, period_v) = periodicity;
    let u_range = match period_u {
        Some(p) => (current_uv.0 - 2.0 * p, current_uv.0 + 2.0 * p),
        None => base_range.0,
    };
    let v_range = match period_v {
        Some(p) => (current_uv.1 - 2.0 * p, current_uv.1 + 2.0 * p),
        None => base_range.1,
    };
    (u_range, v_range)
}

/// True when `uv` lies on an end of a periodic dimension whose window
/// spans a whole period: there the window end is an arbitrary cut
/// through the surface, not a face bound.
pub(crate) fn seed_on_full_period_end(
    range: ((f64, f64), (f64, f64)),
    uv: (f64, f64),
    (period_u, period_v): (Option<f64>, Option<f64>),
) -> bool {
    let on_end = |(lo, hi): (f64, f64), x: f64, period: Option<f64>| {
        period.is_some_and(|p| {
            let slack = p * 1e-9;
            hi - lo >= p - slack && ((x - lo).abs() <= slack || (hi - x).abs() <= slack)
        })
    };
    on_end(range.0, uv.0, period_u) || on_end(range.1, uv.1, period_v)
}
