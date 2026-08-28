//! K-P 2015 bi-6 µ coefficient tables.
//!
//! Source: `docs/journal/polar-kp-math.md` §2.3 (transcribed from
//! [KP15] Appendix). These vectors appear in the BB-coefficient
//! assignment rules (S1)–(S4) of the K-P 2015 bi-6 polar cap
//! construction — see math summary §2.2 for the rules that
//! consume them.
//!
//! Each vector is a pure rational function of `(c, τ)` where
//! `c := cos(2π / n)` and `τ` is the per-valence shape parameter
//! (see [`super::kp15_bi6_reparam::tau_for_valence`]). All return
//! `f64` arrays of fixed length; the size matches the paper's
//! indexing (e.g. µ̇ has 5 entries `µ̇_1..µ̇_5` stored as indices
//! `[0..4]`; callers use `mu_dot[i - 1]` to get the paper's
//! `µ̇_i`).
//!
//! All entries are well-defined and finite for the standard
//! valences `n ∈ {3, 5, 6, 7, 8, 9, 10, >10}` where `c ≠ 0` and
//! `τ > 0`. Divisions by `c`, `τ`, `τ²`, `τ³` appear in the tables
//! and are documented at the relevant call sites.

use super::kp15_bi6_reparam::Bi6Reparams;

/// `µ̇_{1..5}` per [KP15 Appendix] — the coefficient vector used
/// in the (S2) rule that determines `p̀_{21}` from the shared-
/// edge data (math summary §2.2 rule (S2)).
pub fn mu_dot(p: &Bi6Reparams) -> [f64; 5] {
    let c = p.c;
    [2.0 * c / 5.0, 2.0 - c, c, -2.0 * c / 5.0, c / 15.0]
}

/// `µ̈_{1..8}` per [KP15 Appendix] — the coefficient vector used
/// in the (S2) rule that determines `p̀_{12}`.
pub fn mu_ddot(p: &Bi6Reparams) -> [f64; 8] {
    let c = p.c;
    let tau = p.tau;
    let tau2 = tau * tau;
    [
        -4.0 * c * tau2 / 5.0,
        4.0 * (1.0 - c + c * c / 5.0) + 24.0 * c * (c - 1.0) / 5.0 * tau + 4.0 * c * c * tau2 / 5.0,
        24.0 * c / 5.0 * tau - 4.0,
        2.0 * c / 3.0 * (11.0 - 3.0 * c - 6.0 * c * tau),
        -4.0 * c,
        2.0 * c * c,
        -4.0 * c * c / 5.0,
        2.0 * c * c / 15.0,
    ]
}

/// `⃛µ_{1..11}` per [KP15 Appendix] — coefficient vector used in
/// the (S3) rule for `p̀_{22}` (the `n × n` circulant system).
pub fn mu_tdot(p: &Bi6Reparams) -> [f64; 11] {
    let c = p.c;
    let tau = p.tau;
    let tau2 = tau * tau;
    let tau3 = tau2 * tau;
    [
        4.0 * c * tau3 / 25.0,
        4.0 * c / 25.0 * (5.0 - 3.0 * c * tau + 6.0 * (1.0 - c) * tau2 - c * tau3),
        -24.0 * c * tau2 / 25.0,
        4.0 - 2.0 * c + 6.0 * c * (c - 2.0) / 5.0 * tau + 4.0 * c * c / 5.0 * tau2,
        12.0 * c / 5.0 * tau - 4.0,
        2.0 * c / 5.0 * (5.0 - 3.0 * c * tau),
        4.0 * c / 25.0 * (3.0 * c * tau - 5.0),
        2.0 * c / 75.0 * (5.0 - 3.0 * c * tau),
        -3.0 * c / (5.0 * tau),
        6.0 * c * (5.0 - c * tau) / (125.0 * tau2),
        c * (6.0 * c * tau - 5.0) / (125.0 * tau3),
    ]
}

/// `µ̃_{1..10}` per [KP15 Appendix] — coefficient vector used in
/// the (S4) rule for `p̀_{31}` and `ṕ_{31}` (inner CPs near the
/// shared edge). Divisions by `c` appear in entries 2, 3, 10;
/// always safe for `n ∉ {4}` (our polar path is `n ≥ 5` so
/// `c ≠ 0`; `n = 3` gives `c = −1/2 ≠ 0`).
pub fn mu_tilde(p: &Bi6Reparams) -> [f64; 10] {
    let c = p.c;
    let tau = p.tau;
    let tau2 = tau * tau;
    let tau3 = tau2 * tau;
    [
        -11.0 / 20.0,
        (11.0 * c - 10.0) / (8.0 * c),
        5.0 / (4.0 * c),
        1.0 / 8.0 - c / 4.0,
        3.0 * c / 10.0 - 1.0 / 20.0,
        1.0 / 120.0 - c / 20.0,
        3.0 / (8.0 * tau),
        (3.0 * c * tau - 15.0) / (100.0 * tau2),
        (5.0 - 6.0 * c * tau) / (200.0 * tau3),
        -5.0 / (16.0 * c),
    ]
}

/// `µ̄_{1..10}` per [KP15 Appendix] — coefficient vector used in
/// the (S4) rule for `ṕ_{32}`.
pub fn mu_bar(p: &Bi6Reparams) -> [f64; 10] {
    let c = p.c;
    let tau = p.tau;
    let tau2 = tau * tau;
    let tau3 = tau2 * tau;
    [
        1.0,
        5.0 * (2.0 - c) / (2.0 * c),
        -5.0 / c,
        5.0 / 2.0,
        -1.0,
        1.0 / 6.0,
        3.0 / (10.0 * tau) * (c * tau - 5.0),
        3.0 / (125.0 * tau2) * (25.0 - 20.0 * c * tau + c * c * tau2),
        -1.0 / (250.0 * tau3) * (25.0 - 45.0 * c * tau + 6.0 * c * c * tau2),
        5.0 / (4.0 * c),
    ]
}

/// `w₁, w₂, w₃` weights from [KP15 eq (11)] — scalar multipliers
/// in the last three formulas of rule (S4) that determine
/// `p̀_{41}` and `ṕ_{41}`. `w₁` has a `1/τ` factor; safe for
/// `τ > 0`.
pub fn w_scalars(p: &Bi6Reparams) -> (f64, f64, f64) {
    let c = p.c;
    let tau = p.tau;
    (-c / (50.0 * tau), c / 50.0, c / 15.0)
}

#[cfg(test)]
#[path = "kp15_bi6_mu_tests.rs"]
mod tests;
