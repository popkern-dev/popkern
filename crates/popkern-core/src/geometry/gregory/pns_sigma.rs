//! Phase IV PnS — Central σ reparameterization (Appendix B).
//!
//! σ is a bi-4 map `[0, 1]² → R²` that distributes the central
//! quadratic expansion of the n-gon cap across the N sectors in a
//! way that achieves **G² continuity at the pole**. The construction
//! is rotationally symmetric with flip-symmetry of each sector
//! about its diagonal.
//!
//! # Paper reference
//!
//! `docs/papers/14bi234.pdf` Appendix B. Inspired by the conformal
//! map `z^(4/n), z ∈ C` — σ is a smooth polynomial approximation
//! that gives good iso-line distribution in each sector.
//!
//! # What we export
//!
//! For the downstream cap assembly, the only piece of σ that
//! matters is its **order-≤2 Taylor jet at `(0, 0)`** — i.e., the
//! value `σ(0, 0)` plus six partials up to second order. This jet
//! is composed with the undetermined quadratic map `a(s, t)` to
//! produce the central quadratic expansion of each sector's bi-4
//! patch. Higher-order coefficients of σ (`σ_{21}` etc.) affect
//! the cap only via the F_κ optimisation of the remaining free
//! BB-coefs — they're not part of the pinned central jet.
//!
//! Accordingly, [`SigmaCenterJet`] packages the 6 jet values and
//! [`sigma_center_jet`] builds it from the paper defaults.
//!
//! # σ_{01}, σ_{02} derivation
//!
//! The paper gives σ_{00}, σ_{10}, σ_{11}, σ_{20} explicitly. The
//! remaining order-≤2 coefs σ_{01}, σ_{02} follow from flip-symmetry
//! of the sector about its diagonal (reflection about the sector
//! bisector in image space):
//!   σ_{01} = (d_1 cos α, d_1 sin α)
//!   σ_{02} = (σ_{20,x} cos α, σ_{20,x} sin α)
//!
//! These pair-match each sector's σ_{10}, σ_{20} under the flip.
//! Higher-order σ BB-coefs (σ_{21} etc.) affect only the cap's
//! free outer BB-coefs via F_κ minimisation.
//!
//! # Verifying σ against eq (2)
//!
//! An earlier iteration of this module carried a cross-validation
//! test (`central_expansion_satisfies_eq2_g1_constraint`) that
//! plugged `a∘σ` for adjacent sectors into eq (2) and checked
//! `cur_{10}`, `cur_{20}` matched the chain-rule outputs. That
//! test was deleted 2026-04-24 — it incorrectly assumed `nxt`'s
//! central coefs come from `a∘σ^{nxt}` with σ^{nxt} = R_α · σ^{cur}.
//! Per the KK-Peters 2015 paper, sector patches do NOT relate by
//! rotation: on the shared seam `p^k_{i,0} = p^{k+1}_{i,0}`
//! literally, and `p^{k+1}_{i,1}` is independent data coming from
//! the NEXT sector's own β∘b expansion. Eq (2) is satisfied by the
//! WHOLE cap construction (all N sectors' a∘σ + boundary data),
//! not by σ alone — testing it in isolation is ill-posed.

#![allow(dead_code)]

/// Order-≤2 Taylor jet of σ at `(u, v) = (0, 0)` — the central
/// quadratic expansion of the sector-local σ reparameterization.
///
/// All values live in `R²` (per the 2D co-domain of σ). See
/// [`super::pns_sigma`] module docs for the role of this jet.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SigmaCenterJet {
    /// σ(0, 0) — always `(0, 0)` by construction.
    pub s00: [f64; 2],
    /// ∂σ/∂u at `(0, 0)`.
    pub ds_u: [f64; 2],
    /// ∂σ/∂v at `(0, 0)`.
    pub ds_v: [f64; 2],
    /// ∂²σ/∂u² at `(0, 0)`.
    pub ds_uu: [f64; 2],
    /// ∂²σ/∂v² at `(0, 0)`.
    pub ds_vv: [f64; 2],
    /// ∂²σ/∂u∂v at `(0, 0)`.
    pub ds_uv: [f64; 2],
}

/// Paper default scaling `d_1^def = (0.24n − 0.75)/(n − 3.16)`
/// (Appendix B). `d_1` is the x-coordinate of `σ_{10}` (σ's first
/// CP on the sector's lower edge).
pub fn d1_def(n: usize) -> f64 {
    let n_f = n as f64;
    (0.24 * n_f - 0.75) / (n_f - 3.16)
}

/// Paper default scaling `d_2^def = (0.44n − 0.75)/((n − 4) + 2.02·√2)`
/// (Appendix B). `d_2` controls `σ_{11}` along the sector's
/// bisector at the central point.
pub fn d2_def(n: usize) -> f64 {
    let n_f = n as f64;
    (0.44 * n_f - 0.75) / ((n_f - 4.0) + 2.02 * std::f64::consts::SQRT_2)
}

/// Build the order-≤2 jet of σ at `(0, 0)` from explicit scaling
/// parameters `d_1`, `d_2` for valence `n`. See [`sigma_center_jet_default`]
/// for the paper defaults.
///
/// # Construction (from Appendix B + sector flip-symmetry)
///
/// Let `α := 2π/n`, `c := cos α`.
///
/// The paper gives four order-≤2 BB-coefs explicitly:
/// ```text
///   σ_{00} = (0, 0)
///   σ_{10} = (d_1, 0)                      — on the sector's lower edge
///   σ_{11} = (d_2 · cos(α/2), d_2 · sin(α/2))   — on the bisector
///   σ_{20} = ((1 − 4c/3) d_1
///            + (4 cos(α/2) / (3c)) d_2, 0) — on the lower edge
/// ```
///
/// The remaining two order-≤2 coefs follow from flip-symmetry of the
/// sector about its diagonal (which in image space is reflection
/// about the bisector line `angle = α/2`):
/// ```text
///   σ_{01} = reflect_bisector(σ_{10}) = (d_1 · cos α, d_1 · sin α)
///   σ_{02} = reflect_bisector(σ_{20}) = (σ_{20,x} · cos α, σ_{20,x} · sin α)
/// ```
///
/// Reflection across a line through origin at angle `θ` maps
/// `(x, y) ↔ (x cos 2θ + y sin 2θ,  x sin 2θ − y cos 2θ)`. With
/// `θ = α/2` and `(x, y) = (d_1, 0)`, this gives `(d_1 cos α,
/// d_1 sin α)`.
///
/// The bi-4 partials at `(0, 0)` follow from the BB-net:
/// ```text
///   ∂_u σ(0, 0) = 4·(σ_{10} − σ_{00}) = 4·σ_{10}
///   ∂_v σ(0, 0) = 4·σ_{01}
///   ∂²_u σ(0, 0) = 12·(σ_{20} − 2·σ_{10} + σ_{00}) = 12·(σ_{20} − 2·σ_{10})
///   ∂²_v σ(0, 0) = 12·(σ_{02} − 2·σ_{01})
///   ∂²_{uv} σ(0, 0) = 16·(σ_{11} − σ_{10} − σ_{01} + σ_{00})
///                   = 16·(σ_{11} − σ_{10} − σ_{01})
/// ```
pub fn sigma_center_jet(n: usize, d1: f64, d2: f64) -> SigmaCenterJet {
    assert!(n >= 3, "PnS σ: n must be ≥ 3, got {n}");
    let alpha = 2.0 * std::f64::consts::PI / (n as f64);
    let c = alpha.cos();
    let (ha_sin, ha_cos) = (alpha * 0.5).sin_cos();

    // σ_{10} = (d_1, 0)
    let s10 = [d1, 0.0];
    // σ_{11} = (d_2 · cos(α/2), d_2 · sin(α/2))
    let s11 = [d2 * ha_cos, d2 * ha_sin];
    // σ_{20,x} = (1 − 4c/3) d_1 + (4 cos(α/2) / (3c)) d_2
    let s20_x = (1.0 - 4.0 * c / 3.0) * d1 + (4.0 * ha_cos / (3.0 * c)) * d2;
    let s20 = [s20_x, 0.0];
    let alpha_sin = alpha.sin();
    let s01 = [d1 * c, d1 * alpha_sin];
    let s02 = [s20_x * c, s20_x * alpha_sin];
    // σ_{00} = (0, 0)
    let s00 = [0.0, 0.0];

    let ds_u = [4.0 * s10[0], 4.0 * s10[1]];
    let ds_v = [4.0 * s01[0], 4.0 * s01[1]];
    let ds_uu = [
        12.0 * (s20[0] - 2.0 * s10[0]),
        12.0 * (s20[1] - 2.0 * s10[1]),
    ];
    let ds_vv = [
        12.0 * (s02[0] - 2.0 * s01[0]),
        12.0 * (s02[1] - 2.0 * s01[1]),
    ];
    let ds_uv = [
        16.0 * (s11[0] - s10[0] - s01[0]),
        16.0 * (s11[1] - s10[1] - s01[1]),
    ];

    SigmaCenterJet {
        s00,
        ds_u,
        ds_v,
        ds_uu,
        ds_vv,
        ds_uv,
    }
}

/// Build the σ center jet using paper-default scaling parameters
/// [`d1_def`], [`d2_def`]. This is the entry point for the bi-4
/// cap's G² construction at the pole.
pub fn sigma_center_jet_default(n: usize) -> SigmaCenterJet {
    sigma_center_jet(n, d1_def(n), d2_def(n))
}

#[cfg(test)]
#[path = "pns_sigma_tests.rs"]
mod tests;
