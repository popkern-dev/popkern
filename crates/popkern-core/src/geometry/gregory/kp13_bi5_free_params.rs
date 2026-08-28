//! Free-parameter representation of `g△` and the
//! builder that threads the KP13 §3 / KP09a §2.3 pipeline.
//!
//! The KP13 LSQ (eq 4) minimises the boundary-fit residual over the
//! **free BB-coefs** of `g△`. After all C² constraints have been
//! factored in, the free parameters per axis are:
//!
//! | slot          | count | provenance                                  |
//! |---------------|------:|---------------------------------------------|
//! | Q's α, β, γ, δ, ε | 5 | global; ζ pinned to Catmull-Clark limit     |
//! | `m_k` shared mid-edges | n | one per unique inter-sector boundary  |
//! | `g^i_{30}` (paper 230) | n | free depth-0 shared-edge coef, 1/sector |
//! | `g^i_{40}` (paper 140) | n | next depth-0 shared-edge coef, 1/sector |
//! | 122, 023, 032 interior | 3n | KP13 Fig 4's free interior labels, 1/sector each |
//!
//! Total: **6n + 5** scalar DoFs per output axis, matching KP13
//! §5.1's `(6n + 5) × (6n + 5)` linear-system size.
//!
//! Everything else in the `Gtri5` is derived from these via the
//! C² machinery built in I.11b–d + I.11.d.2:
//!
//! * Star (flat 0) = CC limit (pinned, not in the free set)
//! * Near-star 5 (flat 1, 2, 6, 7, 11) = Q identity
//! * `g^i_{03}` (flat 3) = sector `(i+1)`'s `g^{i+1}_{30}` (edge share)
//! * `g^i_{04}` (flat 4) = sector `(i+1)`'s `g^{i+1}_{40}` (edge share)
//! * `g^i_{05}` (flat 5) = `m_{i+1}` (mid-edge share)
//! * `g^i_{21}` (flat 12) = circulant solve (KP09a eq 8, Lemma 1)
//! * `g^i_{12}` (flat 8)  = eq (3) at j = 2 back-fill
//! * `g^i_{31}` (flat 16) = eq (3) at j = 3 (higher-j C¹)
//! * `g^i_{41}` (flat 19) = eq (3) at j = 4 (higher-j C¹)
//!
//! [`apply_free_params`] is the canonical entry point that runs the
//! whole chain and returns a fully-populated `Gtri5` ready for
//! downstream consumption (e.g., the LSQ design-matrix assembly or
//! the I.12 bi-5 guide construction).

#![allow(dead_code)]

use super::kp13_bi5_circulant::{CirculantError, complete_derivative_layers_with_higher_j};
use super::kp13_bi5_connect::paper_flat;
use super::kp13_bi5_gtri5::{Gtri5, IDX_023, IDX_032, IDX_122, IDX_230};
use super::kp13_bi5_star_quadratic::StarQuadratic;
use crate::geometry::Point3;

/// Free BB-coefs of `g△`. `6n + 5` scalar DoFs per axis.
///
/// Invariants enforced by [`apply_free_params`] (not by this struct
/// alone):
///
/// * `m_k.len() == g_30.len() == g_40.len() == g_122.len() ==
///   g_023.len() == g_032.len() == n`.
/// * Q's `ζ` field is ignored by [`apply_free_params`] and replaced
///   by the supplied Catmull-Clark limit.
#[derive(Debug, Clone)]
pub struct Bi5FreeParams {
    /// Star quadratic `Q(x, y)`. Only `alpha..epsilon` are consumed
    /// — `zeta` is overridden by the Catmull-Clark limit passed to
    /// [`apply_free_params`].
    pub q: StarQuadratic,
    /// Per-sector mid-edge vertex `m_k` (paper `050` of sector k =
    /// paper `005` of sector `k-1`). `len() == n`.
    pub m_k: Vec<Point3>,
    /// Per-sector `g^i_{30}` (paper `230`, flat 15). `len() == n`.
    pub g_30: Vec<Point3>,
    /// Per-sector `g^i_{40}` (paper `140`, flat 18). `len() == n`.
    pub g_40: Vec<Point3>,
    /// Per-sector interior `g^i_{22}` (paper `122`, flat 13).
    pub g_122: Vec<Point3>,
    /// Per-sector interior `g^i_{23}` (paper `023`, flat 14).
    pub g_023: Vec<Point3>,
    /// Per-sector interior `g^i_{32}` (paper `032`, flat 17).
    pub g_032: Vec<Point3>,
}

impl Bi5FreeParams {
    /// All-zero free parameters of the correct shape for valence `n`.
    /// Useful as an LSQ initial guess or baseline.
    pub fn zero(n: usize) -> Self {
        assert!(n >= 3, "Bi5FreeParams needs n ≥ 3, got {n}");
        let z = Point3::new(0.0, 0.0, 0.0);
        Self {
            q: StarQuadratic::zero(),
            m_k: vec![z; n],
            g_30: vec![z; n],
            g_40: vec![z; n],
            g_122: vec![z; n],
            g_023: vec![z; n],
            g_032: vec![z; n],
        }
    }

    /// Valence implied by the per-sector vector lengths.
    pub fn valence(&self) -> usize {
        self.m_k.len()
    }

    /// Structural sanity check.
    pub fn validate(&self) -> Result<(), String> {
        let n = self.m_k.len();
        if n < 3 {
            return Err(format!("Bi5FreeParams needs n ≥ 3, got {n}"));
        }
        for (field_name, vec_len) in [
            ("g_30", self.g_30.len()),
            ("g_40", self.g_40.len()),
            ("g_122", self.g_122.len()),
            ("g_023", self.g_023.len()),
            ("g_032", self.g_032.len()),
        ] {
            if vec_len != n {
                return Err(format!(
                    "Bi5FreeParams: {field_name}.len() = {vec_len}, expected {n}"
                ));
            }
        }
        Ok(())
    }
}

/// Error from [`apply_free_params`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ApplyError {
    /// Struct invariants violated (lengths, valence).
    Shape(String),
    /// Circulant solver rejected the valence.
    Circulant(CirculantError),
}

impl std::fmt::Display for ApplyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ApplyError::Shape(s) => write!(f, "apply_free_params shape: {s}"),
            ApplyError::Circulant(e) => write!(f, "apply_free_params circulant: {e}"),
        }
    }
}

impl std::error::Error for ApplyError {}

/// Build a fully-populated `Gtri5` from free parameters + CC star
/// limit.
///
/// Pipeline:
///
/// 1. Pin star to `cc_limit` (KP13 §3 star-point choice, §1.5
///    limit-point formula).
/// 2. Apply Q via [`StarQuadratic::apply_to_all_sectors_pinned`] —
///    this overrides `Q.ζ` with `cc_limit` and propagates the 5
///    remaining DoFs into every sector's near-star region.
/// 3. Write per-sector free coefs: `m_k` (via [`Gtri5::set_mid_edge`]),
///    `g^i_{30}`, `g^i_{40}`, and the three interior coefs.
/// 4. [`complete_derivative_layers_with_higher_j`] — syncs shared-
///    edge CPs at depths 3, 4; solves the circulant for
///    `g^i_{21}`; back-fills `g^i_{12}` (eq 3 at j = 2);
///    back-fills `g^i_{31}`, `g^i_{41}` (eq 3 at j = 3, 4).
///
/// Returns the populated `Gtri5` on success. Only rejects `n ∈ {3,
/// 4}` via [`ApplyError::Circulant`] — all other valences including
/// `n = 6` are solved (n = 6 via pseudo-inverse).
///
/// The resulting `Gtri5` passes [`Gtri5::validate`] (star sharing,
/// mid-edge invariant). Whether the resulting `g△` is globally C²
/// depends on the interpretation of KP09a §2.3 — see the audit
/// discussion in `docs/journal/polar-kp-math.md` §3.2a.
pub fn apply_free_params(params: &Bi5FreeParams, cc_limit: Point3) -> Result<Gtri5, ApplyError> {
    params.validate().map_err(ApplyError::Shape)?;
    let n = params.valence();

    // 1. New Gtri5 with star pinned.
    let mut g = Gtri5::new_with_star_point(n, cc_limit);

    // 2. Apply Q with pinned ζ.
    params.q.apply_to_all_sectors_pinned(&mut g, cc_limit);

    // 3. Per-sector free coefs.
    for (k, &p) in params.m_k.iter().enumerate() {
        g.set_mid_edge(k, p);
    }
    let idx_30 = paper_flat(3, 0);
    let idx_40 = paper_flat(4, 0);
    for k in 0..n {
        g.sector_mut(k).cps[idx_30] = params.g_30[k];
        g.sector_mut(k).cps[idx_40] = params.g_40[k];
        g.sector_mut(k).cps[IDX_122] = params.g_122[k];
        g.sector_mut(k).cps[IDX_023] = params.g_023[k];
        g.sector_mut(k).cps[IDX_032] = params.g_032[k];
    }

    // 4. Full C² derivation pipeline.
    complete_derivative_layers_with_higher_j(&mut g).map_err(ApplyError::Circulant)?;

    // Sanity-check the result before returning.
    debug_assert!(
        g.validate().is_ok(),
        "apply_free_params produced invalid Gtri5"
    );
    // Reference-touch IDX_230 so rust-analyzer sees the import as used;
    // our flat-index constants are available for callers.
    let _ = IDX_230;

    Ok(g)
}

#[cfg(test)]
#[path = "kp13_bi5_free_params_tests.rs"]
mod tests;
