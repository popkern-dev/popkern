//! Phase I — Shared polar K-P infrastructure.
//!
//! Utilities used by **both** K-P 2015 (bi-6) and K-P 2013 (bi-5)
//! polar constructions. Living here rather than inside the
//! per-paper modules prevents cyclic imports and keeps the math
//! spec honest — e.g., the Catmull-Clark limit-point formula
//! below is identical for both papers (math summary §1.5), so
//! there's one canonical implementation.
//!
//! Scope grows across the phase:
//! - I.7 — Catmull-Clark star point.
//! - I.9 / I.14 — F_κ quality-functional primitives (added
//!   later; left here as a place-holder module comment).
//! - I.17 — higher-level shared orchestration.

use crate::geometry::Point3;

/// Weighted sum `Σ λ_i · P_i` coordinate-wise (per the K-P
/// paper convention; see `kp15_bi6_rules` module docs for
/// discussion of non-affine linear combinations of BB-CPs).
#[inline]
pub(super) fn weighted_sum(pairs: &[(f64, Point3)]) -> Point3 {
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for &(lam, p) in pairs {
        out.x += lam * p.x;
        out.y += lam * p.y;
        out.z += lam * p.z;
    }
    out
}

/// Catmull-Clark limit point (star point) of a polar cap's hub
/// per [KP15 §3.6 step 1] and [KP13 §3] (math summary §1.5):
///
/// ```text
///   p_{00} = (1 − n·ν₅ − n·ν₆)·c₇  +  Σ_{k=0..n−1} (ν₅·c^k_5 + ν₆·c^k_6)
///   ν₅ := 1 / (n·(n+5))  +  (1/96  if n = 3  else 0)
///   ν₆ := 4 / (n·(n+5))
/// ```
///
/// `c₇` is the star point of the input cap-net; `c5s[k]` and
/// `c6s[k]` are the two outer-ring nodes on sector `k` (length
/// `n` each).
///
/// Panics if `c5s.len() != n` or `c6s.len() != n`.
pub fn catmull_clark_star_point(n: usize, c7: Point3, c5s: &[Point3], c6s: &[Point3]) -> Point3 {
    assert_eq!(c5s.len(), n, "catmull_clark: c5s length mismatch");
    assert_eq!(c6s.len(), n, "catmull_clark: c6s length mismatch");
    let nf = n as f64;
    let nu5 = 1.0 / (nf * (nf + 5.0)) + if n == 3 { 1.0 / 96.0 } else { 0.0 };
    let nu6 = 4.0 / (nf * (nf + 5.0));
    let center_weight = 1.0 - nf * nu5 - nf * nu6;
    let mut pairs = Vec::with_capacity(1 + 2 * n);
    pairs.push((center_weight, c7));
    for k in 0..n {
        pairs.push((nu5, c5s[k]));
        pairs.push((nu6, c6s[k]));
    }
    weighted_sum(&pairs)
}

#[cfg(test)]
#[path = "kp_polar_shared_tests.rs"]
mod tests;
