//! Chain-curve control-spacing diagnostics and decompression.
//!
//! The constrained LSQ chain fit leaves the curve's speed profile
//! free, and the inherited row parameterization compresses it at
//! chain ends. Compressed end legs shrink the frozen transverse
//! legs of crossing chains (CP0→CP1 distances), thinning the
//! derived cross-derivative bands until fixed-size twist effects
//! rotate their normals O(1). This module measures the end
//! spacing and produces soft retarget rows that restore
//! uniform-scale end legs along the SOLVED leg directions — a
//! speed-profile change only; tangent directions and node pins
//! are unaffected.

use popkern_core::geometry::Point3;

/// End spacing below this fraction of uniform is treated as fit
/// compression and retargeted (a degree-elevated quadratic sits
/// near 1.0, so legitimate geometry stays untouched). Compressed
/// end legs shrink the frozen transverse legs of crossing chains,
/// thinning the derived cross-derivative bands until the
/// fixed-size twist-compatibility gap rotates their normals O(1).
const END_SPACING_COMPRESSION_THRESHOLD: f64 = 0.6;

/// Soft weight for decompression CP targets — strong enough to
/// override the sample mass on one control point, far below the
/// hard-row scale.
const END_SPACING_TARGET_WEIGHT: f64 = 50.0;

/// Detect compressed end control spacing and produce soft retarget
/// rows restoring uniform-scale end legs along the SOLVED leg
/// directions (speed-profile change only; tangent directions and
/// node pins are unaffected).
pub(super) fn end_decompression_targets(
    chain_cps: &[Point3; 10],
) -> Option<Vec<(usize, Point3, f64)>> {
    let total: f64 = (1..10)
        .map(|i| chain_cps[i].distance_to(&chain_cps[i - 1]))
        .sum();
    if total <= 1.0e-300 {
        return None;
    }
    let uniform = total / 9.0;
    let mut targets = Vec::new();
    let first = chain_cps[1] - chain_cps[0];
    let first_length = first.dot(&first).sqrt();
    if first_length < END_SPACING_COMPRESSION_THRESHOLD * uniform && first_length > 1.0e-300 {
        let target = chain_cps[0] + first * (uniform / first_length);
        targets.push((1usize, target, END_SPACING_TARGET_WEIGHT));
    }
    let last = chain_cps[8] - chain_cps[9];
    let last_length = last.dot(&last).sqrt();
    if last_length < END_SPACING_COMPRESSION_THRESHOLD * uniform && last_length > 1.0e-300 {
        let target = chain_cps[9] + last * (uniform / last_length);
        targets.push((8usize, target, END_SPACING_TARGET_WEIGHT));
    }
    if targets.is_empty() {
        None
    } else {
        Some(targets)
    }
}

/// Smallest end control spacing relative to uniform: the lesser
/// of the first and last control-polygon leg lengths divided by
/// one ninth of the polygon's total length. Near 1.0 for uniform
/// speed; ≪ 1 when the fit compresses CPs at the ends (which
/// shrinks the transverse legs of crossing chains and thins the
/// derived cross-derivative bands).
pub(super) fn end_spacing_ratio(chain_cps: &[Point3; 10]) -> f64 {
    let total: f64 = (1..10)
        .map(|i| chain_cps[i].distance_to(&chain_cps[i - 1]))
        .sum();
    if total <= 1.0e-300 {
        return 0.0;
    }
    let uniform = total / 9.0;
    let first = chain_cps[1].distance_to(&chain_cps[0]);
    let last = chain_cps[9].distance_to(&chain_cps[8]);
    first.min(last) / uniform
}
