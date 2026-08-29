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

use popkern_core::geometry::{Point3, Vec3};

use popkern_core::topology::VertexHandle;

use super::chain_solver::dbernstein9;

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

/// Interior node speed below this fraction of the control
/// polygon's total length (= uniform speed on `[0, 1]`) is
/// treated as a compressed cut and retargeted. Single soft pass,
/// interior stations only. Threshold and weight are the
/// deterministic-sweep optimum (0.85/15: G1 1.36e-2; 0.7/5 gave
/// 0.41, 0.9+ or weight 20 regress — over-retargeting damages
/// corner zones; all-station/fixed-point/local-target variants
/// destabilized neighboring chains).
const NODE_SPEED_COMPRESSION_THRESHOLD: f64 = 0.85;

/// Soft weight for derivative retarget rows. Derivative basis
/// values reach ~9, so the effective row mass is ~80x the
/// per-sample position mass at weight 1.
const NODE_SPEED_TARGET_WEIGHT: f64 = 15.0;

/// Measure the solved curve's speed at every INTERIOR node
/// parameter and produce direction-preserving derivative
/// retarget rows where the speed is compressed below uniform.
/// The restriction rows that end at a T-junction cut inherit
/// `C'(cut)` as their frozen corner legs; a compressed cut thins
/// those legs until crossing bands turn tangential and trim
/// normals are noise (observed: tangential fraction 0.988 at a
/// healthy 119-degree corner, normal disagreement 0.295).
pub(super) fn node_speed_decompression_targets(
    chain_cps: &[Point3; 10],
    nodes_local: &[(VertexHandle, f64)],
) -> Vec<(f64, Vec3, f64)> {
    let total: f64 = (1..10)
        .map(|i| chain_cps[i].distance_to(&chain_cps[i - 1]))
        .sum();
    if total <= 1.0e-300 {
        return Vec::new();
    }
    let mut targets = Vec::new();
    for &(_, param) in nodes_local {
        if param <= 1.0e-6 || param >= 1.0 - 1.0e-6 {
            continue;
        }
        let dbasis = dbernstein9(param);
        let mut derivative = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..10 {
            derivative =
                derivative + Vec3::new(chain_cps[i].x, chain_cps[i].y, chain_cps[i].z) * dbasis[i];
        }
        let speed = derivative.dot(&derivative).sqrt();
        if speed >= NODE_SPEED_COMPRESSION_THRESHOLD * total || speed <= 1.0e-300 {
            continue;
        }
        let target = derivative * (total / speed);
        targets.push((param, target, NODE_SPEED_TARGET_WEIGHT));
    }
    targets
}
