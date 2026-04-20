//! Phase H — ring-wide guide-surface assembly.
//!
//! Composes the piecewise Hermite extractors from H.4.c.2–4 into
//! one [`KpGuidePatch`] per sub-patch:
//!
//! - [`extract_central_vertex_hermite`] supplies `c1 = (u = 1,
//!   v = 0) = C` for every sub-patch (all N sub-patches share
//!   the common tangent plane there).
//! - [`extract_hex_corners`] supplies `c0 = V_{2k}`,
//!   `c2 = V_{2k + 2}`, `c3 = M_k` for every sub-patch.
//!
//! Per sub-patch we assemble these into a `[GuideCornerData; 4]`
//! and delegate to [`build_guide_from_hermite`], which turns the
//! 4 × 4 Hermite data into the 16 CPs of a bicubic Bezier
//! [`KpGuidePatch`].
//!
//! # G0 on shared spokes — automatic
//!
//! Adjacent sub-patches `k` and `(k + 1) % N` share spoke `k`.
//! In the guide geometry, this is:
//!
//! * sub-patch `k`'s `u = 1` edge (`C → V_{2k + 2}`, CPs
//!   `grid[3][0..=3]`)
//! * sub-patch `(k + 1)`'s `v = 0` edge (`V_{2k + 2} → C`
//!   reversed, CPs `grid[0..=3][0]`)
//!
//! After assembly, the 4 CPs on each side agree by reversed
//! index — sub-patch `k`'s `grid[3][j]` equals sub-patch
//! `(k + 1)`'s `grid[3 − j][0]`. This is a direct consequence of
//!
//! - `tangent_v` at `c1` of guide `k` = `+ spoke_tangent_at_C`
//! - `tangent_u` at `c1` of guide `k + 1` = `− spoke_tangent_at_C`
//!   (the **exact negative** invariant preserved by
//!   [`extract_central_vertex_hermite`]), and
//! - `tangent_v` at `c2` of guide `k` = `+ spoke_tangent_at_V`
//! - `tangent_u` at `c0` of guide `k + 1` = `− spoke_tangent_at_V`
//!   (the [`extract_outer_corner_hermite`] negation invariant).
//!
//! So shared-spoke G0 holds by construction — no further solve.

use super::kp_guide_build::{GuideCornerData, build_guide_from_hermite};
use super::kp_guide_central_hermite::extract_central_vertex_hermite;
use super::kp_guide_central_twist::extract_central_vertex_hermite_with_twist;
use super::kp_guide_hex_hermite::extract_hex_corners;
use super::kp_guide_outer_twist::extract_hex_corners_with_twist;
use super::kp_guide_patch::KpGuidePatch;
use crate::geometry::gregory::direct_errors::DirectConstructionError;
use crate::geometry::gregory::direct_input::DirectBoundaryInput;
use crate::geometry::gregory::direct_options::DirectConstructionOptions;

/// Build a bicubic guide patch per sub-patch in the vertex-blend
/// ring. Output length equals `input.n()`.
///
/// Corner layout of each returned guide (recap from ADR §2.2):
///
/// ```text
///   grid[0][0] = V_{2k}      grid[0][3] = M_k
///   grid[3][0] = C           grid[3][3] = V_{2k + 2}
/// ```
///
/// (with `c0`, `c1`, `c2`, `c3` at counter-clockwise corners in
/// that order).
///
/// Errors are surfaced from [`extract_central_vertex_hermite`] —
/// specifically [`DirectConstructionError::CornerTangentPlaneSingular`]
/// when the input face tangent planes at `C` are too inconsistent
/// to admit a common normal.
pub fn build_guide_ring(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Result<Vec<KpGuidePatch>, DirectConstructionError> {
    let central = extract_central_vertex_hermite(input, options)?;
    let hex = extract_hex_corners(input);

    let n = input.n();
    debug_assert_eq!(central.len(), n);
    debug_assert_eq!(hex.len(), n);

    let mut guides = Vec::with_capacity(n);
    for k in 0..n {
        let corners: [GuideCornerData; 4] = [hex[k].c0, central[k], hex[k].c2, hex[k].c3];
        guides.push(build_guide_from_hermite(corners));
    }
    Ok(guides)
}

/// Build the ring with **full twist correction** — Shirman-Séquin
/// at `C` + pairwise K-P-style at every outer hex vertex. Use
/// this in production; [`build_guide_ring`] (twist = 0) is kept
/// as a baseline reference and test harness.
///
/// The twist-corrected version produces bicubic guides that are
/// G1-compatible with their neighbours **along the entire shared
/// spoke**, not merely at the endpoints. This is the property
/// H.4.d relies on for row-1 CP agreement between adjacent
/// bidegree-9 sub-patches.
pub fn build_guide_ring_with_twist(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Result<Vec<KpGuidePatch>, DirectConstructionError> {
    let central = extract_central_vertex_hermite_with_twist(input, options)?;
    let hex = extract_hex_corners_with_twist(input);

    let n = input.n();
    debug_assert_eq!(central.len(), n);
    debug_assert_eq!(hex.len(), n);

    let mut guides = Vec::with_capacity(n);
    for k in 0..n {
        let corners: [GuideCornerData; 4] = [hex[k].c0, central[k], hex[k].c2, hex[k].c3];
        guides.push(build_guide_from_hermite(corners));
    }
    Ok(guides)
}

/// Build the guide ring with a caller-supplied central-vertex
/// Hermite per sub-patch — the iterative K-P 2004 §4 entry point.
///
/// Unlike [`build_guide_ring_with_twist`], this skips the
/// Shirman-Séquin solve and uses the provided `central` array
/// directly. The caller is responsible for ensuring `central[k]`
/// has the right position / tangent_u / tangent_v / twist for
/// sub-patch `k` at `C`. Used by the gap-#2 iterative twist
/// optimizer in `popkern-fillet` — after the first call,
/// subsequent iterations only tweak `central[k].twist` and re-
/// invoke this builder.
pub fn build_guide_ring_from_centrals(
    input: &DirectBoundaryInput,
    central: &[GuideCornerData],
) -> Vec<KpGuidePatch> {
    let hex = extract_hex_corners_with_twist(input);
    let n = input.n();
    debug_assert_eq!(central.len(), n);
    debug_assert_eq!(hex.len(), n);
    let mut guides = Vec::with_capacity(n);
    for k in 0..n {
        let corners: [GuideCornerData; 4] = [hex[k].c0, central[k], hex[k].c2, hex[k].c3];
        guides.push(build_guide_from_hermite(corners));
    }
    guides
}

#[cfg(test)]
#[path = "kp_guide_ring_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "kp_guide_ring_gap_tests.rs"]
mod gap_tests;

#[cfg(test)]
#[path = "kp_guide_ring_twist_tests.rs"]
mod twist_tests;
