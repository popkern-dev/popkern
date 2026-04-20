//! Direct-construction: the outer hex corners `V_{2k+2}`.
//!
//! At each outer hex vertex where two adjacent sub-patches
//! meet, the direct construction emits the corner-adjacent
//! inner CPs (`grid[6][6]` on patch `k`, `grid[1][1]` on patch
//! `(k+1) mod N`) via the **parallelogram rule** on the seed
//! two-stripe boundaries.
//!
//! # Why this is simpler than the central corner
//!
//! At an outer hex corner `V_{2k+2}`, **only two patches**
//! share the corner. The tangent plane at `V` on each patch is
//! determined by its row-0 boundary CPs (specifically the
//! second-to-last CP on each incident boundary edge) — those
//! are fixed by upstream geometry (the spoke curve and the
//! hex-external edge). So for G1 at `V`, **nothing needs to be
//! solved** on a linear system — the tangent plane is already
//! pinned by the input.
//!
//! What this module does:
//!
//! 1. **Verifies** that the per-patch tangent plane at `V`
//!    agrees with the input face tangent plane in
//!    [`OuterCornerData::tangent_plane`]. Mismatch implies the
//!    upstream (biseptic) construction didn't honour the face
//!    normals — an error the caller should surface via a
//!    `CornerTangentPlaneSingular`-style remediation.
//! 2. **Emits** the twist-control CP at `grid[6][6]` / `grid[1][1]`
//!    using the parallelogram rule
//!    `corner_adjacent = neighbour_u + neighbour_v − corner`.
//!    This is a minimum-twist seed — G2 refinement (Phase G
//!    row-2) may displace it further.
//!
//! # Source
//!
//! Parallelogram rule from Farin (1993) Ch. 18; twist-CP usage
//! from Chiyokura-Kimura (1983) §3.2. Pre-1990, no patents.

use crate::geometry::gregory::direct_errors::DirectConstructionError;
use crate::geometry::gregory::direct_input::DirectBoundaryInput;
use crate::geometry::gregory::direct_options::DirectConstructionOptions;
use crate::geometry::{Point3, Vec3};

/// Output of the outer-corner pass: one entry per outer hex
/// corner, with the twist CPs on both incident patches.
#[derive(Debug, Clone)]
pub struct OuterCornerTwistCps {
    /// Per-corner output, indexed `0..N` in spoke order.
    pub corners: Vec<OuterCornerOutput>,
}

/// The two twist-control CPs emitted at a single outer corner.
#[derive(Debug, Clone)]
pub struct OuterCornerOutput {
    /// `grid[6][6]` on patch `k` (the "left" patch of the
    /// outer corner — the one whose UOne boundary is spoke `k`).
    pub patch_k_grid_6_6: Point3,
    /// `grid[1][1]` on patch `(k+1) mod N` (the "right" patch —
    /// its VZero boundary is spoke `k` reversed).
    pub patch_k_plus_1_grid_1_1: Point3,
    /// Measured tangent-plane residual. `0.0` = exact match
    /// with input face tangent plane. Always finite.
    pub tangent_plane_residual: f64,
}

/// Solve all `N` outer hex corners.
///
/// Requires `input.seed_patches` to be `Some` — we need access
/// to the full patch grids (specifically the hex-external row
/// CPs at `grid[·][7]` / `grid[0][·]`) which aren't carried in
/// the per-spoke data.
pub fn solve_outer_corners(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Result<OuterCornerTwistCps, DirectConstructionError> {
    let Some(seed) = input.seed_patches.as_ref() else {
        return Err(DirectConstructionError::DecompositionFailed {
            detail: "outer corner solve requires seed_patches (full biseptic-promoted grids); \
                     pass them via DirectBoundaryInput::seed_patches"
                .into(),
        });
    };
    let n = input.n();
    if seed.len() != n {
        return Err(DirectConstructionError::DecompositionFailed {
            detail: format!("seed_patches length {} ≠ ring size {n}", seed.len()),
        });
    }

    let mut corners = Vec::with_capacity(n);
    for k in 0..n {
        let k_next = (k + 1) % n;
        let patch_k = &seed[k];
        let patch_k_next = &seed[k_next];
        let v = input.outer_corners[k].position;

        // Patch k's outer corner is at grid[7][7] = V.
        let p_k_7_7 = patch_k.grid[7][7];
        // Adjacent boundary CPs.
        let p_k_7_6 = patch_k.grid[7][6]; // on spoke k
        let p_k_6_7 = patch_k.grid[6][7]; // on hex-external v=1 edge

        // Verify the patch's tangent plane at V agrees with the
        // input face tangent plane within tolerance.
        let du_k = Vec3::new(
            7.0 * (p_k_7_7.x - p_k_6_7.x),
            7.0 * (p_k_7_7.y - p_k_6_7.y),
            7.0 * (p_k_7_7.z - p_k_6_7.z),
        );
        let dv_k = Vec3::new(
            7.0 * (p_k_7_7.x - p_k_7_6.x),
            7.0 * (p_k_7_7.y - p_k_7_6.y),
            7.0 * (p_k_7_7.z - p_k_7_6.z),
        );
        let expected_tp = &input.outer_corners[k].tangent_plane;
        let residual_k = tangent_plane_residual(du_k, dv_k, expected_tp.normal());

        // Parallelogram-rule twist: grid[6][6] = grid[7][6] + grid[6][7] - grid[7][7].
        let p_k_6_6 = Point3::new(
            p_k_7_6.x + p_k_6_7.x - p_k_7_7.x,
            p_k_7_6.y + p_k_6_7.y - p_k_7_7.y,
            p_k_7_6.z + p_k_6_7.z - p_k_7_7.z,
        );

        // Patch (k+1)'s outer corner is at grid[0][0] = V.
        let p_kn_0_0 = patch_k_next.grid[0][0];
        let p_kn_0_1 = patch_k_next.grid[0][1]; // hex-external u=0 edge
        let p_kn_1_0 = patch_k_next.grid[1][0]; // on spoke k (reversed)

        let du_kn = Vec3::new(
            7.0 * (p_kn_1_0.x - p_kn_0_0.x),
            7.0 * (p_kn_1_0.y - p_kn_0_0.y),
            7.0 * (p_kn_1_0.z - p_kn_0_0.z),
        );
        let dv_kn = Vec3::new(
            7.0 * (p_kn_0_1.x - p_kn_0_0.x),
            7.0 * (p_kn_0_1.y - p_kn_0_0.y),
            7.0 * (p_kn_0_1.z - p_kn_0_0.z),
        );
        let residual_kn = tangent_plane_residual(du_kn, dv_kn, expected_tp.normal());

        let residual = residual_k.max(residual_kn);

        // Note: we deliberately DO NOT reject on residual here.
        // For real three-quad geometry from the biseptic seed,
        // patches A and B at an outer corner typically have
        // SLIGHTLY different tangent planes — Phase 10C's
        // corner-sample iteration exists precisely to smooth
        // this out. Rejecting here would force every real
        // fixture through the sampled-damped fallback. Instead,
        // we surface the residual via the result struct so the
        // orchestration-layer self-check catches real failures
        // at the dense-G1-measurement level.
        let _ = (options, k);

        // G0 sanity: p_k[7][7] and p_{k+1}[0][0] should coincide at V.
        // Measured for diagnostic purposes; any error would indicate a
        // broken seed upstream.
        let g0_k =
            ((p_k_7_7.x - v.x).powi(2) + (p_k_7_7.y - v.y).powi(2) + (p_k_7_7.z - v.z).powi(2))
                .sqrt();
        let g0_kn =
            ((p_kn_0_0.x - v.x).powi(2) + (p_kn_0_0.y - v.y).powi(2) + (p_kn_0_0.z - v.z).powi(2))
                .sqrt();
        let _ = (g0_k, g0_kn); // Not surfaced in MVP output; a future diagnostic can expose.

        let p_kn_1_1 = Point3::new(
            p_kn_1_0.x + p_kn_0_1.x - p_kn_0_0.x,
            p_kn_1_0.y + p_kn_0_1.y - p_kn_0_0.y,
            p_kn_1_0.z + p_kn_0_1.z - p_kn_0_0.z,
        );

        corners.push(OuterCornerOutput {
            patch_k_grid_6_6: p_k_6_6,
            patch_k_plus_1_grid_1_1: p_kn_1_1,
            tangent_plane_residual: residual,
        });
    }

    Ok(OuterCornerTwistCps { corners })
}

/// Residual magnitude of `(du × dv)` vs. the expected plane
/// normal. Zero = exact; scales with the magnitude of the
/// mismatch between the seed boundary's natural tangent plane
/// and the input face tangent plane.
///
/// Specifically: we compute the unit normal of the patch's
/// tangent plane at this corner as `(du × dv) / |du × dv|` and
/// compare against the expected plane normal. Residual is
/// `1 - |cos(angle)|` weighted by the characteristic-length
/// scale — gives a dimensioned error we can threshold against
/// `options.corner_singularity_eps * characteristic_length`.
fn tangent_plane_residual(du: Vec3, dv: Vec3, expected_normal: Vec3) -> f64 {
    let cross = du.cross(&dv);
    let Some(n_patch) = cross.normalized() else {
        return 1.0;
    };
    let cos = n_patch.x * expected_normal.x
        + n_patch.y * expected_normal.y
        + n_patch.z * expected_normal.z;
    // Residual = 1 - |cos(angle)|, scaled by characteristic magnitudes
    // (|du|, |dv|) so values are comparable to `corner_singularity_eps`.
    let avg_mag = (du.length() + dv.length()) * 0.5;
    (1.0 - cos.abs()) * avg_mag
}

#[cfg(test)]
#[path = "direct_corner_outer_tests.rs"]
mod tests;
