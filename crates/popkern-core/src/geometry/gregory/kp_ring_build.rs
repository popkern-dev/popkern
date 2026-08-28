//! Ring-wide bidegree-9 patch construction via the
//! K-P guided-surface pipeline.
//!
//! Combines H.4.c (guide ring) + H.4.d (row-1 frame) + H.4.e
//! (row-2 frame) into a single orchestrator that emits
//! [`RectangularGregoryBiNonic`] per sub-patch.
//!
//! # Shared-spoke CP averaging (G4)
//!
//! The independent-guide construction of H.4.c leaves small
//! residuals on adjacent guides' cross-tangent curves along
//! their shared spokes (measured on truly-3D input at ~0.34 in
//! triple-product units). To guarantee **exact G1 and G2 across
//! shared spokes** on the final bidegree-9 patches regardless of
//! guide-level residual, this module post-processes the derived
//! row-1 / row-2 frames to enforce the shared-CP identity:
//!
//! ```text
//!   row8[j] of sub-patch k  =  col1[9 − j] of sub-patch (k+1)
//! ```
//!
//! (and analogously for row-2: `row7[j] of k = col2[9 − j] of
//! (k+1)`). Implementation: average the two adjacent sub-patches'
//! predictions and write the averaged value back to both.
//!
//! After averaging, the bidegree-9 patch's cross-tangent /
//! cross-second-derivative at the shared spoke matches *exactly*
//! between adjacent sub-patches → G1 and G2 **exact** along the
//! full shared spoke. The residual that H.4.c exposed is
//! absorbed into a small guide-level approximation (the guide's
//! cross-tangent field is slightly redistributed at the shared-
//! spoke interior) which does not affect the continuity between
//! sub-patches.
//!
//! # Boundary construction
//!
//! Input boundaries for each sub-patch are built K-P-consistent
//! with the guide — tangent and curvature-consistent at every
//! corner — so H.4.d / H.4.e's corner 2×2 / 3×3 formulas match
//! the guide to floating-point precision.

use super::direct_errors::DirectConstructionError;
use super::direct_input::DirectBoundaryInput;
use super::direct_options::DirectConstructionOptions;
use super::kp_degree_elevate::elevate_point_curve_to_nine;
use super::kp_guide_patch::KpGuidePatch;
use super::kp_guide_patch_binonic::KpGuidePatchBiNonic;
use super::kp_guide_patch_biquintic::KpGuidePatchBiQuintic;
use super::kp_guide_ring::build_guide_ring_with_twist;
use super::kp_interior_fit::{
    fit_interior_cps_to_binonic_guide, fit_interior_cps_to_biquintic_guide,
    fit_interior_cps_to_guide,
};
use super::kp_row1::{Row1Frame, derive_row1_frame, derive_row1_frame_from_boundary};
use super::kp_row2::{Row2Frame, derive_row2_frame, derive_row2_frame_from_boundary};
use super::rectangular_binonic::RectangularGregoryBiNonic;
use crate::geometry::Point3;
#[cfg(test)]
use crate::geometry::Vec3;

#[path = "kp_ring_build_helpers.rs"]
mod kp_ring_build_helpers;

use kp_ring_build_helpers::{
    assemble_patch, fill_interior_cps_via_coons, pick_hex_boundary, reverse_boundary,
};

/// Build the full N-sub-patch bidegree-9 ring for a vertex blend.
///
/// Pipeline: guide ring (with twist correction) → per-sub-patch
/// boundaries (K-P consistent) → H.4.d row-1 frame → H.4.e
/// row-2 frame → shared-spoke CP averaging → interior CPs left
/// at zero (H.5 territory).
///
/// Returns one [`RectangularGregoryBiNonic`] per sub-patch.
/// Interior CPs (rows 3..=6, cols 3..=6) are zero — the caller
/// (or H.5) fills them via least-squares fit to the guide.
pub fn build_binonic_ring(
    input: &DirectBoundaryInput,
    options: &DirectConstructionOptions,
) -> Result<Vec<RectangularGregoryBiNonic>, DirectConstructionError> {
    let guides = build_guide_ring_with_twist(input, options)?;
    build_binonic_ring_from_guides(&guides, input)
}

/// Build a bidegree-9 ring from a pre-built set of guides plus the
/// input boundary data. Exposed so the iterative K-P guide-twist
/// optimizer (the bi-9 commercial gates gap #2) can rebuild the ring after each
/// twist update without re-running the initial Shirman-Séquin
/// guide-ring construction.
pub fn build_binonic_ring_from_guides(
    guides: &[KpGuidePatch],
    input: &DirectBoundaryInput,
) -> Result<Vec<RectangularGregoryBiNonic>, DirectConstructionError> {
    let n = guides.len();
    if input.n() != n {
        return Err(DirectConstructionError::DecompositionFailed {
            detail: format!("guide count {} ≠ input.n() {}", n, input.n()),
        });
    }

    // Build shared-spoke boundaries ONCE per spoke from the
    // input spoke data (degree-elevated from 7 to 9). This
    // guarantees exact G0 on shared spokes regardless of per-
    // guide twist differences.
    let mut shared_spokes: Vec<[Point3; 10]> = Vec::with_capacity(input.n());
    for s in &input.spokes {
        let Some(elevated) = elevate_point_curve_to_nine(&s.boundary_cps) else {
            return Err(DirectConstructionError::DecompositionFailed {
                detail: format!("spoke {} degree-elevation from 7 to 9 failed", s.edge_id),
            });
        };
        shared_spokes.push(elevated);
    }

    // Per-sub-patch boundaries: shared spokes for row9 + col0,
    // guide-derived hex edges for row0 + col9.
    //
    //   row0 (u = 0, v varies) — left hex edge, guide-derived.
    //   row9 (u = 1, v varies) — spoke k, input-derived.
    //   col0 (v = 0, u varies) — spoke (k - 1) reversed, input-derived.
    //   col9 (v = 1, u varies) — top hex edge, guide-derived.
    let mut boundaries: Vec<[[Point3; 10]; 4]> = Vec::with_capacity(n);
    for k in 0..n {
        let k_prev = (k + n - 1) % n;
        let (row0, col9) = pick_hex_boundary(input, &guides[k], k);
        let row9 = shared_spokes[k];
        let col0 = reverse_boundary(&shared_spokes[k_prev]);
        boundaries.push([row0, row9, col0, col9]);
    }

    // the sphere-fit edge-wiring work Stage 2 dispatch: when boundary injection is
    // active, switch row-1 / row-2 derivation to a boundary-driven
    // path that bypasses the K-P guide's pre-injection partials.
    // The guide's `partial_v(0, 0)` reflects the bicubic fit of the
    // *old* boundary; under injection, the new boundary's tangent
    // 9·(row0[1] − row0[0]) replaces it. Same logic for the other
    // 7 corner tangents and (via linear blend) the cross-derivative
    // profile. Twist `G_uv` is set to 0 (Stage 2 simplification).
    //
    // Empirical caveat (matrix probe in
    // `kp_binonic_boundary_injection.rs`): on cube_mixed and
    // pyramid4_mixed under sphere-fit injection, the resulting
    // bidegree-9 patches still self-intersect — the guide-fitted
    // interior 16 CPs (rows 3..=6, cols 3..=6) are pulled toward
    // a guide that no longer matches the boundary-consistent
    // frame. Stage 3 needs to either rebuild the guide to match
    // the injection OR move row3..=6 / col3..=6 to a boundary-
    // consistent fill (Coons-style or L-M-optimised with row1/2
    // as additional DoF). This dispatch is the foundation for
    // that next stage; today it produces frames that the existing
    // self-intersection gate rejects, falling back to two-stripe
    // — same end-user behaviour as Stage 1.
    let inject = input.hex_boundaries.is_some();

    // H.4.d: row-1 frames.
    let mut r1_frames: Vec<Row1Frame> = if inject {
        boundaries
            .iter()
            .map(|b| derive_row1_frame_from_boundary(&b[0], &b[1], &b[2], &b[3]))
            .collect()
    } else {
        guides
            .iter()
            .zip(boundaries.iter())
            .map(|(g, b)| derive_row1_frame(g, &b[0], &b[1], &b[2], &b[3]))
            .collect()
    };

    // H.4.e: row-2 frames.
    let mut r2_frames: Vec<Row2Frame> = if inject {
        boundaries
            .iter()
            .zip(r1_frames.iter())
            .map(|(b, r1)| derive_row2_frame_from_boundary(&b[0], &b[1], &b[2], &b[3], r1))
            .collect()
    } else {
        guides
            .iter()
            .zip(boundaries.iter())
            .zip(r1_frames.iter())
            .map(|((g, b), r1)| derive_row2_frame(g, &b[0], &b[1], &b[2], &b[3], r1))
            .collect()
    };

    // G4: shared-spoke CP averaging. For each shared spoke k
    // (between sub-patch k and sub-patch (k+1) % N), average
    // row8 of k with col1 of (k+1) (reversed), and analogously
    // for row7 / col2. Guarantees pairwise CP identity on
    // shared spokes → exact G1, G2 across shared spokes on the
    // bidegree-9 patches.
    // Averaging is restricted to j ∈ 1..=8 — the interior of
    // the shared-spoke-adjacent rows/cols. j = 0 (near the
    // C corner) and j = 9 (near the V corner) touch CPs that
    // are on DIFFERENT boundaries (col 0 / col 9 for row8[0] /
    // row8[9]), so averaging them would corrupt the boundary
    // invariants. The corner 2×2 / 3×3 blocks are handled by
    // the H.4.d / H.4.e corner formulas separately.
    for k in 0..n {
        let k_next = (k + 1) % n;
        for j in 1..=8 {
            let a = r1_frames[k].row8[j];
            let b = r1_frames[k_next].col1[9 - j];
            let avg = Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
            r1_frames[k].row8[j] = avg;
            r1_frames[k_next].col1[9 - j] = avg;
        }
        for j in 1..=8 {
            let a = r2_frames[k].row7[j];
            let b = r2_frames[k_next].col2[9 - j];
            let avg = Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
            r2_frames[k].row7[j] = avg;
            r2_frames[k_next].col2[9 - j] = avg;
        }
    }

    // Assemble bidegree-9 patches with frame populated, then
    // fit interior via H.5 least-squares projection to guide.
    //
    // When boundary injection is active, the
    // pre-injection guide doesn't match the boundary-derived
    // frame, so its LSQ fit pulls the interior 16 CPs toward
    // a guide that no longer represents the surface — patches
    // fold. Substitute a transfinite-interpolation fill that
    // smoothly blends the boundary-consistent ring frame.
    let mut patches: Vec<RectangularGregoryBiNonic> = Vec::with_capacity(n);
    for k in 0..n {
        let mut patch = assemble_patch(
            &boundaries[k][0],
            &boundaries[k][1],
            &boundaries[k][2],
            &boundaries[k][3],
            &r1_frames[k],
            &r2_frames[k],
        );
        // Under injection, the K-P guide's
        // boundaries reflect the *pre-injection* deg-5 fit, which
        // disagrees with the post-injection deg-9 boundary by O(mm)
        // on curved arcs. The LSQ fit-to-guide pulls interior CPs
        // toward this stale reference, fighting the boundary-
        // consistent row-1/row-2 frame and producing folded patches.
        // Substitute a Coons transfinite interpolation that fills
        // the 16 strict-interior CPs from the boundary-consistent
        // ring directly — no guide reference, no stale-boundary
        // pull, fold-free by construction.
        if inject {
            fill_interior_cps_via_coons(&mut patch);
        } else {
            fit_interior_cps_to_guide(&mut patch, &guides[k]);
        }
        patches.push(patch);
    }
    Ok(patches)
}

/// Option C — build a bidegree-9 ring from a pre-built bi-quintic
/// guide ring. Same downstream pipeline as the bicubic path, but
/// pulls boundaries from the bi-quintic's higher-degree cross-
/// derivatives and interior-fits against the 36-CP guide via the
/// degree-5-vs-9 cross-Gram.
pub fn build_binonic_ring_from_biquintic_guides(
    guides: &[KpGuidePatchBiQuintic],
    input: &DirectBoundaryInput,
) -> Result<Vec<RectangularGregoryBiNonic>, DirectConstructionError> {
    let n = guides.len();
    if input.n() != n {
        return Err(DirectConstructionError::DecompositionFailed {
            detail: format!("guide count {} ≠ input.n() {}", n, input.n()),
        });
    }
    let mut shared_spokes: Vec<[Point3; 10]> = Vec::with_capacity(input.n());
    for s in &input.spokes {
        let Some(elevated) = elevate_point_curve_to_nine(&s.boundary_cps) else {
            return Err(DirectConstructionError::DecompositionFailed {
                detail: format!("spoke {} degree-elevation from 7 to 9 failed", s.edge_id),
            });
        };
        shared_spokes.push(elevated);
    }
    let mut boundaries: Vec<[[Point3; 10]; 4]> = Vec::with_capacity(n);
    for k in 0..n {
        let k_prev = (k + n - 1) % n;
        let (row0, col9) = pick_hex_boundary(input, &guides[k], k);
        let row9 = shared_spokes[k];
        let col0 = reverse_boundary(&shared_spokes[k_prev]);
        boundaries.push([row0, row9, col0, col9]);
    }
    // the sphere-fit edge-wiring work Stage 2 dispatch — see `build_binonic_ring_from_guides`.
    let inject = input.hex_boundaries.is_some();
    let mut r1_frames: Vec<Row1Frame> = if inject {
        boundaries
            .iter()
            .map(|b| derive_row1_frame_from_boundary(&b[0], &b[1], &b[2], &b[3]))
            .collect()
    } else {
        guides
            .iter()
            .zip(boundaries.iter())
            .map(|(g, b)| derive_row1_frame(g, &b[0], &b[1], &b[2], &b[3]))
            .collect()
    };
    let mut r2_frames: Vec<Row2Frame> = if inject {
        boundaries
            .iter()
            .zip(r1_frames.iter())
            .map(|(b, r1)| derive_row2_frame_from_boundary(&b[0], &b[1], &b[2], &b[3], r1))
            .collect()
    } else {
        guides
            .iter()
            .zip(boundaries.iter())
            .zip(r1_frames.iter())
            .map(|((g, b), r1)| derive_row2_frame(g, &b[0], &b[1], &b[2], &b[3], r1))
            .collect()
    };

    // G4 shared-spoke CP averaging — same as the bicubic path.
    for k in 0..n {
        let k_next = (k + 1) % n;
        for j in 1..=8 {
            let a = r1_frames[k].row8[j];
            let b = r1_frames[k_next].col1[9 - j];
            let avg = Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
            r1_frames[k].row8[j] = avg;
            r1_frames[k_next].col1[9 - j] = avg;
        }
        for j in 1..=8 {
            let a = r2_frames[k].row7[j];
            let b = r2_frames[k_next].col2[9 - j];
            let avg = Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
            r2_frames[k].row7[j] = avg;
            r2_frames[k_next].col2[9 - j] = avg;
        }
    }

    let mut patches: Vec<RectangularGregoryBiNonic> = Vec::with_capacity(n);
    for k in 0..n {
        let mut patch = assemble_patch(
            &boundaries[k][0],
            &boundaries[k][1],
            &boundaries[k][2],
            &boundaries[k][3],
            &r1_frames[k],
            &r2_frames[k],
        );
        // See `build_binonic_ring_from_guides`.
        if inject {
            fill_interior_cps_via_coons(&mut patch);
        } else {
            fit_interior_cps_to_biquintic_guide(&mut patch, &guides[k]);
        }
        patches.push(patch);
    }
    Ok(patches)
}

/// Option C Axis-1 ceiling — build a bidegree-9 ring from a
/// pre-built bi-nonic guide ring. Same pipeline as the bi-quintic
/// path but the row-1/row-2 formulas and interior-fit collapse to
/// near-identity operations (the guide IS the output at matching
/// degree). Still runs through the same ring builder so G4 shared-
/// spoke averaging and boundary consistency invariants are
/// preserved.
pub fn build_binonic_ring_from_binonic_guides(
    guides: &[KpGuidePatchBiNonic],
    input: &DirectBoundaryInput,
) -> Result<Vec<RectangularGregoryBiNonic>, DirectConstructionError> {
    let n = guides.len();
    if input.n() != n {
        return Err(DirectConstructionError::DecompositionFailed {
            detail: format!("guide count {} ≠ input.n() {}", n, input.n()),
        });
    }
    let mut shared_spokes: Vec<[Point3; 10]> = Vec::with_capacity(input.n());
    for s in &input.spokes {
        let Some(elevated) = elevate_point_curve_to_nine(&s.boundary_cps) else {
            return Err(DirectConstructionError::DecompositionFailed {
                detail: format!("spoke {} degree-elevation from 7 to 9 failed", s.edge_id),
            });
        };
        shared_spokes.push(elevated);
    }
    let mut boundaries: Vec<[[Point3; 10]; 4]> = Vec::with_capacity(n);
    for k in 0..n {
        let k_prev = (k + n - 1) % n;
        let (row0, col9) = pick_hex_boundary(input, &guides[k], k);
        let row9 = shared_spokes[k];
        let col0 = reverse_boundary(&shared_spokes[k_prev]);
        boundaries.push([row0, row9, col0, col9]);
    }
    // the sphere-fit edge-wiring work Stage 2 dispatch — see `build_binonic_ring_from_guides`.
    let inject = input.hex_boundaries.is_some();
    let mut r1_frames: Vec<Row1Frame> = if inject {
        boundaries
            .iter()
            .map(|b| derive_row1_frame_from_boundary(&b[0], &b[1], &b[2], &b[3]))
            .collect()
    } else {
        guides
            .iter()
            .zip(boundaries.iter())
            .map(|(g, b)| derive_row1_frame(g, &b[0], &b[1], &b[2], &b[3]))
            .collect()
    };
    let mut r2_frames: Vec<Row2Frame> = if inject {
        boundaries
            .iter()
            .zip(r1_frames.iter())
            .map(|(b, r1)| derive_row2_frame_from_boundary(&b[0], &b[1], &b[2], &b[3], r1))
            .collect()
    } else {
        guides
            .iter()
            .zip(boundaries.iter())
            .zip(r1_frames.iter())
            .map(|((g, b), r1)| derive_row2_frame(g, &b[0], &b[1], &b[2], &b[3], r1))
            .collect()
    };
    for k in 0..n {
        let k_next = (k + 1) % n;
        for j in 1..=8 {
            let a = r1_frames[k].row8[j];
            let b = r1_frames[k_next].col1[9 - j];
            let avg = Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
            r1_frames[k].row8[j] = avg;
            r1_frames[k_next].col1[9 - j] = avg;
        }
        for j in 1..=8 {
            let a = r2_frames[k].row7[j];
            let b = r2_frames[k_next].col2[9 - j];
            let avg = Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
            r2_frames[k].row7[j] = avg;
            r2_frames[k_next].col2[9 - j] = avg;
        }
    }
    let mut patches: Vec<RectangularGregoryBiNonic> = Vec::with_capacity(n);
    for k in 0..n {
        let mut patch = assemble_patch(
            &boundaries[k][0],
            &boundaries[k][1],
            &boundaries[k][2],
            &boundaries[k][3],
            &r1_frames[k],
            &r2_frames[k],
        );
        // See `build_binonic_ring_from_guides`.
        if inject {
            fill_interior_cps_via_coons(&mut patch);
        } else {
            fit_interior_cps_to_binonic_guide(&mut patch, &guides[k]);
        }
        patches.push(patch);
    }
    Ok(patches)
}

#[cfg(test)]
#[path = "kp_ring_build_tests.rs"]
mod tests;
