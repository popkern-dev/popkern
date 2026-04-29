//! Phase H — ring-wide bidegree-9 patch construction via the
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
use super::kp_guide_trait::KpGuide;
use super::kp_interior_fit::{
    fit_interior_cps_to_binonic_guide, fit_interior_cps_to_biquintic_guide,
    fit_interior_cps_to_guide,
};
use super::kp_row1::{Row1Frame, derive_row1_frame, derive_row1_frame_from_boundary};
use super::kp_row2::{Row2Frame, derive_row2_frame, derive_row2_frame_from_boundary};
use super::rectangular_binonic::RectangularGregoryBiNonic;
use crate::geometry::{Point3, Vec3};

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
/// optimizer (Phase H.7 gap #2) can rebuild the ring after each
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

    // Phase 6c-D Stage 2 dispatch: when boundary injection is
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
    // Phase 6c-D Stage 2: when boundary injection is active, the
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
        // Phase 6c-D Stage 3a: under injection, the K-P guide's
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
    // Phase 6c-D Stage 2 dispatch — see `build_binonic_ring_from_guides`.
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
        // Phase 6c-D Stage 3a — see `build_binonic_ring_from_guides`.
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
    // Phase 6c-D Stage 2 dispatch — see `build_binonic_ring_from_guides`.
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
        // Phase 6c-D Stage 3a — see `build_binonic_ring_from_guides`.
        if inject {
            fill_interior_cps_via_coons(&mut patch);
        } else {
            fit_interior_cps_to_binonic_guide(&mut patch, &guides[k]);
        }
        patches.push(patch);
    }
    Ok(patches)
}

/// Phase 6c-D Stage 3 scaffolding — transfinite-interpolation
/// (Coons-style) fill for the bidegree-9 patch's strict-interior
/// 4×4 block (`grid[3..=6][3..=6]`).
///
/// With injection, the K-P guide doesn't represent the surface
/// the patch should approximate (boundary moved; guide didn't),
/// so the LSQ interior fit pulls CPs toward an inconsistent
/// reference and patches fold. This function instead fills the
/// 16 interior CPs by linearly blending the surrounding ring
/// frame (rings at index 2 and 7 along each axis), Coons style:
///
/// ```text
///   P[i][j] = (1-α)·P[2][j] + α·P[7][j]
///           + (1-β)·P[i][2] + β·P[i][7]
///           − (1-α)(1-β)·P[2][2] − α(1-β)·P[7][2]
///           − (1-α)β·P[2][7] − αβ·P[7][7]
/// ```
///
/// where `α = (i-2)/5`, `β = (j-2)/5`. Standard transfinite
/// interpolation on the inner [2..=7]² rectangle. Guarantees
/// fold-freeness: the 4×4 interior is a smooth bilinear blend
/// of the boundary-consistent ring.
///
/// Trade-off: the surface in the deep interior loses any "pull"
/// from the original cap-shape guide; it becomes a pure
/// frame-driven shape. Acceptable for cap surfaces because the
/// cap's outer-boundary fidelity dominates the user-visible
/// quality bar; deep interior is not reachable by the cap-vs-
/// neighbour edge measurements.
///
/// **Stage 3a status:** Default interior-fill under injection.
/// Eliminates the LSQ-against-stale-guide pull and produces
/// patches that the existing `self_intersects_with_tolerance(1e-6)`
/// gate accepts on pyramid4_mixed (4 binonic ship). cube_mixed
/// remains rejected by the quality regression guard
/// (`worst_ring_coplanarity` ratio > 1.2× pre-rung); Stage 3b will
/// address that.
fn fill_interior_cps_via_coons(patch: &mut RectangularGregoryBiNonic) {
    for i in 3..=6 {
        let alpha = ((i - 2) as f64) / 5.0;
        for j in 3..=6 {
            let beta = ((j - 2) as f64) / 5.0;
            let p2j = patch.grid[2][j];
            let p7j = patch.grid[7][j];
            let pi2 = patch.grid[i][2];
            let pi7 = patch.grid[i][7];
            let p22 = patch.grid[2][2];
            let p72 = patch.grid[7][2];
            let p27 = patch.grid[2][7];
            let p77 = patch.grid[7][7];
            let mix = |c2j: f64,
                       c7j: f64,
                       ci2: f64,
                       ci7: f64,
                       c22: f64,
                       c72: f64,
                       c27: f64,
                       c77: f64|
             -> f64 {
                (1.0 - alpha) * c2j + alpha * c7j + (1.0 - beta) * ci2 + beta * ci7
                    - (1.0 - alpha) * (1.0 - beta) * c22
                    - alpha * (1.0 - beta) * c72
                    - (1.0 - alpha) * beta * c27
                    - alpha * beta * c77
            };
            patch.grid[i][j] = Point3::new(
                mix(p2j.x, p7j.x, pi2.x, pi7.x, p22.x, p72.x, p27.x, p77.x),
                mix(p2j.y, p7j.y, pi2.y, pi7.y, p22.y, p72.y, p27.y, p77.y),
                mix(p2j.z, p7j.z, pi2.z, pi7.z, p22.z, p72.z, p27.z, p77.z),
            );
        }
    }
}

fn reverse_boundary(curve: &[Point3; 10]) -> [Point3; 10] {
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for i in 0..10 {
        out[i] = curve[9 - i];
    }
    out
}

/// Phase 6c-D Stage 1c: pick the row0 / col9 boundary CPs for
/// sub-patch `k` — either from `input.hex_boundaries[k]` (when
/// the caller has injected substituted analytic boundary CPs) or
/// synthesised from the guide's corner Hermite data via
/// [`k_p_boundary`].
///
/// When `input.hex_boundaries` is `Some` and indexable at `k`, the
/// injected `row0` / `col9` arrays are returned verbatim — exact
/// `sphere ∩ cylinder` great-circle CPs (deg-9 from
/// `BlendBoundary.hex_boundary_cps_deg9`). When `None`, the
/// generic `k_p_boundary` fallback synthesises from corner
/// position + tangent + curvature, the original Phase H behaviour.
///
/// See `docs/decisions/2026-04-29-phase-6c-d-stage-1-boundary-injection.md`.
fn pick_hex_boundary<G: KpGuide>(
    input: &DirectBoundaryInput,
    guide: &G,
    k: usize,
) -> ([Point3; 10], [Point3; 10]) {
    if let Some(ref hex) = input.hex_boundaries
        && let Some(entry) = hex.get(k)
    {
        return (entry.row0, entry.col9);
    }
    (
        k_p_boundary(guide, true, false),
        k_p_boundary(guide, false, true),
    )
}

fn k_p_boundary<G: KpGuide>(guide: &G, is_u: bool, at_end: bool) -> [Point3; 10] {
    let (ua, va, ub, vb) = match (is_u, at_end) {
        (true, false) => (0.0, 0.0, 0.0, 1.0),
        (true, true) => (1.0, 0.0, 1.0, 1.0),
        (false, false) => (0.0, 0.0, 1.0, 0.0),
        (false, true) => (0.0, 1.0, 1.0, 1.0),
    };
    let p_start = guide.point_at(ua, va);
    let p_end = guide.point_at(ub, vb);
    let (t_start, t_end, tt_start, tt_end) = if is_u {
        (
            guide.partial_v(ua, va),
            guide.partial_v(ub, vb),
            guide.second_partial_vv(ua, va),
            guide.second_partial_vv(ub, vb),
        )
    } else {
        (
            guide.partial_u(ua, va),
            guide.partial_u(ub, vb),
            guide.second_partial_uu(ua, va),
            guide.second_partial_uu(ub, vb),
        )
    };
    let mut cps = [Point3::new(0.0, 0.0, 0.0); 10];
    cps[0] = p_start;
    cps[9] = p_end;
    cps[1] = offset(p_start, t_start, 1.0 / 9.0);
    cps[8] = offset(p_end, t_end, -1.0 / 9.0);
    cps[2] = offset2(p_start, t_start, 2.0 / 9.0, tt_start, 1.0 / 72.0);
    cps[7] = offset2(p_end, t_end, -2.0 / 9.0, tt_end, 1.0 / 72.0);
    for i in 3..=6 {
        let t = ((i - 2) as f64) / 5.0;
        cps[i] = Point3::new(
            cps[2].x + t * (cps[7].x - cps[2].x),
            cps[2].y + t * (cps[7].y - cps[2].y),
            cps[2].z + t * (cps[7].z - cps[2].z),
        );
    }
    cps
}

#[inline]
fn offset(p: Point3, v: Vec3, s: f64) -> Point3 {
    Point3::new(p.x + v.x * s, p.y + v.y * s, p.z + v.z * s)
}

#[inline]
fn offset2(p: Point3, v1: Vec3, s1: f64, v2: Vec3, s2: f64) -> Point3 {
    Point3::new(
        p.x + v1.x * s1 + v2.x * s2,
        p.y + v1.y * s1 + v2.y * s2,
        p.z + v1.z * s1 + v2.z * s2,
    )
}

fn assemble_patch(
    row0: &[Point3; 10],
    row9: &[Point3; 10],
    col0: &[Point3; 10],
    col9: &[Point3; 10],
    r1: &Row1Frame,
    r2: &Row2Frame,
) -> RectangularGregoryBiNonic {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    grid[0].copy_from_slice(row0);
    grid[9].copy_from_slice(row9);
    grid[1].copy_from_slice(&r1.row1);
    grid[8].copy_from_slice(&r1.row8);
    grid[2].copy_from_slice(&r2.row2);
    grid[7].copy_from_slice(&r2.row7);
    for i in 0..10 {
        grid[i][0] = col0[i];
        grid[i][9] = col9[i];
        grid[i][1] = r1.col1[i];
        grid[i][8] = r1.col8[i];
        grid[i][2] = r2.col2[i];
        grid[i][7] = r2.col7[i];
    }
    RectangularGregoryBiNonic { grid }
}

#[cfg(test)]
#[path = "kp_ring_build_tests.rs"]
mod tests;
