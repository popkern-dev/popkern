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
use super::kp_guide_patch_biquintic::KpGuidePatchBiQuintic;
use super::kp_guide_ring::build_guide_ring_with_twist;
use super::kp_guide_trait::KpGuide;
use super::kp_interior_fit::{fit_interior_cps_to_biquintic_guide, fit_interior_cps_to_guide};
use super::kp_row1::{Row1Frame, derive_row1_frame};
use super::kp_row2::{Row2Frame, derive_row2_frame};
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
        let row0 = k_p_boundary(&guides[k], true, false);
        let row9 = shared_spokes[k];
        let col0 = reverse_boundary(&shared_spokes[k_prev]);
        let col9 = k_p_boundary(&guides[k], false, true);
        boundaries.push([row0, row9, col0, col9]);
    }

    // H.4.d: row-1 frames.
    let mut r1_frames: Vec<Row1Frame> = guides
        .iter()
        .zip(boundaries.iter())
        .map(|(g, b)| derive_row1_frame(g, &b[0], &b[1], &b[2], &b[3]))
        .collect();

    // H.4.e: row-2 frames.
    let mut r2_frames: Vec<Row2Frame> = guides
        .iter()
        .zip(boundaries.iter())
        .zip(r1_frames.iter())
        .map(|((g, b), r1)| derive_row2_frame(g, &b[0], &b[1], &b[2], &b[3], r1))
        .collect();

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
        fit_interior_cps_to_guide(&mut patch, &guides[k]);
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
        let row0 = k_p_boundary(&guides[k], true, false);
        let row9 = shared_spokes[k];
        let col0 = reverse_boundary(&shared_spokes[k_prev]);
        let col9 = k_p_boundary(&guides[k], false, true);
        boundaries.push([row0, row9, col0, col9]);
    }
    let mut r1_frames: Vec<Row1Frame> = guides
        .iter()
        .zip(boundaries.iter())
        .map(|(g, b)| derive_row1_frame(g, &b[0], &b[1], &b[2], &b[3]))
        .collect();
    let mut r2_frames: Vec<Row2Frame> = guides
        .iter()
        .zip(boundaries.iter())
        .zip(r1_frames.iter())
        .map(|((g, b), r1)| derive_row2_frame(g, &b[0], &b[1], &b[2], &b[3], r1))
        .collect();

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
        fit_interior_cps_to_biquintic_guide(&mut patch, &guides[k]);
        patches.push(patch);
    }
    Ok(patches)
}

fn reverse_boundary(curve: &[Point3; 10]) -> [Point3; 10] {
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for i in 0..10 {
        out[i] = curve[9 - i];
    }
    out
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
