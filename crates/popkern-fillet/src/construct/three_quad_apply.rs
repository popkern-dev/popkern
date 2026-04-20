//! §5.3 Round 8 C5 — apply a three-quad decomposition to a BRep.
//!
//! Given a BRep with a hexagonal blend face (from the existing
//! `insert_vertex_blend` pipeline), apply the commercial-style
//! three-quad decomposition:
//!   1. Build the [`ThreeQuadDecomposition`] from the hex boundary.
//!   2. Poke the hex face into three quadrilateral sub-faces
//!      via [`poke_hex_face_into_three_quads`].
//!   3. Assign each sub-patch's `RectangularGregory` surface to
//!      the corresponding sub-face.
//!   4. Return three [`CornerPatchInfo`] entries — one per
//!      sub-face — replacing the single entry the caller would
//!      otherwise have emitted.
//!
//! The helper is a **pure function on the BRep + boundary**; it
//! doesn't know about `FilletOptions` flags or dispatch — the
//! caller decides when to invoke it. R8 C6 wires the dispatch.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    DirectBoundaryInput, DirectConstructionOptions, GregoryPatch, OuterCornerData,
    RectangularGregoryTwoStripe, Side, Spoke, SpokeBoundaryData, SpokeEnd, TangentPlane,
    build_binonic_ring,
};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::poke_hex::poke_n_gon_face_into_n_quads;
use super::three_quad_promotions::{
    promote_decomp_to_biquintic_full_with_spoke_g1, promote_decomp_to_biseptic_with_spoke_g1,
    promote_decomp_to_two_stripe_with_spoke_g1,
};
use crate::blend::boundary::BlendBoundary;
use crate::blend::gregory_three_quad::build_n_quad_decomposition_with_options;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod, CornerPatchInfo};
use crate::vertex_classifier::VertexBlendClassification;

/// §5.3 item #10 — why the three-quad dispatch failed. Callers
/// use this to distinguish an honest topology / geometry failure
/// (caller should record the reason + fall through to polar) from
/// a self-intersection remediation (caller must ALSO emit the
/// `CornerRemediationOverlap` warning).
#[derive(Debug)]
pub(crate) enum ThreeQuadFailure {
    /// The Gregory sub-patches self-intersect; remediate to polar.
    Overlap,
    /// Anything else — boundary shape, decomposition rejection,
    /// poke topology failure. Error string included for logging.
    Other(#[allow(dead_code)] String),
}

impl From<String> for ThreeQuadFailure {
    fn from(s: String) -> Self {
        ThreeQuadFailure::Other(s)
    }
}

/// Apply a three-quad decomposition to an existing hexagonal blend
/// face: build the decomposition, poke the face, assign Gregory
/// surfaces, return one [`CornerPatchInfo`] per sub-face.
///
/// # Errors
///
/// Returns an error string if:
/// * The boundary isn't a valid hexagon (6 alternating segments).
/// * The decomposition builder rejects the input.
/// * The poke topology surgery fails.
#[allow(dead_code)] // wired into dispatch in R8 C6
#[allow(clippy::too_many_arguments)]
pub(crate) fn apply_three_quad_decomposition_to_blend(
    brep: &mut BRep,
    blend_face: FaceHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    vertex: VertexHandle,
    n_sides: usize,
    classification: VertexBlendClassification,
    use_deg5: bool,
    use_biquintic_full: bool,
    use_biseptic: bool,
    use_two_stripe: bool,
    use_kp_binonic: bool,
) -> Result<Vec<CornerPatchInfo>, ThreeQuadFailure> {
    let mut decomp = build_n_quad_decomposition_with_options(
        boundary,
        fillet_surfaces,
        use_deg5,
        classification,
    )
    .ok_or_else(|| ThreeQuadFailure::Other("boundary not decomposable into N quads".to_string()))?;

    // §5.3 item #3 extended — promote deg5 sub-patches to the
    // biquintic-full form and run the Karciauskas-Peters spoke-
    // G1 joint solver over the N-quad ring. Biseptic (degree-7)
    // opt-in runs a second pass: degree-elevate each biquintic-
    // full sub-patch to biseptic and re-solve on the wider
    // 4-DOF-per-spoke space, yielding 4–6× further residual
    // reduction (limited by the corner-adjacent `B_1^n + B_{n-1}^n`
    // out-of-span structure — Phase 10 two-stripe / Peters work
    // closes this fully).
    if use_biquintic_full && use_deg5 {
        promote_decomp_to_biquintic_full_with_spoke_g1(&mut decomp.sub_patches);
        // Conservative ladder: only attempt the next rung if the
        // prior one's promotion succeeded AND its solver didn't
        // trigger a fold-revert. This guarantees the scorecard
        // never regresses across rungs.
        if use_biseptic
            && promote_decomp_to_biseptic_with_spoke_g1(&mut decomp.sub_patches)
            && use_two_stripe
        {
            promote_decomp_to_two_stripe_with_spoke_g1(&mut decomp.sub_patches);
            // Phase H — K-P bidegree-9 guided-surface construction.
            // Top rung: extract the DirectBoundaryInput from the
            // two-stripe-promoted decomp, build the bicubic guide
            // ring, derive the row-1 / row-2 frames + interior fit,
            // and replace each sub-patch with a binonic variant.
            // Failure reverts to two-stripe silently (the ladder's
            // gate contract — never regress).
            if use_kp_binonic {
                promote_decomp_to_kp_binonic(&mut decomp.sub_patches);
            }
        }
    }

    // §5.3 item #10 — self-intersection screen before committing
    // the sub-patches to the BRep. If any sub-patch folds or
    // curls, abort so the caller falls through to polar-NURBS
    // (which is single-patch and can't produce the same overlap
    // mode). Signalled as `ThreeQuadFailure::Overlap` so the
    // caller knows to emit the `CornerRemediationOverlap` warning.
    const SELF_INT_TOL: f64 = 1e-6;
    for patch in &decomp.sub_patches {
        if patch.self_intersects_with_tolerance(SELF_INT_TOL) {
            return Err(ThreeQuadFailure::Overlap);
        }
    }

    let poke = poke_n_gon_face_into_n_quads(brep, blend_face, decomp.central_point)
        .map_err(|e| ThreeQuadFailure::Other(format!("poke failed: {e}")))?;
    if poke.sub_faces.len() != decomp.sub_patches.len() {
        return Err(ThreeQuadFailure::Other(format!(
            "poke/decomposition size mismatch: poke produced {} sub-faces, \
             decomposition has {} sub-patches",
            poke.sub_faces.len(),
            decomp.sub_patches.len()
        )));
    }

    let mut patches = Vec::with_capacity(decomp.sub_patches.len());
    for (sub_face, sub_patch) in poke.sub_faces.iter().zip(decomp.sub_patches.iter()) {
        brep.faces[*sub_face].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Gregory(sub_patch.clone()),
            orientation: true,
        });
        patches.push(CornerPatchInfo {
            vertex,
            n_sides,
            method: CornerMethod::Gregory,
            achieved_continuity: ContinuityLevel::G1,
            classification,
            face: *sub_face,
            bridge_fallbacks: boundary.bridge_fallbacks.clone(),
            overlap_remediated: false,
        });
    }
    Ok(patches)
}

/// Phase H — attempt to promote every sub-patch from
/// `RectangularTwoStripe` to `RectangularBiNonic` via the
/// Karciauskas-Peters guided-surface construction. No-op unless
/// every sub-patch is already a two-stripe (the upstream ladder
/// is required; the K-P extraction reads the two-stripe grid).
///
/// Failure modes revert to two-stripe silently:
///   - Non-two-stripe patches present (caller out of order).
///   - `extract_direct_boundary_input` returns `None`
///     (normalise failure on a boundary-tangent — happens on
///     degenerate inputs).
///   - `build_binonic_ring` returns `Err` (e.g. guide-ring
///     build rejected the input; upstream ladder remains the
///     best available result).
///   - Any resulting binonic patch self-intersects.
fn promote_decomp_to_kp_binonic(sub_patches: &mut [GregoryPatch]) {
    let n = sub_patches.len();
    if n == 0 {
        return;
    }
    // Every sub-patch must already be two-stripe — we read the
    // grid as the source of the K-P DirectBoundaryInput.
    let mut two_stripe: Vec<RectangularGregoryTwoStripe> = Vec::with_capacity(n);
    for patch in sub_patches.iter() {
        let GregoryPatch::RectangularTwoStripe(ts) = patch else {
            return;
        };
        two_stripe.push((**ts).clone());
    }

    // Ring-topology spokes: consistent with the other promotion
    // rungs — spoke `k` shares sub-patch `k`'s UOne with sub-patch
    // `(k+1) % n`'s VZero (reversed).
    let spokes: Vec<Spoke> = (0..n)
        .map(|k| Spoke {
            a: SpokeEnd {
                patch_idx: k,
                side: Side::UOne,
                reversed: false,
            },
            b: SpokeEnd {
                patch_idx: (k + 1) % n,
                side: Side::VZero,
                reversed: true,
            },
        })
        .collect();

    let Some(input) = extract_direct_boundary_input(&two_stripe, &spokes) else {
        return;
    };
    let char_len = estimate_characteristic_length(&two_stripe);
    let options = DirectConstructionOptions::with_characteristic_length(char_len);
    let Ok(binonic) = build_binonic_ring(&input, &options) else {
        return;
    };

    // Self-intersection gate — never commit a fold.
    for patch in &binonic {
        let gp = GregoryPatch::RectangularBiNonic(Box::new(patch.clone()));
        if gp.self_intersects_with_tolerance(1e-6) {
            return;
        }
    }

    for (idx, patch) in binonic.into_iter().enumerate() {
        sub_patches[idx] = GregoryPatch::RectangularBiNonic(Box::new(patch));
    }
}

/// Build a [`DirectBoundaryInput`] from the biseptic-promoted
/// two-stripe seed patches + ring topology. Extracts the per-
/// spoke boundary CPs, corner xtan values, and input-face
/// tangent planes (approximated from the patch normals at C
/// and V — upstream biseptic promotion already wrote these
/// consistent with the true input face planes, so using the
/// patch's computed normal is equivalent).
fn extract_direct_boundary_input(
    patches: &[RectangularGregoryTwoStripe],
    spokes: &[Spoke],
) -> Option<DirectBoundaryInput> {
    let n = patches.len();
    if n == 0 || spokes.len() != n {
        return None;
    }
    // Central vertex C: patch 0's grid[7][0] (u=1, v=0) per the
    // three-quad convention.
    let central_vertex = patches[0].grid[7][0];

    let mut spoke_data = Vec::with_capacity(n);
    let mut outer_corners = Vec::with_capacity(n);
    let mut central_vertex_tps = Vec::with_capacity(n);

    for (k, spoke) in spokes.iter().enumerate() {
        let patch_a = &patches[spoke.a.patch_idx];
        let patch_b = &patches[spoke.b.patch_idx];

        // Spoke boundary CPs: for the standard three-quad (A UOne
        // not reversed), these are patch_a.grid[7][0..=7].
        let mut boundary_cps = [Point3::new(0.0, 0.0, 0.0); 8];
        for (j, cp) in boundary_cps.iter_mut().enumerate() {
            *cp = match spoke.a.side {
                Side::UOne => patch_a.grid[7][j],
                Side::UZero => patch_a.grid[0][j],
                Side::VOne => patch_a.grid[j][7],
                Side::VZero => patch_a.grid[j][0],
            };
        }

        // xtan at corner C (spoke-t = 0) and V (spoke-t = 1) on
        // each side. Side A: partial in the cross-boundary
        // direction at the spoke's side. Side B: same, but with
        // reversal mapping spoke-t → B's local parameter via
        // t_local = 1 - t when b.reversed is true.
        let (ua_c, va_c) = uv_at_spoke_t(spoke.a.side, 0.0, spoke.a.reversed);
        let (ua_v, va_v) = uv_at_spoke_t(spoke.a.side, 1.0, spoke.a.reversed);
        let xtan_a_at_c = cross_tangent_on_side(patch_a, spoke.a.side, ua_c, va_c);
        let xtan_a_at_v = cross_tangent_on_side(patch_a, spoke.a.side, ua_v, va_v);

        let (ub_c, vb_c) = uv_at_spoke_t(spoke.b.side, 0.0, spoke.b.reversed);
        let (ub_v, vb_v) = uv_at_spoke_t(spoke.b.side, 1.0, spoke.b.reversed);
        let xtan_b_at_c = cross_tangent_on_side(patch_b, spoke.b.side, ub_c, vb_c);
        let xtan_b_at_v = cross_tangent_on_side(patch_b, spoke.b.side, ub_v, vb_v);

        spoke_data.push(SpokeBoundaryData {
            boundary_cps,
            xtan_a_at_c,
            xtan_b_at_c,
            xtan_a_at_v,
            xtan_b_at_v,
            // Radii not used by the current blend-coeffs formula
            // (constant-α form). Populate with a plausible placeholder
            // so the struct validates; the K-K-C 1991 variable-
            // radius extension will consume these later.
            radius_at_c: 1.0,
            radius_at_v: 1.0,
            edge_id: format!("spoke_{k}"),
        });

        // Outer corner position: the shared (u=1, v=1) of patch A.
        let v_pos = patch_a.grid[7][7];
        // Tangent plane at V from patch A's normal there.
        let du_v = patch_a.partial_u(1.0, 1.0);
        let dv_v = patch_a.partial_v(1.0, 1.0);
        let n_v = du_v.cross(&dv_v).normalized()?;
        let tps_v = TangentPlane::from_normal(n_v)?;
        outer_corners.push(OuterCornerData {
            position: v_pos,
            tangent_plane: tps_v,
        });

        // Central-vertex tangent plane from patch A's normal at C.
        let du_c = patch_a.partial_u(1.0, 0.0);
        let dv_c = patch_a.partial_v(1.0, 0.0);
        let n_c = du_c.cross(&dv_c).normalized()?;
        let tps_c = TangentPlane::from_normal(n_c)?;
        central_vertex_tps.push(tps_c);
    }

    Some(DirectBoundaryInput {
        central_vertex,
        spokes: spoke_data,
        outer_corners,
        central_vertex_tps,
        seed_patches: Some(patches.to_vec()),
    })
}

/// Map a spoke parameter `t ∈ [0, 1]` to the `(u, v)` on a
/// side-indexed patch boundary, honouring the `reversed` flag
/// that three-quad uses for its VZero B-ends (spoke-t = 0 ↔
/// patch-u = 1 when reversed).
fn uv_at_spoke_t(side: Side, t: f64, reversed: bool) -> (f64, f64) {
    let t_local = if reversed { 1.0 - t } else { t };
    match side {
        Side::UZero => (0.0, t_local),
        Side::UOne => (1.0, t_local),
        Side::VZero => (t_local, 0.0),
        Side::VOne => (t_local, 1.0),
    }
}

/// Partial derivative in the cross-boundary direction at the
/// given side. For U-sides that's `∂S/∂u`; for V-sides `∂S/∂v`.
fn cross_tangent_on_side(
    patch: &RectangularGregoryTwoStripe,
    side: Side,
    u: f64,
    v: f64,
) -> popkern_core::geometry::Vec3 {
    match side {
        Side::UZero | Side::UOne => patch.partial_u(u, v),
        Side::VZero | Side::VOne => patch.partial_v(u, v),
    }
}

/// Estimate a characteristic length for tolerance scaling — the
/// max bounding-box extent across the seed patches. Populates
/// `DirectConstructionOptions::characteristic_length` so all
/// ratio-based tolerances scale with the geometry's natural
/// size.
fn estimate_characteristic_length(patches: &[RectangularGregoryTwoStripe]) -> f64 {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for patch in patches {
        for i in 0..8 {
            for j in 0..8 {
                let p = patch.grid[i][j];
                for (k, v) in [p.x, p.y, p.z].iter().enumerate() {
                    if *v < min[k] {
                        min[k] = *v;
                    }
                    if *v > max[k] {
                        max[k] = *v;
                    }
                }
            }
        }
    }
    let dx = max[0] - min[0];
    let dy = max[1] - min[1];
    let dz = max[2] - min[2];
    (dx * dx + dy * dy + dz * dz).sqrt().max(1e-9)
}

#[cfg(test)]
#[path = "three_quad_apply_tests.rs"]
mod tests;
