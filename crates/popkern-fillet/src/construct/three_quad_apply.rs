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

use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::poke_hex::poke_n_gon_face_into_n_quads;
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
) -> Result<Vec<CornerPatchInfo>, ThreeQuadFailure> {
    let decomp = build_n_quad_decomposition_with_options(
        boundary,
        fillet_surfaces,
        use_deg5,
        classification,
    )
    .ok_or_else(|| ThreeQuadFailure::Other("boundary not decomposable into N quads".to_string()))?;

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

#[cfg(test)]
#[path = "three_quad_apply_tests.rs"]
mod tests;
