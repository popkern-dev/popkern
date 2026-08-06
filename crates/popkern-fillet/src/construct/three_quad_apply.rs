//! §5.3 Round 8 C5 — apply a three-quad decomposition to a BRep.
//!
//! Applies an N-sector Gregory decomposition to an existing blend
//! face: build patches, poke topology, assign surfaces, and return
//! one [`CornerPatchInfo`] per sub-face.

#[path = "three_quad_apply_helpers.rs"]
pub(super) mod three_quad_apply_helpers;
#[path = "three_quad_blend_seam_reanchor.rs"]
mod three_quad_blend_seam_reanchor;
#[path = "three_quad_face_matching.rs"]
mod three_quad_face_matching;
#[path = "three_quad_kp_boundary_relock.rs"]
mod three_quad_kp_boundary_relock;
#[path = "three_quad_kp_promote.rs"]
mod three_quad_kp_promote;
#[path = "three_quad_planned_seam_authoring.rs"]
mod three_quad_planned_seam_authoring;
#[path = "three_quad_reproject.rs"]
mod three_quad_reproject;
#[path = "three_quad_seam_pcurves.rs"]
mod three_quad_seam_pcurves;
#[path = "three_quad_sphere_cap.rs"]
mod three_quad_sphere_cap;
#[path = "three_quad_trace.rs"]
mod three_quad_trace;

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    DirectBoundaryInput, DirectConstructionOptions, GregoryPatch, OuterCornerData,
    RectangularGregoryTwoStripe, Side, Spoke, SpokeBoundaryData, SpokeEnd, TangentPlane,
};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};
pub(crate) use three_quad_apply_helpers::{
    ThreeQuadApplyResult, ThreeQuadFailure, author_cap_face_boundary_pcurves,
    author_model_level_support_face_boundary_seams,
};
use three_quad_apply_helpers::{
    decomp_self_intersects, macro_family_from_promotion, max_outer_boundary_residual,
    record_support_cell_nonselection, sub_patch_provenance,
};
use three_quad_face_matching::boundary_segment_to_sub_face_permutation;
use three_quad_kp_promote::promote_decomp_to_kp_binonic;
pub(crate) use three_quad_kp_promote::{extract_direct_boundary_input, worst_ring_coplanarity};
use three_quad_reproject::reproject_loop_pcurves_onto_gregory;
use three_quad_seam_pcurves::author_boundary_owned_seam_pcurves;
use three_quad_sphere_cap::{
    cap_spheres_for_subfaces, severe_radius_spread,
    sphere_cap_boundary_degeneracy_in_boundary_order, sphere_cap_boundary_seams_are_nondegenerate,
};
use three_quad_trace::trace_cap_subface_neighbors;

use super::poke_hex::poke_n_gon_face_into_expected_n_quads;
use super::promotion_ladder::PromotionResult;
use super::sphere_fit_edge_wire_cap_identity::semantic_cap_indices_for_subfaces;
use crate::blend::boundary::{BlendBoundary, BoundarySegment};
use crate::blend::gregory_three_quad::build_n_quad_decomposition_with_options;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod, CornerPatchInfo};
use crate::vertex_classifier::VertexBlendClassification;

/// Apply a three-quad decomposition to an existing hexagonal blend
/// face: build the decomposition, poke the face, assign Gregory
/// surfaces, return one [`CornerPatchInfo`] per sub-face. Errors on
/// a non-2N boundary, decomposition rejection, or poke failure.
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
    use_kp_polar: bool,
) -> Result<ThreeQuadApplyResult, ThreeQuadFailure> {
    let mut decomp = build_n_quad_decomposition_with_options(
        boundary,
        fillet_surfaces,
        use_deg5,
        classification,
    )
    .or_else(|| {
        if !use_deg5 {
            return None;
        }
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace §DirectApex.4/gregory-nquad] degree-5 decomposition rejected boundary; \
                 retrying validated cubic decomposition"
            );
        }
        build_n_quad_decomposition_with_options(boundary, fillet_surfaces, false, classification)
    })
    .ok_or_else(|| ThreeQuadFailure::Other("boundary not decomposable into N quads".to_string()))?;

    // Phase H.8 step-D — relaxed-continuity carry-through. `Some(c)`
    // after K-P promotion means the binonic shipped but worst-case
    // coplanarity `c` exceeded the 5e-2 commercial gate.
    let mut relaxed_coplanarity: Option<f64> = None;

    // Capture or fit hex-side boundary CPs so K-P can measure exact
    // live-boundary ownership instead of guide-only Hermite synthesis.
    let fitted_hex_boundary_cps_deg9 = if boundary.hex_boundary_cps_deg9.is_empty() {
        super::three_quad_boundary_fit::fit_hex_boundary_cps_deg9(boundary).unwrap_or_default()
    } else {
        Vec::new()
    };
    let hex_boundary_data: Option<&[([Point3; 10], [Point3; 10])]> =
        if !boundary.hex_boundary_cps_deg9.is_empty() {
            Some(&boundary.hex_boundary_cps_deg9)
        } else if !fitted_hex_boundary_cps_deg9.is_empty() {
            Some(&fitted_hex_boundary_cps_deg9)
        } else {
            None
        };
    let unpromoted_constrained_decomp = decomp.clone();
    let mut diagnostics = Vec::new();

    let mut promotion_report = super::promotion_ladder::run_promotion_ladder(
        &mut decomp,
        fillet_surfaces,
        Some((boundary, &*brep)),
        &mut relaxed_coplanarity,
        use_deg5,
        use_biquintic_full,
        use_biseptic,
        use_two_stripe,
        use_kp_binonic,
        use_kp_polar,
        |patches, fillets| promote_decomp_to_kp_binonic(patches, fillets, hex_boundary_data),
    );
    // PK-240 fold-triggered balanced-spoke retry (see the module
    // doc); on decline the straight result stands unchanged.
    super::three_quad_balanced_retry::apply_fold_triggered_balanced_retry(
        &mut decomp,
        &mut promotion_report,
        &mut relaxed_coplanarity,
        &mut diagnostics,
        boundary,
        brep,
        fillet_surfaces,
        classification,
        super::three_quad_balanced_retry::BalancedRetryLadderFlags {
            use_deg5,
            use_biquintic_full,
            use_biseptic,
            use_two_stripe,
            use_kp_binonic,
            use_kp_polar,
        },
        |patches, fillets| promote_decomp_to_kp_binonic(patches, fillets, hex_boundary_data),
    );
    diagnostics.extend(std::mem::take(&mut promotion_report.diagnostics));
    let mut macro_family = macro_family_from_promotion(promotion_report.result);

    const PROMOTED_BOUNDARY_RESIDUAL_GATE: f64 = 5.0e-3;
    let promoted_boundary_residual = max_outer_boundary_residual(&decomp, boundary);
    let unpromoted_boundary_residual =
        max_outer_boundary_residual(&unpromoted_constrained_decomp, boundary);
    if promoted_boundary_residual.is_some_and(|residual| residual > PROMOTED_BOUNDARY_RESIDUAL_GATE)
        && unpromoted_boundary_residual
            .is_some_and(|residual| residual <= PROMOTED_BOUNDARY_RESIDUAL_GATE)
    {
        diagnostics.push(
            super::three_quad_kp_diagnostics::outer_boundary_residual_exceeded(
                promoted_boundary_residual.unwrap_or_default(),
                unpromoted_boundary_residual,
                PROMOTED_BOUNDARY_RESIDUAL_GATE,
            ),
        );
        decomp = unpromoted_constrained_decomp.clone();
        record_support_cell_nonselection(
            &mut diagnostics,
            "promoted_boundary_residual_exceeded",
            n_sides,
            decomp.sub_patches.len(),
        );
        promotion_report.result = PromotionResult::Standard;
        macro_family = None;
    }

    let sphere_cap_available = boundary.sphere_fit_substituted
        && !boundary.per_fillet_caps.is_empty()
        && boundary.per_fillet_caps.len() == decomp.sub_patches.len();
    let sphere_cap_unavailable_reason = (!sphere_cap_available)
        .then(|| "per-fillet sphere caps unavailable for this decomposition".to_string());
    let sphere_cap_boundary_rejection = if sphere_cap_available {
        sphere_cap_boundary_degeneracy_in_boundary_order(brep, boundary, fillet_surfaces)
    } else {
        None
    };
    if let Some(rejection) = &sphere_cap_boundary_rejection {
        diagnostics.push(rejection.to_diagnostic());
    }
    let sphere_cap_boundary_preflight =
        sphere_cap_available && sphere_cap_boundary_rejection.is_none();
    let sphere_cap_rejection_reason = || {
        sphere_cap_boundary_rejection
            .as_ref()
            .map(|rejection| rejection.describe())
            .or_else(|| sphere_cap_unavailable_reason.clone())
            .unwrap_or_else(|| "available".to_string())
    };

    // Self-intersection screen. N ≥ 5 bypass: bi-5 parallelogram
    // folds on narrow 72°/60° sectors by design (K-P 2015 motivates
    // the polar construction precisely for this regime); downstream
    // polar / KP13 bi-5 paths use only shared-edge CPs. The same
    // containment applies when the planned sphere-fit boundary carries
    // explicit per-fillet analytic sphere caps: the Gregory surface is
    // only an intermediate used to drive the topological sub-face split
    // and will not be committed.
    const SELF_INT_TOL: f64 = 1e-6;
    let mut intermediate_self_intersects = false;
    if decomp.sub_patches.len() < 5 {
        let promoted_self_intersects = decomp_self_intersects(&decomp, SELF_INT_TOL);
        if promoted_self_intersects {
            intermediate_self_intersects = true;
            if !sphere_cap_boundary_preflight {
                if !decomp_self_intersects(&unpromoted_constrained_decomp, SELF_INT_TOL) {
                    decomp = unpromoted_constrained_decomp.clone();
                    promotion_report.result = PromotionResult::Standard;
                    macro_family = None;
                } else {
                    return Err(ThreeQuadFailure::Other(format!(
                        "cap surface candidate planner rejected all candidates: \
                             sphere caps rejected ({sphere_reason}); promoted constrained patch \
                             rejected (self-intersection); unpromoted constrained patch rejected \
                             (self-intersection)",
                        sphere_reason = sphere_cap_rejection_reason()
                    )));
                }
            }
        }
    }

    if matches!(
        promotion_report.result,
        PromotionResult::ExpandedKp20Support | PromotionResult::Kp20SupportCell
    ) && let Some(result) =
        super::three_quad_expanded_kp20_apply::try_apply_kp20_specialized_patchwork(
            super::three_quad_expanded_kp20_apply::ExpandedKp20ApplyInput {
                brep,
                blend_face,
                central_point: decomp.central_point,
                support_surfaces: &decomp.sub_patches,
                support_cell_floor_surfaces: promotion_report
                    .kp20_support_cell_floor_patches
                    .as_deref(),
                support_cell_provenance: promotion_report.kp20_support_cell_provenance.as_deref(),
                promotion_result: promotion_report.result,
                vertex,
                n_sides,
                classification,
                bridge_fallbacks: &boundary.bridge_fallbacks,
                relaxed_coplanarity,
                diagnostics: std::mem::take(&mut diagnostics),
            },
        )
        .map_err(ThreeQuadFailure::Other)?
    {
        return Ok(result);
    }

    let sub_faces_raw: Vec<FaceHandle> =
        if matches!(promotion_report.result, PromotionResult::Macro4n(_)) {
            super::poke_4n::poke_n_gon_face_into_4n_quads(brep, blend_face, decomp.central_point)
                .map_err(|e| ThreeQuadFailure::Other(format!("4n poke failed: {e}")))?
                .sub_faces
        } else {
            poke_n_gon_face_into_expected_n_quads(brep, blend_face, decomp.central_point, n_sides)
                .map_err(|e| ThreeQuadFailure::Other(format!("poke failed: {e}")))?
                .sub_faces
        };
    if sub_faces_raw.len() != decomp.sub_patches.len() {
        return Err(ThreeQuadFailure::Other(format!(
            "poke/decomposition size mismatch: poke produced {} sub-faces, \
             decomposition has {} sub-patches",
            sub_faces_raw.len(),
            decomp.sub_patches.len()
        )));
    }
    if crate::direct_apex_trace::enabled() {
        trace_cap_subface_neighbors(brep, "raw", &sub_faces_raw);
    }
    // Re-order by cylinder identity so sub-face order matches boundary order.
    let sub_faces = match boundary_segment_to_sub_face_permutation(
        brep,
        &sub_faces_raw,
        boundary,
        fillet_surfaces,
    ) {
        Some(perm) => perm
            .into_iter()
            .map(|i| sub_faces_raw[i])
            .collect::<Vec<_>>(),
        None => sub_faces_raw,
    };
    if crate::direct_apex_trace::enabled() {
        trace_cap_subface_neighbors(brep, "ordered", &sub_faces);
    }
    brep.normalize_edge_representatives();
    let semantic_cap_indices = semantic_cap_indices_for_subfaces(brep, boundary, &sub_faces);

    // Write substituted boundary curves into BRep edges when available.
    if boundary.sphere_fit_substituted {
        let wired = super::sphere_fit_edge_wire::wire_substituted_boundary_to_edges(
            brep, boundary, &sub_faces,
        );
        if std::env::var_os("POPKERN_KP_TRACE").is_some() {
            eprintln!(
                "[wire_substituted_boundary_to_edges] wrote {wired} of {} segments to BRep edges",
                boundary.segments.len()
            );
        }
    } else {
        // Re-anchor the mitred seam ring onto the blend boundary the
        // promoted cap interpolates (no-op when already anchored).
        let reanchor = three_quad_blend_seam_reanchor::reanchor_blend_seams_to_boundary(
            brep, boundary, &sub_faces,
        );
        diagnostics.push(three_quad_blend_seam_reanchor::reanchor_diagnostic(
            &reanchor, n_sides,
        ));
    }

    // Support-frame rewriting is diagnostic enrichment, not accepted
    // construction: on cube_mixed it can satisfy a local cap/support
    // metric while breaking ring G0 or forcing DirectApex fallback. Keep it out
    // of the construction path until support continuation is rebuilt as
    // a postcondition-gated surface primitive.
    let support_frames_enabled = false;

    // Always compute per-fillet sphere data when sphere-fit substitution
    // succeeded, even if the sphere-cap rung is not selected.
    let cap_sphere_per_subface: Vec<Option<(Point3, f64)>> = if boundary.sphere_fit_substituted
        && !boundary.per_fillet_caps.is_empty()
        && boundary.per_fillet_caps.len() == sub_faces.len()
    {
        let mut matched = cap_spheres_for_subfaces(brep, &sub_faces, &boundary.per_fillet_caps);
        for (k, entry) in matched.iter_mut().enumerate() {
            if entry.is_none() {
                *entry = Some(boundary.per_fillet_caps[k]);
            }
        }
        matched
    } else {
        vec![None; sub_faces.len()]
    };

    // Per-fillet sphere-cap rung: only replace Gregory sub-faces when
    // every planned seam is non-degenerate on the sphere and neighbor.
    let sphere_cap_requested = boundary.sphere_fit_substituted
        && !boundary.per_fillet_caps.is_empty()
        && boundary.per_fillet_caps.len() == sub_faces.len()
        && (intermediate_self_intersects
            || severe_radius_spread(&boundary.per_fillet_caps)
            || std::env::var("POPKERN_SPHERE_CAP_RUNG").is_ok_and(|v| v == "1"));
    let sphere_cap_rung = sphere_cap_requested
        && sphere_cap_boundary_preflight
        && sphere_cap_boundary_seams_are_nondegenerate(
            brep,
            boundary,
            fillet_surfaces,
            &cap_sphere_per_subface,
        );

    let mut patches = Vec::with_capacity(decomp.sub_patches.len());
    let brep_before_surface_assignment = brep.clone();
    for (k, sub_face) in sub_faces.iter().enumerate() {
        let semantic_cap_index = semantic_cap_indices
            .get(k)
            .and_then(|cap_index| *cap_index)
            .unwrap_or(k);
        let Some(sub_patch) = decomp.sub_patches.get(semantic_cap_index) else {
            *brep = brep_before_surface_assignment.clone();
            return Err(ThreeQuadFailure::Other(format!(
                "semantic cap index {semantic_cap_index} has no Gregory sub-patch"
            )));
        };
        // Internal cap-to-cap edges may lack a 3-D curve, so keep the old
        // surface as the reprojection fallback before replacing it.
        let old_surface = brep
            .faces
            .get(*sub_face)
            .and_then(|f| f.surface.as_ref())
            .map(|sof| sof.geometry.clone());
        let new_geometry = if sphere_cap_rung {
            let chosen =
                cap_sphere_per_subface[k].unwrap_or(boundary.per_fillet_caps[semantic_cap_index]);
            SurfaceGeometry::Sphere {
                center: chosen.0,
                radius: chosen.1,
            }
        } else {
            SurfaceGeometry::Gregory(sub_patch.clone())
        };
        brep.faces[*sub_face].surface = Some(SurfaceOnFace {
            geometry: new_geometry,
            orientation: true,
        });
        if sphere_cap_rung
            && let Err(err) =
                three_quad_planned_seam_authoring::author_planned_sphere_cap_boundary_seams(
                    brep,
                    *sub_face,
                    boundary,
                    fillet_surfaces,
                    semantic_cap_index,
                )
        {
            *brep = brep_before_surface_assignment.clone();
            return Err(ThreeQuadFailure::Other(err));
        }
        // Re-fit each boundary pcurve into the new surface's UV space.
        reproject_loop_pcurves_onto_gregory(
            brep,
            *sub_face,
            old_surface.as_ref(),
            boundary.sphere_fit_substituted,
        );
        if !sphere_cap_rung {
            diagnostics.extend(author_boundary_owned_seam_pcurves(
                brep,
                *sub_face,
                boundary,
                fillet_surfaces,
                k,
                Some(semantic_cap_index),
            ));
        }
        if support_frames_enabled && !sphere_cap_rung {
            let surface_before_support_frame = brep.faces[*sub_face].surface.clone();
            super::boundary_support_frames::apply_face_owned_support_frames(brep, *sub_face);
            let support_frame_rewrite_survived = matches!(
                brep.faces.get(*sub_face).and_then(|f| f.surface.as_ref()),
                Some(SurfaceOnFace {
                    geometry: SurfaceGeometry::Gregory(patch),
                    ..
                }) if !patch.self_intersects_with_tolerance(SELF_INT_TOL)
            );
            if !support_frame_rewrite_survived {
                brep.faces[*sub_face].surface = surface_before_support_frame;
            }
        }
        three_quad_reproject::author_missing_outer_loop_pcurves_from_vertices(brep, *sub_face);
        super::three_quad_support_pcurve_completion::complete_missing_cap_support_seams(
            brep, *sub_face,
        );
        if let Err(err) = three_quad_reproject::verify_outer_loop_pcurves(brep, *sub_face) {
            *brep = brep_before_surface_assignment.clone();
            return Err(ThreeQuadFailure::Other(format!(
                "cap surface assignment left invalid trim pcurves: {err}"
            )));
        }
        if crate::direct_apex_trace::enabled() {
            trace_cap_subface_neighbors(brep, "post-assign", &[*sub_face]);
        }
        patches.push(CornerPatchInfo {
            vertex,
            n_sides,
            method: CornerMethod::Gregory,
            achieved_continuity: ContinuityLevel::G1,
            classification,
            face: *sub_face,
            bridge_fallbacks: boundary.bridge_fallbacks.clone(),
            overlap_remediated: false,
            relaxed_continuity_coplanarity: relaxed_coplanarity,
            cap_sphere: cap_sphere_per_subface[k],
            sub_patch_provenance: sub_patch_provenance(macro_family, n_sides, semantic_cap_index),
            setback_cap: false,
            exact_boundary_cap: false,
        });
    }
    // PK-240: a balanced fan curves the interior spokes; re-derive
    // each interior seam's edge triple from the committed binonic
    // sides (straight fans keep their placeholder seams untouched).
    if decomp.balanced_spokes.is_some() && !sphere_cap_rung {
        let _ = super::three_quad_spoke_seam_curves::
            author_interior_spoke_seam_curves_from_committed_caps(brep, &sub_faces);
    }
    if crate::direct_apex_trace::enabled() {
        trace_cap_subface_neighbors(brep, "final", &sub_faces);
    }
    Ok(ThreeQuadApplyResult {
        patches,
        diagnostics,
    })
}
#[cfg(test)]
#[path = "three_quad_apply_tests.rs"]
mod tests;
