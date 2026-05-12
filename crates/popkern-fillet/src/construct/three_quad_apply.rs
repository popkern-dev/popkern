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

#[path = "three_quad_face_matching.rs"]
mod three_quad_face_matching;
#[path = "three_quad_kp_promote.rs"]
mod three_quad_kp_promote;
#[path = "three_quad_planned_seam_authoring.rs"]
mod three_quad_planned_seam_authoring;
#[path = "three_quad_reproject.rs"]
mod three_quad_reproject;
#[path = "three_quad_seam_pcurves.rs"]
mod three_quad_seam_pcurves;

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    DirectBoundaryInput, DirectConstructionOptions, GregoryPatch, OuterCornerData,
    RectangularGregoryTwoStripe, Side, Spoke, SpokeBoundaryData, SpokeEnd, TangentPlane,
};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};
use three_quad_face_matching::{
    boundary_segment_to_sub_face_permutation, primary_cylinder_identity_of_face,
};
use three_quad_kp_promote::promote_decomp_to_kp_binonic;
pub(crate) use three_quad_kp_promote::{extract_direct_boundary_input, worst_ring_coplanarity};
use three_quad_reproject::reproject_loop_pcurves_onto_gregory;
use three_quad_seam_pcurves::author_boundary_owned_seam_pcurves;

use super::poke_hex::poke_n_gon_face_into_n_quads;
use crate::blend::boundary::{BlendBoundary, BoundarySegment};
use crate::blend::gregory_three_quad::build_n_quad_decomposition_with_options;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod, CornerPatchInfo, CornerSubPatchProvenance};
use crate::vertex_classifier::VertexBlendClassification;

/// §5.3 item #10 — why the three-quad dispatch failed. Callers
/// use this to distinguish an honest topology / geometry failure
/// (caller should record the reason + fall through to polar) from
/// a self-intersection remediation (caller must ALSO emit the
/// `CornerRemediationOverlap` warning).
#[derive(Debug)]
pub(crate) enum ThreeQuadFailure {
    /// The Gregory sub-patches self-intersect; remediate to polar.
    #[allow(dead_code)]
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
/// Returns an error string if the boundary isn't a valid hexagon
/// (6 alternating segments), the decomposition builder rejects
/// the input, or the poke topology surgery fails.
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
) -> Result<Vec<CornerPatchInfo>, ThreeQuadFailure> {
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
        if crate::a2_trace::enabled() {
            eprintln!(
                "[a2-trace §A2.4/gregory-nquad] degree-5 decomposition rejected boundary; \
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

    // Phase 6c-D Stage 1d: capture the substituted hex-side
    // boundary CPs from `boundary` so the K-P bi-nonic ring
    // builder can use them via `pick_hex_boundary` instead of
    // `k_p_boundary` Hermite synthesis. When the boundary's
    // `hex_boundary_cps_deg9` is empty (sphere-fit substitution
    // didn't run / declined), the closure passes `None` and
    // existing Hermite-blend behaviour is preserved.
    let hex_boundary_data: Option<&[([Point3; 10], [Point3; 10])]> =
        if boundary.hex_boundary_cps_deg9.is_empty() {
            None
        } else {
            Some(&boundary.hex_boundary_cps_deg9)
        };
    let sphere_cap_available_before_promotion = boundary.sphere_fit_substituted
        && !boundary.per_fillet_caps.is_empty()
        && boundary.per_fillet_caps.len() == decomp.sub_patches.len();
    let unpromoted_sphere_cap_decomp =
        sphere_cap_available_before_promotion.then(|| decomp.clone());

    let mut is_4n_path = super::promotion_ladder::run_promotion_ladder(
        &mut decomp,
        fillet_surfaces,
        &mut relaxed_coplanarity,
        use_deg5,
        use_biquintic_full,
        use_biseptic,
        use_two_stripe,
        use_kp_binonic,
        use_kp_polar,
        |patches, fillets| promote_decomp_to_kp_binonic(patches, fillets, hex_boundary_data),
    );

    let sphere_cap_available = boundary.sphere_fit_substituted
        && !boundary.per_fillet_caps.is_empty()
        && boundary.per_fillet_caps.len() == decomp.sub_patches.len();
    let sphere_cap_boundary_rejection = if sphere_cap_available {
        sphere_cap_boundary_degeneracy_in_boundary_order(brep, boundary, fillet_surfaces)
    } else {
        Some("per-fillet sphere caps unavailable for this decomposition".to_string())
    };
    let sphere_cap_boundary_preflight = sphere_cap_boundary_rejection.is_none();

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
                match unpromoted_sphere_cap_decomp.clone() {
                    Some(unpromoted) if !decomp_self_intersects(&unpromoted, SELF_INT_TOL) => {
                        decomp = unpromoted;
                        is_4n_path = false;
                    }
                    Some(_) => {
                        return Err(ThreeQuadFailure::Other(format!(
                            "cap surface candidate planner rejected all candidates: \
                             sphere caps rejected ({sphere_reason}); promoted constrained patch \
                             rejected (self-intersection); unpromoted constrained patch rejected \
                             (self-intersection)",
                            sphere_reason = sphere_cap_boundary_rejection
                                .as_deref()
                                .unwrap_or("available")
                        )));
                    }
                    None => {
                        return Err(ThreeQuadFailure::Other(format!(
                            "cap surface candidate planner rejected all candidates: \
                             sphere caps rejected ({sphere_reason}); promoted constrained patch \
                             rejected (self-intersection); no unpromoted constrained candidate \
                             was available",
                            sphere_reason = sphere_cap_boundary_rejection
                                .as_deref()
                                .unwrap_or("available")
                        )));
                    }
                }
            }
        }
    }

    let sub_faces_raw: Vec<FaceHandle> = if is_4n_path {
        super::poke_4n::poke_n_gon_face_into_4n_quads(brep, blend_face, decomp.central_point)
            .map_err(|e| ThreeQuadFailure::Other(format!("4n poke failed: {e}")))?
            .sub_faces
    } else {
        poke_n_gon_face_into_n_quads(brep, blend_face, decomp.central_point)
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
    if crate::a2_trace::enabled() {
        trace_cap_subface_neighbors(brep, "raw", &sub_faces_raw);
    }
    // Phase 6c-D Stage 4.7e: re-order `sub_faces` so each
    // `sub_faces[k]` is paired with the `sub_patches[k]` whose
    // row 0 (= FilletEnd[k] in the boundary's segment order)
    // matches that sub_face's actual geometric cyl-twin.  The
    // poke produces sub_faces in BRep half-edge traversal
    // order, which generally does NOT match the boundary's
    // segments[2k] order — without reordering, sub_face[k] is
    // adjacent to fillet[(k+rot) % N]'s cylinder for some
    // unpredictable rotation, and downstream cap surfaces
    // describe the wrong fillet's arc on each edge.
    //
    // We match by cylinder identity (origin on axis, axis dir,
    // radius) — robust to BRep face surface normalisations.
    // Fixtures where the poke happens to align with the
    // boundary order (e.g. cube_uniform with isotropic radii)
    // produce the identity permutation and behave unchanged.
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
    if crate::a2_trace::enabled() {
        trace_cap_subface_neighbors(brep, "ordered", &sub_faces);
    }
    brep.normalize_edge_representatives();

    // Phase 6c-D: write Phase-6c-A-substituted boundary curves
    // into BRep edges so the pcurve refit
    // (`reproject_loop_pcurves_onto_gregory`) propagates them to
    // the cap-vs-fillet seam pcurves. Gated on the boundary's
    // success flag set in `compute_blend_boundary` — when
    // substitution declined, edges are untouched and the legacy
    // hermite-cubic synthesis path runs unchanged.
    //
    // **Known limitation (Phase 6c-D blocker):** the rational-
    // Gregory dispatch this branch sits in fires from the *second*
    // `blend_insert` invocation per vertex, where the substitution
    // typically declines (the BRep face state has been mutated
    // between calls). The substituted boundary from the *first*
    // invocation never reaches us because the first invocation's
    // `face_has_2n_outer_half_edges` check returns `false`. Net:
    // `boundary.sphere_fit_substituted` is `false` here today even
    // for fixtures where Phase 6c-A's math fully applies (e.g.,
    // `cube_uniform`). Remediation = Phase 6c-C (apex-region
    // topology rebuild) so face handles in the boundary stay valid
    // through the surgery.
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
    }

    let support_frames_enabled = hex_boundary_data.is_some();

    // Always compute the per-fillet sphere `(centre, radius)` for each
    // sub-face when sphere-fit substitution succeeded — independent of
    // whether the sphere-cap rung will be applied.  Stashed on
    // `CornerPatchInfo.cap_sphere` so downstream Phase 6c-A reroute
    // and Phase 6c-D Gregory cap rung can identify cap faces and
    // their nominal sphere parameters even when the surface is a
    // Gregory patch.
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

    // Phase 6c-A continuation §4.2 option (c) — per-fillet
    // sphere-cap rung. When the boundary was sphere-fit-substituted,
    // *replace* each sub-face's Gregory surface with the per-fillet
    // `Sphere(cap_centre_k, r_k)` only if every required planned seam
    // can exist as a non-degenerate curve on that sphere and its
    // neighbour surface. The validated seam authoring primitive below
    // remains the commit gate; this preflight prevents a known invalid
    // topology/surface pairing from reaching mutation at all.
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
    for (k, (sub_face, sub_patch)) in sub_faces.iter().zip(decomp.sub_patches.iter()).enumerate() {
        // CAPTURE THE OLD SURFACE BEFORE REPLACING. Internal cap-to-
        // cap edges (created by `mev`/`mef` during the poke) have
        // `edge.curve = None`; we can't re-project from a 3-D curve
        // that doesn't exist. Instead, we sample
        // `old_surface(old_pcurve(t))` to reconstruct the 3-D path
        // of the edge, then project those points onto the new Gregory.
        let old_surface = brep
            .faces
            .get(*sub_face)
            .and_then(|f| f.surface.as_ref())
            .map(|sof| sof.geometry.clone());
        let new_geometry = if sphere_cap_rung {
            let chosen = cap_sphere_per_subface[k].unwrap_or(boundary.per_fillet_caps[k]);
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
                    k,
                )
        {
            *brep = brep_before_surface_assignment.clone();
            return Err(ThreeQuadFailure::Other(err));
        }
        // Re-fit each boundary pcurve into the new surface's UV space.
        // Uses the edge's 3-D curve when present; falls back to
        // `old_surface(old_pcurve)` when not (internal cap-to-cap
        // edges created by `mev`/`mef`).
        reproject_loop_pcurves_onto_gregory(
            brep,
            *sub_face,
            old_surface.as_ref(),
            boundary.sphere_fit_substituted,
        );
        if support_frames_enabled && !sphere_cap_rung {
            author_boundary_owned_seam_pcurves(brep, *sub_face, boundary, fillet_surfaces, k);
        }
        if support_frames_enabled && !sphere_cap_rung {
            super::boundary_support_frames::apply_face_owned_support_frames(brep, *sub_face);
            let Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Gregory(patch),
                ..
            }) = brep.faces.get(*sub_face).and_then(|f| f.surface.as_ref())
            else {
                *brep = brep_before_surface_assignment.clone();
                return Err(ThreeQuadFailure::Other(
                    "support-frame rewrite lost Gregory cap surface".to_string(),
                ));
            };
            if patch.self_intersects_with_tolerance(SELF_INT_TOL) {
                *brep = brep_before_surface_assignment.clone();
                return Err(ThreeQuadFailure::Other(format!(
                    "cap surface candidate planner rejected constrained patch after \
                     support-frame rewrite (self-intersection); sphere caps rejected ({})",
                    sphere_cap_boundary_rejection
                        .as_deref()
                        .unwrap_or("not selected for this candidate")
                )));
            }
        }
        if let Err(err) = three_quad_reproject::verify_outer_loop_pcurves(brep, *sub_face) {
            *brep = brep_before_surface_assignment.clone();
            return Err(ThreeQuadFailure::Other(format!(
                "cap surface assignment left invalid trim pcurves: {err}"
            )));
        }
        if crate::a2_trace::enabled() {
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
            sub_patch_provenance: sub_patch_provenance(is_4n_path, n_sides, k),
        });
    }
    if crate::a2_trace::enabled() {
        trace_cap_subface_neighbors(brep, "final", &sub_faces);
    }
    Ok(patches)
}

fn sub_patch_provenance(
    is_4n_path: bool,
    n_sides: usize,
    sub_patch_index: usize,
) -> Option<CornerSubPatchProvenance> {
    is_4n_path
        .then(|| CornerSubPatchProvenance::kp20(sub_patch_index / 4, n_sides, sub_patch_index % 4))
}

fn severe_radius_spread(per_fillet_caps: &[(Point3, f64)]) -> bool {
    let mut min_radius = f64::INFINITY;
    let mut max_radius = 0.0_f64;
    for &(_, radius) in per_fillet_caps {
        if radius <= 0.0 {
            return false;
        }
        min_radius = min_radius.min(radius);
        max_radius = max_radius.max(radius);
    }
    min_radius.is_finite() && max_radius / min_radius >= 4.0
}

fn decomp_self_intersects(
    decomp: &crate::blend::gregory_three_quad::NQuadDecomposition,
    tolerance: f64,
) -> bool {
    decomp
        .sub_patches
        .iter()
        .any(|patch| patch.self_intersects_with_tolerance(tolerance))
}

fn cap_spheres_for_subfaces(
    brep: &BRep,
    sub_faces: &[FaceHandle],
    per_fillet_caps: &[(Point3, f64)],
) -> Vec<Option<(Point3, f64)>> {
    let mut used = vec![false; per_fillet_caps.len()];
    let mut assigned = vec![None; sub_faces.len()];

    for (face_index, &face) in sub_faces.iter().enumerate() {
        let Some((cyl_origin, cyl_axis, cyl_radius)) =
            primary_cylinder_identity_of_face(brep, face)
        else {
            continue;
        };
        let Some(axis_unit) = cyl_axis.normalized() else {
            continue;
        };
        let mut best: Option<(usize, f64)> = None;
        for (cap_index, &(centre, radius)) in per_fillet_caps.iter().enumerate() {
            if used[cap_index] {
                continue;
            }
            let radius_error = (radius - cyl_radius).abs();
            if radius_error > cyl_radius.max(radius) * 0.05 + 1.0e-9 {
                continue;
            }
            let v = centre - cyl_origin;
            let along = v.dot(&axis_unit);
            let projected = cyl_origin + axis_unit * along;
            let axis_error = centre.distance_to(&projected);
            let score = axis_error + radius_error;
            if best
                .map(|(_, best_score)| score < best_score)
                .unwrap_or(true)
            {
                best = Some((cap_index, score));
            }
        }
        if let Some((cap_index, _)) = best {
            used[cap_index] = true;
            assigned[face_index] = Some(per_fillet_caps[cap_index]);
            if crate::a2_trace::enabled() {
                eprintln!(
                    "[a2-trace] sphere-cap assignment face={face:?} -> cap[{cap_index}] r={:.3}",
                    per_fillet_caps[cap_index].1
                );
            }
        } else if crate::a2_trace::enabled() {
            eprintln!(
                "[a2-trace] sphere-cap assignment face={face:?} unmatched cyl_radius={cyl_radius:.3}"
            );
        }
    }

    assigned
}

fn sphere_cap_boundary_degeneracy_in_boundary_order(
    brep: &BRep,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<String> {
    if boundary.segments.len() != 2 * boundary.per_fillet_caps.len() {
        return Some(format!(
            "boundary segment count {} does not match per-fillet cap count {}",
            boundary.segments.len(),
            boundary.per_fillet_caps.len()
        ));
    }
    for (k, &(center, radius)) in boundary.per_fillet_caps.iter().enumerate() {
        for segment in [&boundary.segments[2 * k], &boundary.segments[2 * k + 1]] {
            if let Some(reason) =
                sphere_cap_segment_degeneracy(brep, fillet_surfaces, center, radius, segment)
            {
                if crate::a2_trace::enabled() {
                    eprintln!(
                        "[a2-trace] sphere-cap preflight declined boundary[{k}] segment={}: {reason}",
                        boundary_segment_kind(segment)
                    );
                }
                return Some(format!(
                    "boundary[{k}] {}: {reason}",
                    boundary_segment_kind(segment)
                ));
            }
        }
    }
    None
}

fn sphere_cap_boundary_seams_are_nondegenerate(
    brep: &BRep,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    cap_sphere_per_subface: &[Option<(Point3, f64)>],
) -> bool {
    if boundary.segments.len() != 2 * cap_sphere_per_subface.len() {
        return false;
    }
    for (k, cap) in cap_sphere_per_subface.iter().enumerate() {
        let Some((center, radius)) = cap else {
            return false;
        };
        for segment in [&boundary.segments[2 * k], &boundary.segments[2 * k + 1]] {
            if let Some(reason) =
                sphere_cap_segment_degeneracy(brep, fillet_surfaces, *center, *radius, segment)
            {
                if crate::a2_trace::enabled() {
                    eprintln!(
                        "[a2-trace] sphere-cap rung declined subface[{k}] segment={}: {reason}",
                        boundary_segment_kind(segment)
                    );
                }
                return false;
            }
        }
    }
    true
}

fn sphere_cap_segment_degeneracy(
    brep: &BRep,
    fillet_surfaces: &[FilletSurfaceResult],
    sphere_center: Point3,
    sphere_radius: f64,
    segment: &BoundarySegment,
) -> Option<String> {
    match segment {
        BoundarySegment::FaceEdge { face, .. } => {
            let surface = brep.faces.get(*face)?.surface.as_ref()?;
            if let SurfaceGeometry::Plane { origin, normal, .. } = &surface.geometry {
                let normal = normal.normalized()?;
                let plane_offset = (sphere_center - *origin).dot(&normal).abs();
                let circle_radius_sq = sphere_radius * sphere_radius - plane_offset * plane_offset;
                let tol = sphere_radius.max(1.0) * 1.0e-8;
                if circle_radius_sq <= tol * tol {
                    return Some(format!(
                        "sphere-plane boundary is tangent or empty \
                         (sphere radius {sphere_radius:.6e}, plane offset {plane_offset:.6e})"
                    ));
                }
            }
        }
        BoundarySegment::FilletEnd { fillet_index, .. } => {
            let Some(fillet) = fillet_surfaces.get(*fillet_index) else {
                return Some(format!("fillet index {fillet_index} missing"));
            };
            if let SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } = &fillet.surface
            {
                let Some(axis_unit) = axis.normalized() else {
                    return Some("fillet cylinder axis is degenerate".to_string());
                };
                let to_center = sphere_center - *origin;
                let along = to_center.dot(&axis_unit);
                let projected = *origin + axis_unit * along;
                let axis_error = sphere_center.distance_to(&projected);
                let radius_error = (sphere_radius - *radius).abs();
                let tol = sphere_radius.max(*radius).max(1.0) * 1.0e-5;
                if axis_error > tol || radius_error > tol {
                    return Some(format!(
                        "sphere-cylinder tangent circle is not well-defined \
                         (axis error {axis_error:.3e}, radius error {radius_error:.3e})"
                    ));
                }
            }
        }
    }
    None
}

fn boundary_segment_kind(segment: &BoundarySegment) -> &'static str {
    match segment {
        BoundarySegment::FilletEnd { .. } => "FilletEnd",
        BoundarySegment::FaceEdge { .. } => "FaceEdge",
    }
}

fn trace_cap_subface_neighbors(brep: &BRep, label: &str, faces: &[FaceHandle]) {
    let kind = |face: FaceHandle| -> &'static str {
        use popkern_core::geometry::surface::SurfaceGeometry;
        match brep
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry)
        {
            Some(SurfaceGeometry::Plane { .. }) => "Plane",
            Some(SurfaceGeometry::Cylinder { .. }) => "Cyl",
            Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
            Some(SurfaceGeometry::Cone { .. }) => "Cone",
            Some(SurfaceGeometry::Torus { .. }) => "Torus",
            Some(SurfaceGeometry::Nurbs(_)) => "Nurbs",
            Some(SurfaceGeometry::Gregory(_)) => "Gregory",
            None => "None",
        }
    };
    for (idx, &face) in faces.iter().enumerate() {
        let Some(face_data) = brep.faces.get(face) else {
            continue;
        };
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut entries = Vec::new();
        for _ in 0..(brep.half_edges.len() + 1) {
            let Some(he) = brep.half_edges.get(cursor) else {
                break;
            };
            let twin_face = brep
                .half_edges
                .get(he.twin)
                .and_then(|t| brep.loops.get(t.loop_))
                .map(|l| l.face);
            let twin_desc = twin_face
                .map(|f| format!("{f:?}:{}", kind(f)))
                .unwrap_or_else(|| "None".to_string());
            entries.push(format!("{:?}->{twin_desc}", he.edge));
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
        eprintln!(
            "[a2-trace] K-P poke {label} subface[{idx}] {face:?}: {}",
            entries.join(", ")
        );
    }
}

#[cfg(test)]
#[path = "three_quad_apply_tests.rs"]
mod tests;
