//! The construct-stage driver: apply computed fillet/chamfer
//! geometry to a BRep (anchor splits, face-trim splits, consumed
//! merges, surface assignment, vertex blends).

use super::*;
use crate::blend_warnings::trace_smooth_blend_mitred_fallback as trace_mitred;
use crate::error::{BlendFailure, BlendWarning};
use crate::types::VertexMode;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use route_decision::{ApexRouteDecisionInput, apex_route_decision};

/// Apply computed fillet/chamfer geometry to a BRep. 3-phase: A perp splits, B multi-chord face splits, C kef merge + surface assign. `excluded_direct_apex_vertices` = PK-250 retry exclusions. `cancel` = §17 cooperative cancellation, polled between stages and corners. `corner_reporter` = per-corner progress inside the dominant vertex-blend loop (PK-276).
#[allow(clippy::too_many_arguments)]
pub fn apply_blend_geometry(
    brep: &BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
    vertex_mode: &VertexMode,
    use_rational_gregory: bool,
    use_deg5_gregory: bool,
    use_biquintic_full_gregory: bool,
    use_biseptic_gregory: bool,
    use_two_stripe_gregory: bool,
    use_kp_binonic_gregory: bool,
    use_kp_polar_gregory: bool,
    use_direct_apex_construction: bool,
    requested_continuity: crate::types::ContinuityLevel,
    excluded_direct_apex_vertices: &std::collections::BTreeSet<VertexHandle>,
    corner_cap: crate::types::CornerCapPreference,
    cancel: Option<&crate::cancellation::CancellationToken>,
    // `+ '_` keeps the trait object's pointee lifetime independent of
    // the reborrow so callers can reuse their reporter after this call.
    mut reporter: Option<&mut (dyn crate::progress::ProgressReporter + '_)>,
) -> Result<ApplyBlendResult, ApplyBlendError> {
    let mut result = brep.clone();
    post_pipeline::determinism_trace_mod::pre_apply_trace(brep, edge_geometries);
    let mut failures: Vec<(EdgeHandle, BlendFailure)> = Vec::new();
    let mut warnings: Vec<(EdgeHandle, BlendWarning)> = Vec::new();
    let mut corner_patches: Vec<crate::types::CornerPatchInfo> = Vec::new();
    let mut repair_report = crate::types::BlendRepairReport::default();
    let mut unresolved_corners: Vec<corner_resolution::UnresolvedCorner> = Vec::new();
    // Per-vertex (boundary, sub_faces) for the post-sliver_absorb wire pass.
    type BoundaryFaceJobs = Vec<(crate::blend::BlendBoundary, Vec<FaceHandle>)>;
    let mut late_wire_jobs: BoundaryFaceJobs = Vec::new();
    let mut support_seam_jobs: BoundaryFaceJobs = Vec::new();
    // §6c-direct-apex dispatch.d — faces authored by `try_direct_apex_construction`,
    // protected so sliver_absorb/cap_insert/late re-wiring don't re-author the cap.
    let mut direct_apex_built_faces: Vec<popkern_core::topology::FaceHandle> = Vec::new();
    let mut pipeline_timing = pipeline_timing::ApplyBlendPipelineTiming::default();
    // A sub-patch family flag keeps the vertex-blend ladder reachable on
    // the corners that are its reaching fixtures: direct-apex promotion
    // stands down for it, and so does the setback rung's oblique-corner
    // eligibility (PK-597).
    let explicit_subpatch_family = use_biquintic_full_gregory
        || use_biseptic_gregory
        || use_two_stripe_gregory
        || use_kp_binonic_gregory
        || use_kp_polar_gregory;
    // PK-603: a stated setback/fitted preference keeps direct apex off
    // the corners its rung can reach, so the rung runs there first.
    let excluded_direct_apex_vertices: &std::collections::BTreeSet<VertexHandle> =
        &excluded_direct_apex_vertices
            .union(&corner_cap_preference::rung_claimed_corners(
                brep,
                vertex_mode,
                edge_geometries,
                corner_cap,
                tolerance,
                explicit_subpatch_family,
            ))
            .copied()
            .collect();
    // PK-127: promotion decided pre-Phase-A; anchors + dispatch share it.
    let mut promoted_apexes = direct_apex_promotion::promoted_apexes_for_mode(
        brep,
        vertex_mode,
        edge_geometries,
        direct_apex_promotion::PromotionFlags {
            use_rational_gregory,
            use_kp_binonic_gregory,
            explicit_subpatch_family_requested: explicit_subpatch_family,
        },
        tolerance,
    );
    promoted_apexes.retain(|vertex| !excluded_direct_apex_vertices.contains(vertex));
    let (mut split_state, ms) = pipeline_timing::timed(|| {
        perpendicular_splits::run_perpendicular_splits(
            brep,
            &mut result,
            edge_geometries,
            tolerance,
            perpendicular_splits::DirectApexRouting {
                global: use_direct_apex_construction,
                promoted: &promoted_apexes,
                excluded: excluded_direct_apex_vertices,
            },
            &mut failures,
        )
    });
    pipeline_timing.perpendicular_split_ms = ms;
    direct_apex_built_faces.extend(split_state.closed_fillet_protected.iter().copied());
    post_pipeline::determinism_trace("post_perp_splits", &result);
    let (trim_outcome, ms) = pipeline_timing::timed(|| {
        face_trim_pipeline::apply_face_trim_pipeline(
            brep,
            &mut result,
            edge_geometries,
            tolerance,
            &split_state.perp_data,
            &mut split_state.descendants,
            &mut failures,
        )
    });
    repair_report.region_resolved_trims = trim_outcome.region_or_curved_trims;
    pipeline_timing.face_trim_ms = ms;
    post_pipeline::determinism_trace("post_face_trim", &result);

    // Hole-rim terminal discs (PK-185/186) join the protected set.
    direct_apex_built_faces.extend(hole_rim_terminal_disc::author_reflex_end_discs(
        &mut result,
        edge_geometries,
        &split_state.perp_data,
        &trim_outcome.fillet_faces,
        tolerance,
        &mut failures,
    ));
    post_pipeline::determinism_trace("post_reflex_discs", &result);
    // Transactional: ANY apply-phase failure fails the whole
    // operation. The prior all-failed-only gate let partial applies
    // proceed with the failures silently discarded — measured: a
    // 6-spec seam-crossing chain returned Ok with ONE torus band and
    // five blends missing, body validating clean (the same silent
    // wrong-shape class the Phase-2 gate closed one layer up).
    if !failures.is_empty() {
        return Err(ApplyBlendError::Failures(failures));
    }
    apply_cancel_poll(cancel)?; // §17 stage-boundary poll (post-trim).

    // Vertex blend handling. `BTreeMap` keyed on VertexHandle for
    // deterministic insertion order across runs. The closure reports whether
    // a cancel poll fired (early return happens outside the timing wrapper).
    let (blend_cancelled, ms) = pipeline_timing::timed(|| {
        if *vertex_mode == VertexMode::Mitred {
            unresolved_corners.extend(corner_resolution::requested_mitre_corners(
                &result,
                edge_geometries,
            ));
        }
        if *vertex_mode == VertexMode::Smooth && edge_geometries.len() > 1 {
            // Closed edges count ONCE at their closure vertex — a
            // closure vertex is not a corner and must not enter the
            // vertex-blend dispatch as one (PK-296).
            let vertex_edges = crate::analysis::vertex_incidences(
                edge_geometries
                    .iter()
                    .map(|eg| (eg.analysis.vertex_start, eg.analysis.vertex_end)),
            );
            // §5.3 item #5 — classify on the ORIGINAL `brep` (pre-Phase
            // A-through-C). Post-modification, original edge handles are
            // stale and `classify_vertex_blend` silently defaults to
            // `SameSideConvex`, masking `OppositeSide` saddles.
            type VC = crate::vertex_classifier::VertexBlendClassification;
            let precomputed: std::collections::BTreeMap<VertexHandle, VC> = vertex_edges
                .iter()
                .map(|(v, eis)| {
                    let edges: Vec<_> = eis
                        .iter()
                        .filter_map(|&i| edge_geometries.get(i).map(|g| g.analysis.edge))
                        .collect();
                    let c = crate::vertex_classifier::classify_vertex_blend(
                        brep, *v, &edges, tolerance,
                    )
                    .classification;
                    (*v, c)
                })
                .collect();
            // PK-278 parallel solve phase (contract: presolve module docs).
            let cap_presolve = setback_cap_presolve::presolve_setback_caps(
                &result,
                &vertex_edges,
                edge_geometries,
                &precomputed,
                tolerance,
                explicit_subpatch_family,
                apply_cancel_poll(cancel).is_err(),
            );
            // PK-276 per-corner progress: corners are the dominant
            // apply cost, so each completed corner ticks the reporter.
            let total_corners = vertex_edges
                .values()
                .filter(|edge_indices| edge_indices.len() >= 2)
                .count();
            let mut corner_index = 0usize;
            for (vertex, edge_indices) in &vertex_edges {
                // §17 per-corner poll: corner construction dominates apply
                // wall-clock — bounds a cancelled run's stretch to one corner.
                if apply_cancel_poll(cancel).is_err() {
                    return true;
                }
                if edge_indices.len() >= 2 {
                    let corner_started = popkern_core::time::Stopwatch::start();
                    let pre = precomputed.get(vertex).copied().unwrap_or(
                        crate::vertex_classifier::VertexBlendClassification::SameSideConvex,
                    );
                    // PK-192: N=3 corners are claimed by the direct-apex
                    // route (PK-099) and never reach the PnS dispatch
                    // (PK-192 measured AGAINST releasing them:
                    // pns_low_valence_gates.rs).
                    // PK-258: the serving-retry exclusion binds the per-apex
                    // flag directly — pruning promotion alone is a no-op when
                    // the global flag is on (measured: byte-identical rerun).
                    let use_direct_apex_for_this_apex = !excluded_direct_apex_vertices
                        .contains(vertex)
                        && (use_direct_apex_construction || promoted_apexes.contains(vertex));
                    let direct_apex_outcome = if use_direct_apex_for_this_apex {
                        Some(direct_apex_dispatch::try_direct_apex_construction(
                            &mut result,
                            *vertex,
                            edge_indices,
                            edge_geometries,
                            &split_state.descendants,
                            tolerance,
                            !use_direct_apex_construction,
                        ))
                    } else {
                        None
                    };
                    let mut direct_apex_succeeded = false;
                    let mut direct_apex_fell_back_to_vertex_blend = false;
                    if let Some(outcome) = direct_apex_outcome {
                        match outcome {
                            Ok(success) => {
                                repair_report.apex_routes.push(apex_route_decision(
                                    ApexRouteDecisionInput {
                                        vertex: *vertex,
                                        edge_indices,
                                        edge_geometries,
                                        route: crate::types::ApexConstructionRoute::DirectApex,
                                        built_faces: success.built_faces.clone(),
                                        corner_patch_count: success.corner_patches.len(),
                                        exact_boundary_sector_count: success
                                            .exact_boundary_sector_count,
                                        exact_boundary_prepared_curve_count: success
                                            .exact_boundary_prepared_curve_count,
                                        exact_boundary_inner_spoke_count: success
                                            .exact_boundary_inner_spoke_count,
                                        exact_boundary_ribbon_count: success
                                            .exact_boundary_ribbon_count,
                                        exact_boundary_solved_surface_count: success
                                            .exact_boundary_solved_surface_count,
                                        exact_boundary_pcurve_pair_count: success
                                            .exact_boundary_pcurve_pair_count,
                                        exact_boundary_performance: success
                                            .exact_boundary_performance
                                            .clone(),
                                        exact_boundary_report: success
                                            .exact_boundary_report
                                            .clone()
                                            .with_route_label("direct_apex"),
                                        fallback_reason: None,
                                    },
                                ));
                                repair_report.apex_diagnostics.extend(success.diagnostics);
                                direct_apex_built_faces.extend(success.built_faces);
                                corner_patches.extend(success.corner_patches);
                                direct_apex_succeeded = true;
                            }
                            Err(error) => {
                                direct_apex_fell_back_to_vertex_blend = true;
                                let warn_edge = edge_geometries[edge_indices[0]].analysis.edge;
                                let reason = error.reason;
                                repair_report.apex_diagnostics.extend(error.diagnostics);
                                repair_report.apex_routes.push(
                                    route_decision::direct_apex_fallback_route(
                                        *vertex,
                                        edge_indices,
                                        edge_geometries,
                                        reason.clone(),
                                    ),
                                );
                                warnings.push((
                                    warn_edge,
                                    BlendWarning::DirectApexFellBackToVertexBlend {
                                        vertex: *vertex,
                                        reason,
                                    },
                                ));
                                // Fall through to the vertex-blend tail below.
                            }
                        }
                    }
                    if direct_apex_succeeded {
                        apply_types::emit_corner_progress(
                            &mut reporter,
                            *vertex,
                            &mut corner_index,
                            total_corners,
                            corner_started,
                        );
                        continue;
                    }
                    if !use_direct_apex_for_this_apex {
                        repair_report.apex_routes.push(
                            route_decision::vertex_blend_dispatch_route(
                                *vertex,
                                edge_indices,
                                edge_geometries,
                            ),
                        );
                    }

                    // PK-516: the star the corner owned before the vertex
                    // blend touched it — the territory the exit gate
                    // measures if no corner blend is built here.
                    let corner_star = corner_resolution::record_unresolved_corner(
                        &result,
                        *vertex,
                        edge_geometries[edge_indices[0]].analysis.edge,
                        corner_resolution::UnresolvedCornerCause::AllOrNothingRollback(Vec::new()),
                    );
                    match blend_insert::insert_vertex_blend(
                        &mut result,
                        *vertex,
                        edge_indices,
                        edge_geometries,
                        tolerance,
                        requested_continuity,
                        use_rational_gregory,
                        use_deg5_gregory,
                        use_biquintic_full_gregory,
                        use_biseptic_gregory,
                        use_two_stripe_gregory,
                        use_kp_binonic_gregory,
                        use_kp_polar_gregory,
                        Some(pre),
                        cap_presolve.get(vertex),
                        corner_cap,
                        explicit_subpatch_family,
                    ) {
                        Ok((infos, boundary_opt, diagnostics)) => {
                            failures.extend(corner_cap_preference::strict_decline_failure(
                                brep,
                                &diagnostics,
                                edge_geometries[edge_indices[0]].analysis.edge,
                            ));
                            if infos.is_empty() {
                                let verdicts =
                                    corner_resolution::route_decline_verdicts(&diagnostics);
                                unresolved_corners.extend(corner_star.map(|mut corner| {
                                    corner.cause = corner_resolution::UnresolvedCornerCause::
                                        AllOrNothingRollback(verdicts);
                                    corner
                                }));
                            }
                            repair_report.apex_diagnostics.extend(diagnostics);
                            route_decision::record_vertex_blend_serve_on_route(
                                &result,
                                &mut repair_report,
                                *vertex,
                                infos.iter().map(|info| info.face).collect(),
                                direct_apex_fell_back_to_vertex_blend,
                            );
                            for emit in [
                                bridge_warnings::emit_bridge_fallback_warnings,
                                bridge_warnings::emit_tangent_corner_warning,
                                bridge_warnings::emit_overlap_remediation_warning,
                                bridge_warnings::emit_corner_continuity_relaxed_warning,
                                bridge_warnings::emit_empty_blend_mitred_warning,
                            ] {
                                emit(
                                    &infos,
                                    *vertex,
                                    edge_indices,
                                    edge_geometries,
                                    &mut warnings,
                                );
                            }
                            if let Some(b) = boundary_opt {
                                let face_handles: Vec<_> = infos.iter().map(|cp| cp.face).collect();
                                if !face_handles.is_empty() {
                                    support_seam_jobs.push((b.clone(), face_handles.clone()));
                                }
                                if b.sphere_fit_substituted {
                                    late_wire_jobs.push((b, face_handles));
                                }
                            }
                            corner_patches.extend(infos);
                        }
                        Err(err) => {
                            trace_mitred(*vertex, &err);
                            unresolved_corners.extend(corner_star.map(|mut corner| {
                                corner.cause =
                                    corner_resolution::UnresolvedCornerCause::SmoothBlendDeclined(
                                        err.clone(),
                                    );
                                corner
                            }));
                            let edge_idx = edge_indices[0];
                            if edge_idx < edge_geometries.len() {
                                warnings.push((
                                    edge_geometries[edge_idx].analysis.edge,
                                    BlendWarning::VertexModeFallback {
                                        requested: VertexMode::Smooth,
                                        actual: VertexMode::Mitred,
                                        vertex: Some(*vertex),
                                    },
                                ));
                            }
                        }
                    }
                    // Corner done (served, fell back, or mitred).
                    apply_types::emit_corner_progress(
                        &mut reporter,
                        *vertex,
                        &mut corner_index,
                        total_corners,
                        corner_started,
                    );
                }
                post_pipeline::determinism_trace(&format!("post_vertex_{vertex:?}"), &result);
            }
        }
        false
    });
    pipeline_timing.vertex_blend_ms = ms;
    if blend_cancelled {
        return Err(ApplyBlendError::Cancelled);
    }
    // A strict corner-cap decline fails the operation (transactional).
    if !failures.is_empty() {
        return Err(ApplyBlendError::Failures(failures));
    }
    // Fillet end runoff: trim strips at planar end faces they
    // overshoot and cut the end faces' sharp corners back to the
    // runoff curve (see end_runoff module docs).
    let setback_blend_present = corner_patches.iter().any(|cp| cp.setback_cap);
    end_runoff_report::apply_and_report(&mut result, setback_blend_present, &mut repair_report);
    post_pipeline::determinism_trace("post_vertex_blend", &result);
    post_pipeline::determinism_face_trace("post_vertex_blend", &result);
    let (_, ms) = pipeline_timing::timed(|| {
        post_pipeline::determinism_trace("pre_post_pipeline", &result);
        let setback_blend_present = corner_patches.iter().any(|cp| cp.setback_cap);
        post_pipeline::run_post_blend_pipeline(post_pipeline::PostBlendPipelineInputs {
            source_brep: brep,
            setback_blend_present,
            result: &mut result,
            edge_geometries,
            tolerance,
            late_wire_jobs: &late_wire_jobs,
            support_seam_jobs: &support_seam_jobs,
            direct_apex_built_faces: &direct_apex_built_faces,
            warnings: &mut warnings,
            corner_patches: &mut corner_patches,
            repair_report: &mut repair_report,
        })
    });
    pipeline_timing.post_pipeline_ms = ms;
    repair_report
        .apex_diagnostics
        .push(pipeline_timing.diagnostic(
            edge_geometries.len(),
            corner_patches.len(),
            warnings.len(),
            direct_apex_built_faces.len(),
            late_wire_jobs.len(),
            support_seam_jobs.len(),
        ));
    Ok(ApplyBlendResult {
        brep: result,
        warnings,
        corner_patches,
        repair_report,
        unresolved_corners,
    })
}
