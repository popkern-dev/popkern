//! Topology construction — apply fillet/chamfer geometry to the BRep.
mod apex_face_lift;
mod apply_types;
pub use apply_types::{ApplyBlendResult, EdgeGeometry};
mod apex_vertex_coherence;
pub(crate) mod blend_insert;
pub(crate) mod blend_seam_corner;
mod blend_seam_corner_geometry;
mod blend_seam_corner_march;
mod blend_seam_corner_repair;
mod blend_seam_corner_swept;
mod blend_seam_corner_swept_exit;
pub(crate) mod blend_topo;
pub(crate) mod blend_topo_curved;
mod boundary_support_frames;
mod bridge_warnings;
pub(crate) mod cap_insert;
pub(crate) mod cap_outer_acceptance;
pub(crate) mod closed_fillet;
pub(crate) mod corner_arc_boundaries;
pub(crate) mod corner_cap_geometry;
mod corner_chord_settle;
pub(crate) mod corner_patch_support_completion;
mod corner_wall_carve;
pub(crate) mod cyl_plane_intersect;
mod direct_apex_dispatch;
pub(crate) mod edge_reconcile;
mod end_runoff;
mod end_runoff_report;
mod exact_boundary_cap_trim_fit;
mod exact_boundary_cap_trim_restore;
mod exact_boundary_cap_trim_surface;
mod exact_boundary_export_repair;
pub(crate) mod final_declared_accuracy;
mod gregory_pcurve_clamp;
pub(crate) mod junction_corner_solve;
pub(crate) mod junction_remnant_absorb;
mod junction_surface_join;
mod multi_chord_anchor_transplant;
pub(crate) mod multi_chord_region_seeding;
pub(crate) mod periodic_branch_finalize;
mod seam_self_fit;
mod setback_commit_gate;
pub(crate) use post_pipeline::determinism_trace_mod::trace_trim_endpoints as post_pipeline_trace_trim_endpoints;
pub(crate) mod pns_setback_candidate;
pub(crate) mod rim_chain_chamfer;
pub(crate) mod setback_vertex_topology;
pub(crate) use exact_boundary_cap_trim_restore::restore_exact_boundary_cap_pcurves as restore_exact_boundary_cap_pcurves_post_phase5;
pub(crate) use gregory_pcurve_clamp::clamp_bounded_pcurves as clamp_bounded_pcurves_post_phase5;
mod direct_apex_promotion;
pub(crate) mod helpers;
mod helpers_pcurve;
mod hole_rim_terminal_disc;
mod kp13_bi5_promote;
pub mod kp20_paper_border_extract;
mod kp20_polar_bi6_split;
mod kp20_polar_comparison;
mod kp20_polar_diagnostics;
mod kp20_polar_dispatch;
mod kp20_polar_expanded_support_binding;
mod kp20_polar_expanded_support_diagnostics;
mod kp20_polar_expanded_support_lane_shape;
mod kp20_polar_expanded_support_topology;
mod kp20_polar_expanded_support_topology_poke;
mod kp20_polar_fold_summary;
mod kp20_polar_multipatch_apply;
mod kp20_polar_multipatch_evaluate;
mod kp20_polar_multipatch_exact_guided;
mod kp20_polar_multipatch_inner;
mod kp20_polar_multipatch_inner_shoulder;
mod kp20_polar_multipatch_interface;
mod kp20_polar_multipatch_inventory;
mod kp20_polar_multipatch_labels;
mod kp20_polar_multipatch_strip;
mod kp20_polar_multipatch_strip_optimize;
mod kp20_polar_multipatch_topology_split;
mod kp20_polar_multipatch_validate;
mod kp20_polar_support_cell_inventory;
pub(crate) mod kp20_polar_support_cell_shared_trim_quality;
mod kp20_polar_support_image;
mod kp20_polar_support_image_ownership;
mod kp20_polar_support_image_params;
mod kp20_polar_support_image_separated;
mod kp20_polar_support_lift;
mod kp20_polar_support_macro_construction;
mod kp20_polar_support_reparameterization;
mod kp20_polar_support_strip_reauthor;
mod kp20_polar_support_surface;
mod kp20_polar_support_surface_split;
mod kp20_polar_support_transition_strip;
mod kp20_polar_support_two_lane;
mod kp_bi5_macro_apply;
mod kp_binonic_optimizer;
mod kp_biquintic_optimizer;
mod kp_exact_boundary_binonic_optimizer;
mod kp_lm_solver;
pub(crate) mod kp_polar_bi6_apply;
mod kp_polar_bi6_diagnostics;
pub(crate) mod kp_polar_bi6_promote;
mod kp_polar_deep_fill;
mod kp_polar_deep_trace;
mod kp_polar_lm_optimizer;
mod kp_polar_weighted_f5;
mod kp_twist_optimizer;
pub(crate) mod merge;
mod multi_chord;
mod multi_chord_crossings;
mod multi_chord_region;
mod multi_chord_region_separators;
mod multi_chord_span_geometry;
pub(crate) mod parallel_seam_weld;
pub(crate) mod pcurve_attach;
mod pcurve_attach_compute;
mod pcurve_measured_repair;
mod pcurve_refit;
mod perpendicular_splits;
mod phase_split;
mod pipeline_timing;
pub(crate) mod poke_4n;
pub(crate) mod poke_hex;
pub(crate) mod post_pipeline;
mod promotion_ladder;
mod route_decision;
mod sliver_absorb;
mod sliver_absorb_classifiers;
mod sphere_fit_edge_wire;
mod sphere_fit_edge_wire_cap_identity;
pub(crate) mod stale_vertex_settle;
mod strip_end_arcs;
pub(crate) mod strip_end_surface;
pub(crate) mod strip_orientation;
pub(crate) mod strip_pocket_expel;
pub(crate) mod strip_rail_coparameterization;
mod surface_driven_edge_repair;
pub(crate) mod surface_fix;
mod terminal_arc_author;
mod terminal_arc_fit_trim;
pub(crate) mod terminal_arc_trim;
mod terminal_chord_trim;
mod terminal_cone_trim;
mod terminal_conic_fit_trim;
mod terminal_corner_recognition;
mod terminal_obstructed_crossing;
mod terminal_obstructed_end;
pub(crate) mod three_quad_apply;
mod three_quad_apply_edge_triple;
mod three_quad_balanced_retry;
mod three_quad_boundary_fit;
pub(crate) mod three_quad_boundary_pcurve_diagnostics;
mod three_quad_expanded_kp20_apply;
mod three_quad_kp_diagnostics;
pub(crate) mod three_quad_promotions;
mod three_quad_spoke_seam_curves;
mod three_quad_support_pcurve_completion;
pub(crate) mod topology_repair;
mod torus_end_section_caps;
use crate::blend_warnings::trace_smooth_blend_mitred_fallback as trace_mitred;
use crate::error::{BlendFailure, BlendWarning};
use crate::types::VertexMode;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use route_decision::{ApexRouteDecisionInput, apex_route_decision};

/// Apply computed fillet/chamfer geometry to a BRep. 3-phase: A perp splits, B multi-chord face splits, C kef merge + surface assign. `excluded_direct_apex_vertices` = PK-250 retry exclusions.
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
) -> Result<ApplyBlendResult, Vec<(EdgeHandle, BlendFailure)>> {
    let mut result = brep.clone();
    if post_pipeline::determinism_trace_enabled() {
        post_pipeline::determinism_trace("pre_apply_input_brep", brep);
        for (i, eg) in edge_geometries.iter().enumerate() {
            eprintln!(
                "dtrace: pre_apply eg[{i}] edge={:?} f1={:?} f2={:?} surf={:016x} trim1={:016x} trim2={:016x}",
                eg.analysis.edge,
                eg.analysis.face_1,
                eg.analysis.face_2,
                post_pipeline::determinism_value_hash(&eg.surface),
                post_pipeline::determinism_value_hash(&eg.trim_3d_1),
                post_pipeline::determinism_value_hash(&eg.trim_3d_2),
            );
        }
    }
    let mut failures: Vec<(EdgeHandle, BlendFailure)> = Vec::new();
    let mut warnings: Vec<(EdgeHandle, BlendWarning)> = Vec::new();
    let mut corner_patches: Vec<crate::types::CornerPatchInfo> = Vec::new();
    let mut repair_report = crate::types::BlendRepairReport::default();
    // Phase 6c-D Stage 4.4: per-vertex (boundary, sub_faces) for a
    // post-`sliver_absorb` wire pass (sphere-fit substitution only).
    let mut late_wire_jobs: Vec<(
        crate::blend::BlendBoundary,
        Vec<popkern_core::topology::FaceHandle>,
    )> = Vec::new();
    let mut support_seam_jobs: Vec<(
        crate::blend::BlendBoundary,
        Vec<popkern_core::topology::FaceHandle>,
    )> = Vec::new();
    // §6c-direct-apex dispatch.d — faces authored by `try_direct_apex_construction`,
    // protected so sliver_absorb/cap_insert/late re-wiring don't re-author the cap.
    let mut direct_apex_built_faces: Vec<popkern_core::topology::FaceHandle> = Vec::new();
    let mut pipeline_timing = pipeline_timing::ApplyBlendPipelineTiming::default();
    // PK-127: promotion decided pre-Phase-A; anchors + dispatch share it.
    let mut promoted_apexes = direct_apex_promotion::promoted_apexes_for_mode(
        brep,
        vertex_mode,
        edge_geometries,
        direct_apex_promotion::PromotionFlags {
            use_rational_gregory,
            use_kp_binonic_gregory,
            explicit_legacy_surface_family: use_biquintic_full_gregory
                || use_biseptic_gregory
                || use_two_stripe_gregory
                || use_kp_binonic_gregory
                || use_kp_polar_gregory,
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
            use_direct_apex_construction,
            &promoted_apexes,
            &mut failures,
        )
    });
    pipeline_timing.perpendicular_split_ms = ms;
    direct_apex_built_faces.extend(split_state.closed_fillet_protected.iter().copied());
    post_pipeline::determinism_trace("post_perp_splits", &result);
    let (trim_outcome, ms) = pipeline_timing::timed(|| {
        phase_split::apply_face_trim_pipeline(
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
        return Err(failures);
    }

    // Vertex blend handling. `BTreeMap` keyed on VertexHandle for
    // deterministic insertion order across runs.
    let (_, ms) = pipeline_timing::timed(|| {
        if *vertex_mode == VertexMode::Smooth && edge_geometries.len() > 1 {
            let mut vertex_edges: std::collections::BTreeMap<VertexHandle, Vec<usize>> =
                std::collections::BTreeMap::new();
            for (idx, eg) in edge_geometries.iter().enumerate() {
                vertex_edges
                    .entry(eg.analysis.vertex_start)
                    .or_default()
                    .push(idx);
                vertex_edges
                    .entry(eg.analysis.vertex_end)
                    .or_default()
                    .push(idx);
            }
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
            for (vertex, edge_indices) in &vertex_edges {
                if edge_indices.len() >= 2 {
                    let pre = precomputed.get(vertex).copied().unwrap_or(
                        crate::vertex_classifier::VertexBlendClassification::SameSideConvex,
                    );
                    // PK-192: N=3 corners are claimed by the direct-apex
                    // route (PK-099), so they never reach the PnS
                    // dispatch; releasing them lets PnS bid first.
                    // CAUTION (measured 2026-07-28): when PnS then
                    // declines, the blend does NOT return here — it
                    // ships a creased legacy-Gregory cap (cap-cap g0
                    // 2.7e-2..1.9e-1, cos_min ~0.002, vs Coons ~1e-15 /
                    // 0.9999 with the switch off). Tolerable only
                    // because this is a recorded A/B probe switch,
                    // never a production default — PK-192 decided
                    // AGAINST releasing N=3 (pns_low_valence_gates.rs).
                    let released_for_pns = crate::blend::pns_promote::pns_first_at_low_valence()
                        && edge_indices.len() == 3;
                    // PK-258: the serving-retry exclusion binds the per-apex
                    // flag directly — pruning promotion alone is a no-op when
                    // the global flag is on (measured: byte-identical rerun).
                    let use_direct_apex_for_this_apex = !released_for_pns
                        && !excluded_direct_apex_vertices.contains(vertex)
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
                    let mut direct_apex_fell_back_to_legacy = false;
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
                                direct_apex_fell_back_to_legacy = true;
                                let warn_edge = edge_geometries[edge_indices[0]].analysis.edge;
                                let reason = error.reason;
                                repair_report.apex_diagnostics.extend(error.diagnostics);
                                repair_report
                                .apex_routes
                                .push(apex_route_decision(ApexRouteDecisionInput {
                                vertex: *vertex,
                                edge_indices,
                                edge_geometries,
                                route:
                                    crate::types::ApexConstructionRoute::DirectApexFallbackToLegacy,
                                built_faces: Vec::new(),
                                corner_patch_count: 0,
                                exact_boundary_sector_count: 0,
                                exact_boundary_prepared_curve_count: 0,
                                exact_boundary_inner_spoke_count: 0,
                                exact_boundary_ribbon_count: 0,
                                exact_boundary_solved_surface_count: 0,
                                exact_boundary_pcurve_pair_count: 0,
                                exact_boundary_performance: Default::default(),
                                exact_boundary_report:
                                    crate::types::ExactBoundaryCapReport::not_applicable(
                                        "direct_apex_fallback_to_legacy",
                                    ),
                                fallback_reason: Some(reason.clone()),
                            }));
                                warnings.push((
                                    warn_edge,
                                    BlendWarning::DirectApexFellBackToLegacy {
                                        vertex: *vertex,
                                        reason,
                                    },
                                ));
                                // Fall through to legacy below.
                            }
                        }
                    }
                    if direct_apex_succeeded {
                        continue;
                    }
                    if !use_direct_apex_for_this_apex {
                        repair_report.apex_routes.push(apex_route_decision(
                            ApexRouteDecisionInput {
                                vertex: *vertex,
                                edge_indices,
                                edge_geometries,
                                route: crate::types::ApexConstructionRoute::Legacy,
                                built_faces: Vec::new(),
                                corner_patch_count: 0,
                                exact_boundary_sector_count: 0,
                                exact_boundary_prepared_curve_count: 0,
                                exact_boundary_inner_spoke_count: 0,
                                exact_boundary_ribbon_count: 0,
                                exact_boundary_solved_surface_count: 0,
                                exact_boundary_pcurve_pair_count: 0,
                                exact_boundary_performance: Default::default(),
                                exact_boundary_report:
                                    crate::types::ExactBoundaryCapReport::not_applicable("legacy"),
                                fallback_reason: None,
                            },
                        ));
                    }

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
                    ) {
                        Ok((infos, boundary_opt, diagnostics)) => {
                            repair_report.apex_diagnostics.extend(diagnostics);
                            route_decision::record_legacy_serve_on_route(
                                &result,
                                &mut repair_report,
                                *vertex,
                                infos.iter().map(|info| info.face).collect(),
                                direct_apex_fell_back_to_legacy,
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
                }
                post_pipeline::determinism_trace(&format!("post_vertex_{vertex:?}"), &result);
            }
        }
    });
    pipeline_timing.vertex_blend_ms = ms;
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
    })
}
