//! Topology construction — apply fillet/chamfer geometry to the BRep. 3-phase single-pass (A: perp splits, B: multi-chord face splits, C: kef merge + surface assign) then surface-fix, sliver-absorb (§13.1(2)), and §8.3 cap insertion post-passes.

mod apex_face_lift;
mod apex_vertex_coherence;
pub(crate) mod blend_insert;
pub(crate) mod blend_topo;
mod boundary_support_frames;
mod bridge_warnings;
pub(crate) mod cap_insert;
pub(crate) mod cap_outer_acceptance;
pub(crate) mod closed_fillet;
pub(crate) mod corner_patch_support_completion;
pub(crate) mod cyl_plane_intersect;
mod direct_apex_dispatch;
pub(crate) mod edge_reconcile;
mod exact_boundary_cap_trim_fit;
mod exact_boundary_cap_trim_restore;
mod exact_boundary_cap_trim_surface;
mod exact_boundary_export_repair;
mod gregory_pcurve_clamp;
/// Re-export of the late-pipeline clamp pass for use by `fillet_edges`
/// after `fix_stale_pcurves` has run.
pub fn clamp_bounded_pcurves_post_phase5(brep: &mut popkern_core::topology::BRep) -> usize {
    gregory_pcurve_clamp::clamp_bounded_pcurves(brep)
}

/// Restore exact-boundary cap trim authority after generic Phase 5 pcurve repair.
pub fn restore_exact_boundary_cap_pcurves_post_phase5(
    brep: &mut popkern_core::topology::BRep,
    repair_report: &mut crate::types::BlendRepairReport,
) {
    exact_boundary_cap_trim_restore::restore_exact_boundary_cap_pcurves(brep, repair_report);
}
#[cfg(test)]
mod apex_termination_spike;
pub(crate) mod helpers;
mod kp13_bi5_promote;
mod kp20_polar_diagnostics;
mod kp20_polar_dispatch;
mod kp_bi5_macro_apply;
mod kp_binonic_optimizer;
mod kp_biquintic_optimizer;
mod kp_exact_boundary_binonic_optimizer;
mod kp_lm_solver;
pub(crate) mod kp_polar_bi6_apply;
mod kp_polar_bi6_diagnostics;
pub(crate) mod kp_polar_bi6_promote;
mod kp_polar_deep_trace;
mod kp_polar_lm_optimizer;
mod kp_twist_optimizer;
pub(crate) mod merge;
mod mixed_radius_direct_apex;
mod multi_chord;
pub(crate) mod pcurve_attach;
mod pcurve_measured_repair;
mod pcurve_refit;
mod perpendicular_splits;
mod phase_split;
pub(crate) mod poke_4n;
pub(crate) mod poke_hex;
mod post_pipeline;
mod promotion_ladder;
mod route_decision;
mod sliver_absorb;
mod sliver_absorb_classifiers;
mod sphere_fit_edge_wire;
mod sphere_fit_edge_wire_cap_identity;
mod surface_driven_edge_repair;
pub(crate) mod surface_fix;
pub(crate) mod three_quad_apply;
mod three_quad_apply_edge_triple;
mod three_quad_boundary_fit;
pub(crate) mod three_quad_boundary_pcurve_diagnostics;
mod three_quad_kp_diagnostics;
pub(crate) mod three_quad_promotions;
mod three_quad_support_pcurve_completion;
pub(crate) mod topology_repair;

use crate::analysis::EdgeAnalysis;
use crate::error::{BlendFailure, BlendWarning};
use crate::types::VertexMode;
use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use route_decision::{ApexRouteDecisionInput, apex_route_decision};

/// Geometry computed for one fillet/chamfer edge.
#[derive(Debug)]
pub struct EdgeGeometry {
    /// The edge analysis.
    pub analysis: EdgeAnalysis,
    /// Fillet or chamfer surface.
    pub surface: SurfaceGeometry,
    /// Trim curve on face 1 (3D).
    pub trim_3d_1: NurbsCurve,
    /// Trim curve on face 2 (3D).
    pub trim_3d_2: NurbsCurve,
    /// Pcurve on face 1.
    pub pcurve_1: NurbsCurve,
    /// Pcurve on face 2.
    pub pcurve_2: NurbsCurve,
    /// Average fillet/chamfer radius (kept for legacy callers).
    pub effective_radius: f64,
    /// Radius at the edge's start vertex (constant ⇒ `effective_radius`; variable ⇒ `start`).
    pub radius_at_start: f64,
    /// Radius at the edge's end vertex.
    pub radius_at_end: f64,
    /// §18 telemetry: walker-path Newton iterations.
    pub iterations_used: Option<usize>,
    /// §18 telemetry: line-search rejections.
    pub step_rejections_used: Option<usize>,
    /// §8.3 start-endpoint cap (`None` when `cap_kind == Open`).
    pub start_cap: Option<crate::cap::CapGeometry>,
    /// §8.3 end-endpoint cap (`None` when `cap_kind == Open`).
    pub end_cap: Option<crate::cap::CapGeometry>,
    /// §7 continuity downgrade record — `Some((achieved, reason))` when the engine fell
    /// back; caller emits [`BlendWarning::ContinuityDowngraded`](crate::error::BlendWarning::ContinuityDowngraded).
    pub continuity_downgrade: Option<(crate::types::ContinuityLevel, String)>,
}

/// Result of applying blend geometry, including warnings.
pub struct ApplyBlendResult {
    /// The modified BRep.
    pub brep: BRep,
    /// Warnings generated during vertex mode handling.
    pub warnings: Vec<(EdgeHandle, BlendWarning)>,
    /// §5.3 per-corner patch reports — one entry per vertex where a blend patch was inserted.
    pub corner_patches: Vec<crate::types::CornerPatchInfo>,
    /// Structured diagnostics for repair and hardening actions.
    pub repair_report: crate::types::BlendRepairReport,
}

/// Apply computed fillet/chamfer geometry to a BRep. 3-phase: A perp splits, B multi-chord face splits, C kef merge + surface assign.
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
) -> Result<ApplyBlendResult, Vec<(EdgeHandle, BlendFailure)>> {
    let mut result = brep.clone();
    let mut failures: Vec<(EdgeHandle, BlendFailure)> = Vec::new();
    let mut warnings: Vec<(EdgeHandle, BlendWarning)> = Vec::new();
    let mut corner_patches: Vec<crate::types::CornerPatchInfo> = Vec::new();
    let mut repair_report = crate::types::BlendRepairReport::default();
    // Phase 6c-D Stage 4.4: per-vertex (boundary, sub_faces) for
    // a post-`sliver_absorb` wire pass. Populated by
    // `insert_vertex_blend` when sphere-fit substitution fired.
    let mut late_wire_jobs: Vec<(
        crate::blend::BlendBoundary,
        Vec<popkern_core::topology::FaceHandle>,
    )> = Vec::new();
    let mut support_seam_jobs: Vec<(
        crate::blend::BlendBoundary,
        Vec<popkern_core::topology::FaceHandle>,
    )> = Vec::new();
    // §6c-direct-apex dispatch.d — face handles authored by `try_direct_apex_construction`.
    // These faces are protected from `sliver_absorb` /
    // `cap_insert` / late re-wiring so the post-loop stages
    // don't re-author or absorb the direct-apex-built cap.
    let mut direct_apex_built_faces: Vec<popkern_core::topology::FaceHandle> = Vec::new();

    let mut split_state = perpendicular_splits::run_perpendicular_splits(
        brep,
        &mut result,
        edge_geometries,
        tolerance,
        use_direct_apex_construction,
        &mut failures,
    );
    phase_split::apply_face_trim_pipeline(
        brep,
        &mut result,
        edge_geometries,
        tolerance,
        &split_state.perp_data,
        &mut split_state.descendants,
        &mut failures,
    );

    if !failures.is_empty() && failures.len() == edge_geometries.len() {
        return Err(failures);
    }

    // Vertex blend handling. `BTreeMap` keyed on VertexHandle for
    // deterministic insertion order across runs.
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
                let c =
                    crate::vertex_classifier::classify_vertex_blend(brep, *v, &edges, tolerance)
                        .classification;
                (*v, c)
            })
            .collect();
        for (vertex, edge_indices) in &vertex_edges {
            if edge_indices.len() >= 2 {
                let pre = precomputed
                    .get(vertex)
                    .copied()
                    .unwrap_or(crate::vertex_classifier::VertexBlendClassification::SameSideConvex);

                // §6c-direct-apex dispatch.b/d — when use_direct_apex_construction is set,
                // try the direct-apex planner + builder for this apex.
                //   * On success, the BRep already has the direct-apex-built
                //     cap.  We collect the cap face handles into
                //     `direct_apex_built_faces` so the post-loop protected
                //     set (sliver_absorb / cap_insert /
                //     late_wire_jobs) leaves them alone, then skip
                //     the legacy insert_vertex_blend call entirely.
                //   * On direct-apex error (plan or build), fall back to
                //     legacy with a typed
                //     BlendWarning::DirectApexFellBackToLegacy recording the
                //     reason.  `construct_apex`'s snapshot has
                //     already rolled the BRep back to its pre-call
                //     state, so legacy can run cleanly.
                //
                // DirectApex's failure modes:
                //   * plan_apex_resolution returns ApexPlanError
                //     (e.g. SphereCylinderNoIntersection on extreme
                //     radii, SingularCapCentre on nearly-coplanar
                //     faces).
                //   * construct_apex returns ApexBuildError —
                //     today most commonly Gregory cap surface
                //     assignment (deferred to direct-apex cap-surface).
                let use_direct_apex_for_this_apex = use_direct_apex_construction
                    || mixed_radius_direct_apex::should_promote_mixed_radius_gregory_apex(
                        edge_indices,
                        edge_geometries,
                        use_rational_gregory,
                        use_kp_binonic_gregory,
                    );
                let direct_apex_outcome = if use_direct_apex_for_this_apex {
                    Some(direct_apex_dispatch::try_direct_apex_construction(
                        &mut result,
                        *vertex,
                        edge_indices,
                        edge_geometries,
                        &split_state.descendants,
                        tolerance,
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
                    repair_report
                        .apex_routes
                        .push(apex_route_decision(ApexRouteDecisionInput {
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
                        }));
                }

                match blend_insert::insert_vertex_blend(
                    &mut result,
                    *vertex,
                    edge_indices,
                    edge_geometries,
                    tolerance,
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
                        if direct_apex_fell_back_to_legacy {
                            let cap_faces = infos.iter().map(|info| info.face).collect::<Vec<_>>();
                            if let Some(diagnostic) =
                                cap_outer_acceptance::cap_support_continuity_diagnostic(
                                    &result, &cap_faces, 5.0e-3, 5.0e-2,
                                )
                            {
                                repair_report.apex_diagnostics.push(diagnostic);
                            }
                        }
                        for emit in [
                            bridge_warnings::emit_bridge_fallback_warnings,
                            bridge_warnings::emit_tangent_corner_warning,
                            bridge_warnings::emit_overlap_remediation_warning,
                            bridge_warnings::emit_corner_continuity_relaxed_warning,
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
                            if std::env::var_os("POPKERN_APEX_PATH_A_TRACE").is_some() {
                                eprintln!(
                                    "[apex-path-a/dispatch] vertex={vertex:?} \
                                     sphere_fit_substituted={} -> late_wire_jobs {}",
                                    b.sphere_fit_substituted,
                                    if b.sphere_fit_substituted {
                                        "PUSHED"
                                    } else {
                                        "skipped"
                                    },
                                );
                            }
                            let face_handles: Vec<_> = infos.iter().map(|cp| cp.face).collect();
                            if !face_handles.is_empty() {
                                support_seam_jobs.push((b.clone(), face_handles.clone()));
                            }
                            if b.sphere_fit_substituted {
                                late_wire_jobs.push((b, face_handles));
                            }
                        } else if std::env::var_os("POPKERN_APEX_PATH_A_TRACE").is_some() {
                            eprintln!(
                                "[apex-path-a/dispatch] vertex={vertex:?} \
                                 boundary_opt=None -> late_wire_jobs skipped"
                            );
                        }
                        corner_patches.extend(infos);
                    }
                    Err(_) => {
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
        }
    }

    post_pipeline::run_post_blend_pipeline(post_pipeline::PostBlendPipelineInputs {
        source_brep: brep,
        result: &mut result,
        edge_geometries,
        tolerance,
        late_wire_jobs: &late_wire_jobs,
        support_seam_jobs: &support_seam_jobs,
        direct_apex_built_faces: &direct_apex_built_faces,
        warnings: &mut warnings,
        corner_patches: &mut corner_patches,
        repair_report: &mut repair_report,
    });

    Ok(ApplyBlendResult {
        brep: result,
        warnings,
        corner_patches,
        repair_report,
    })
}
