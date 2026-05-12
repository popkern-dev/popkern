//! Topology construction — apply fillet/chamfer geometry to the BRep. 3-phase single-pass (A: perp splits, B: multi-chord face splits, C: kef merge + surface assign) then surface-fix, sliver-absorb (§13.1(2)), and §8.3 cap insertion post-passes.

mod a2_dispatch;
mod apex_face_lift;
mod apex_vertex_coherence;
pub(crate) mod blend_insert;
pub(crate) mod blend_topo;
mod boundary_support_frames;
mod bridge_warnings;
pub(crate) mod cap_insert;
pub(crate) mod closed_fillet;
pub(crate) mod cyl_plane_intersect;
pub(crate) mod edge_reconcile;
mod gregory_pcurve_clamp;

/// Re-export of the late-pipeline clamp pass for use by `fillet_edges`
/// after `fix_stale_pcurves` has run.
pub fn clamp_bounded_pcurves_post_phase5(brep: &mut popkern_core::topology::BRep) -> usize {
    gregory_pcurve_clamp::clamp_bounded_pcurves(brep)
}
#[cfg(test)]
mod apex_termination_spike;
pub(crate) mod helpers;
mod kp13_bi5_promote;
mod kp20_polar_dispatch;
mod kp_bi5_macro_apply;
mod kp_binonic_optimizer;
mod kp_biquintic_optimizer;
mod kp_lm_solver;
pub(crate) mod kp_polar_bi6_apply;
pub(crate) mod kp_polar_bi6_promote;
mod kp_polar_deep_trace;
mod kp_polar_lm_optimizer;
mod kp_twist_optimizer;
pub(crate) mod merge;
mod multi_chord;
pub(crate) mod pcurve_attach;
mod pcurve_refit;
mod perpendicular_splits;
mod phase_split;
pub(crate) mod poke_4n;
pub(crate) mod poke_hex;
mod post_pipeline;
mod promotion_ladder;
mod sliver_absorb;
mod sliver_absorb_classifiers;
mod sphere_fit_edge_wire;
pub(crate) mod surface_fix;
pub(crate) mod three_quad_apply;
mod three_quad_apply_edge_triple;
pub(crate) mod three_quad_promotions;
pub(crate) mod topology_repair;

use crate::analysis::EdgeAnalysis;
use crate::error::{BlendFailure, BlendWarning};
use crate::types::VertexMode;
use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

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
    use_a2_apex_construction: bool,
) -> Result<ApplyBlendResult, Vec<(EdgeHandle, BlendFailure)>> {
    let mut result = brep.clone();
    let mut failures: Vec<(EdgeHandle, BlendFailure)> = Vec::new();
    let mut warnings: Vec<(EdgeHandle, BlendWarning)> = Vec::new();
    let mut corner_patches: Vec<crate::types::CornerPatchInfo> = Vec::new();
    // Phase 6c-D Stage 4.4: per-vertex (boundary, sub_faces) for
    // a post-`sliver_absorb` wire pass. Populated by
    // `insert_vertex_blend` when sphere-fit substitution fired.
    let mut late_wire_jobs: Vec<(
        crate::blend::BlendBoundary,
        Vec<popkern_core::topology::FaceHandle>,
    )> = Vec::new();
    // §6c-A2.3.d — face handles authored by `try_a2_apex_construction`.
    // These faces are protected from `sliver_absorb` /
    // `cap_insert` / late re-wiring so the post-loop stages
    // don't re-author or absorb the A2-built cap.
    let mut a2_built_faces: Vec<popkern_core::topology::FaceHandle> = Vec::new();

    let mut split_state = perpendicular_splits::run_perpendicular_splits(
        brep,
        &mut result,
        edge_geometries,
        tolerance,
        use_a2_apex_construction,
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

                // §6c-A2.3.b/d — when use_a2_apex_construction is set,
                // try the A2 planner + builder for this apex.
                //   * On success, the BRep already has the A2-built
                //     cap.  We collect the cap face handles into
                //     `a2_built_faces` so the post-loop protected
                //     set (sliver_absorb / cap_insert /
                //     late_wire_jobs) leaves them alone, then skip
                //     the legacy insert_vertex_blend call entirely.
                //   * On A2 error (plan or build), fall back to
                //     legacy with a typed
                //     BlendWarning::A2FellBackToLegacy recording the
                //     reason.  `construct_apex`'s snapshot has
                //     already rolled the BRep back to its pre-call
                //     state, so legacy can run cleanly.
                //
                // A2's failure modes:
                //   * plan_apex_resolution returns ApexPlanError
                //     (e.g. SphereCylinderNoIntersection on extreme
                //     radii, SingularCapCentre on nearly-coplanar
                //     faces).
                //   * construct_apex returns ApexBuildError —
                //     today most commonly Gregory cap surface
                //     assignment (deferred to A2.4).
                let a2_outcome = if use_a2_apex_construction {
                    Some(a2_dispatch::try_a2_apex_construction(
                        &mut result,
                        *vertex,
                        edge_indices,
                        edge_geometries,
                        &split_state.descendants,
                        tolerance,
                    ))
                } else {
                    None // flag not set — legacy path runs unchanged
                };
                let mut a2_succeeded = false;
                if let Some(outcome) = a2_outcome {
                    match outcome {
                        Ok(success) => {
                            a2_built_faces.extend(success.built_faces);
                            corner_patches.extend(success.corner_patches);
                            a2_succeeded = true;
                        }
                        Err(reason) => {
                            let warn_edge = edge_geometries[edge_indices[0]].analysis.edge;
                            warnings.push((
                                warn_edge,
                                BlendWarning::A2FellBackToLegacy {
                                    vertex: *vertex,
                                    reason,
                                },
                            ));
                            // Fall through to legacy below.
                        }
                    }
                }
                if a2_succeeded {
                    continue;
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
                    Ok((infos, boundary_opt)) => {
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
                        // Phase 6c-D Stage 4.4: capture the blend
                        // boundary + per-vertex sub-face handles so a
                        // post-`sliver_absorb` wire pass can
                        // re-write substituted boundary curves to
                        // any cap-loop edges that were *added* by
                        // the absorb step. The wire originally runs
                        // mid-pipeline (inside three_quad_apply)
                        // when the cap face's outer loop is the
                        // raw poke output; sliver_absorb expands
                        // that loop with absorbed Plane fragments,
                        // leaving 8/12 cap-vs-cylinder edges on
                        // pyramid4_mixed unwired.
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
                            if b.sphere_fit_substituted {
                                let face_handles: Vec<_> = infos.iter().map(|cp| cp.face).collect();
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
        a2_built_faces: &a2_built_faces,
        warnings: &mut warnings,
        corner_patches: &mut corner_patches,
    });

    Ok(ApplyBlendResult {
        brep: result,
        warnings,
        corner_patches,
    })
}
