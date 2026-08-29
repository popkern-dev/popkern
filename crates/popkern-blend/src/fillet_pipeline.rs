//! The `fillet_edges` pipeline body — validation through final
//! repair, extracted from `lib.rs` (which keeps the public entry
//! points and re-exports).
//!
//! Phase structure and the repair-pass ordering are documented inline;
//! every pass here runs on the post-apply body, and several are
//! order-sensitive (see the PK-148/PK-141/PK-139 notes below).

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::apex_construction::repeated_apex_partition_topology::repair_stale_shell_anchors;
use crate::error::{BlendError, BlendFailure, BlendWarning};
use crate::geometry_compute::{build_overlap_inputs, compute_fillet_geometry};
use crate::naming::assign_fillet_names;
use crate::types::{BlendResult, FilletOptions, FilletSpec};
use crate::{
    analysis, cancellation, collapsed_skip, construct, continuity, cusp_subdivide,
    interference_seam, overlap, progress, residual_scaffolding, surface_boundary_evidence,
    validate,
};

/// One pipeline attempt with the given direct-apex promotion
/// exclusions (empty on the first pass). The serving-boundary retry
/// wrapper lives in [`crate::fillet_serving_retry`].
#[allow(clippy::too_many_arguments)]
pub(crate) fn fillet_edges_attempt(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: &mut Option<&mut TrackingContext>,
    cancel: Option<&cancellation::CancellationToken>,
    reporter: &mut Option<&mut dyn progress::ProgressReporter>,
    excluded_direct_apex_vertices: &std::collections::BTreeSet<VertexHandle>,
) -> Result<BlendResult, BlendError> {
    macro_rules! emit {
        ($phase:expr, $start:expr, $count:expr) => {
            if let Some(r) = reporter.as_deref_mut() {
                r.on_phase(progress::ProgressEvent {
                    phase: $phase,
                    duration: $start.elapsed(),
                    entity_count: $count,
                });
            }
        };
    }
    macro_rules! check_cancel {
        () => {
            if let Some(tok) = cancel
                && tok.is_cancelled()
            {
                return Err(BlendError::Cancelled);
            }
        };
    }
    // §7 request-implies-capability: a G2 request dispatches the
    // G2-capable construction ladder (rational Gregory rungs + the
    // KP setback/polar cap, whose commit postcondition now gates the
    // honest transverse curvature residual — PK-187). Without this a
    // caller asking for G2 through the default options silently got
    // the flag-gated ladder skipped: the strips upgraded to the G2
    // sweep but every N ≥ 5 corner shipped a G1-at-best cap. The
    // ladder flags stay opt-in under G1 requests, where the analytic
    // paths are the commercial default.
    // PK-189 profile validation and resolution. Conic rho is a
    // strict open-interval parameter; conic sections are tangent-
    // continuous by construction, so a simultaneous G2 request is a
    // contradiction the caller must resolve, not a silent downgrade.
    // BlendProfile::CurvatureContinuous is the G2 sweep spelled as a
    // profile — it implies a G2 continuity request below.
    if let crate::types::BlendProfile::Conic { rho } = options.profile {
        if !(rho.is_finite() && 0.0 < rho && rho < 1.0) {
            return Err(BlendError::InvalidParameter {
                parameter: "options.profile.rho".to_string(),
                value: format!("{rho}"),
                reason: "conic shoulder factor must lie strictly between 0 and 1 \
                         (0.5 = parabolic; the circular value is \
                         cos(θ/2)/(1+cos(θ/2)) for section opening angle θ)"
                    .to_string(),
            });
        }
        if options.continuity == crate::types::ContinuityLevel::G2 {
            return Err(BlendError::InvalidParameter {
                parameter: "options.profile".to_string(),
                value: format!("Conic {{ rho: {rho} }}"),
                reason: "conic profiles are tangent-continuous (G1) by construction; \
                         curvature-continuous blends use \
                         BlendProfile::CurvatureContinuous or ContinuityLevel::G2 \
                         with the circular profile"
                    .to_string(),
            });
        }
    }
    let g2_ladder_options;
    let wants_g2_construction = options.continuity == crate::types::ContinuityLevel::G2
        || options.profile == crate::types::BlendProfile::CurvatureContinuous;
    let options: &FilletOptions = if wants_g2_construction {
        g2_ladder_options = FilletOptions {
            continuity: crate::types::ContinuityLevel::G2,
            use_rational_gregory: true,
            use_deg5_gregory: true,
            use_biquintic_full_gregory: true,
            use_biseptic_gregory: true,
            use_two_stripe_gregory: true,
            use_kp_binonic_gregory: true,
            use_kp_polar_gregory: true,
            ..options.clone()
        };
        &g2_ladder_options
    } else {
        options
    };
    // PHASE 0: Validate.
    let t0 = popkern_core::time::Stopwatch::start();
    validate::validate_fillet_specs(brep, specs)?;
    emit!(progress::FilletPhase::Validate, t0, specs.len());
    check_cancel!();
    // Canonicalize radius laws AFTER validation (errors report the
    // user's spelling); see [`FilletSpec::canonicalized_all`].
    let specs: Vec<FilletSpec> = FilletSpec::canonicalized_all(specs);
    let specs = specs.as_slice();
    // §13.4 soft-supports: shadow-relax tolerance for low-precision input.
    let (tolerance, soft_warn) = options.relax_tolerance(tolerance);
    let mut prepass_warnings: Vec<BlendWarning> = soft_warn.into_iter().collect();
    // PHASE 0b: §13.1(5) collapsed-edge soft-skip + §13.1(1) cusp subdivide.
    let (skipped, collapsed_warnings) =
        collapsed_skip::skip_collapsed_edges(brep, specs, &tolerance);
    prepass_warnings.extend(collapsed_warnings);
    let (healed_brep, healed_specs) =
        cusp_subdivide::subdivide_cusp_edges(brep, &skipped, &tolerance);
    let brep = &healed_brep;
    let specs: &[FilletSpec] = &healed_specs;
    // Capture pre-operation faces for naming.
    let pre_faces: std::collections::HashSet<FaceHandle> = brep.faces.keys().collect();
    // PHASE 1: Analyze edges.
    let t1 = popkern_core::time::Stopwatch::start();
    let edges: Vec<EdgeHandle> = specs.iter().map(|s| s.edge).collect();
    let analyses = analysis::analyze_edges(brep, &edges, options.propagate_tangent)
        .map_err(BlendError::NoValidEdges)?;
    emit!(progress::FilletPhase::Analyze, t1, analyses.len());
    check_cancel!();
    // PHASE 2: Compute geometry for each edge.
    let t2 = popkern_core::time::Stopwatch::start();
    let mut edge_geometries = Vec::new();
    let mut failures = Vec::new();
    for (i, spec) in specs.iter().enumerate() {
        if i >= analyses.len() {
            continue;
        }
        // Per-edge cancellation point — lets the user kill a long
        // NURBS surface intersection mid-stride.
        check_cancel!();
        // §18 telemetry: per-edge Instant bracketing + reporter event.
        let edge_start = popkern_core::time::Stopwatch::start();
        let compute_result = compute_fillet_geometry(
            &analyses[i],
            &spec.radius,
            brep,
            &tolerance,
            options.continuity,
            options.cap_kind,
            options.profile,
        );
        if let Some(r) = reporter.as_deref_mut() {
            let (succeeded, newton_iterations, step_rejections) = match &compute_result {
                Ok(g) => (true, g.iterations_used, g.step_rejections_used),
                Err(_) => (false, None, None),
            };
            r.on_edge(progress::EdgeTelemetry {
                edge: spec.edge,
                duration: edge_start.elapsed(),
                newton_iterations,
                step_rejections,
                succeeded,
            });
        }
        match compute_result {
            Ok(geom) => edge_geometries.push(geom),
            Err(failure) => failures.push((spec.edge, failure)),
        }
    }
    construct::post_pipeline_trace_trim_endpoints(&edge_geometries);
    emit!(
        progress::FilletPhase::ComputeGeometry,
        t2,
        edge_geometries.len()
    );
    check_cancel!();
    // Transactional: any per-edge failure fails the whole op.
    if !failures.is_empty() {
        return Err(BlendError::PartialFailure {
            failed: failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures,
        });
    }
    // PHASE 3: Detect overlaps (warnings only — don't block the operation).
    let t3 = popkern_core::time::Stopwatch::start();
    let mut overlap_warnings = Vec::new();
    let (overlap_trims, _face_trim_map) = build_overlap_inputs(&edge_geometries, brep);
    let overlaps = overlap::detect_overlaps(&overlap_trims, &tolerance);
    for ov in &overlaps {
        overlap_warnings.push(BlendWarning::FilletOverlap {
            fillet_a: ov.fillet_a,
            fillet_b: ov.fillet_b,
        });
    }
    // Strip interference: complete parallel-cylinder collisions by
    // seam trimming (PK-159); every other configuration follows the
    // caller's interference policy — transactional decline (PK-140,
    // default) or drop-the-pair partial service (PK-185).
    interference_seam::complete_or_decline_interference(
        &mut edge_geometries,
        brep,
        &tolerance,
        options.interference_policy,
        &mut overlap_warnings,
    )?;
    emit!(
        progress::FilletPhase::DetectOverlaps,
        t3,
        overlap_warnings.len()
    );
    check_cancel!();
    if edge_geometries.is_empty() {
        return Err(BlendError::PartialFailure {
            failed: failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures,
        });
    }

    // PHASE 4: Apply topology.
    let t4 = popkern_core::time::Stopwatch::start();
    let apply_args = (&edge_geometries, &tolerance, &options.vertex_mode);
    let mut apply_result = construct::apply_blend_geometry(
        brep,
        apply_args.0,
        apply_args.1,
        apply_args.2,
        options.use_rational_gregory,
        options.use_deg5_gregory,
        options.use_biquintic_full_gregory,
        options.use_biseptic_gregory,
        options.use_two_stripe_gregory,
        options.use_kp_binonic_gregory,
        options.use_kp_polar_gregory,
        options.use_direct_apex_construction,
        options.continuity,
        excluded_direct_apex_vertices,
        cancel,
        // PK-276: the apply phase ticks the reporter once per
        // vertex-blend corner (the dominant cost inside ApplyTopology).
        reporter.as_deref_mut(),
    )
    .map_err(|apply_error| apply_error.into_blend_error(specs.len(), brep))?;
    let stale_anchor_failure = |reason: String| BlendError::PartialFailure {
        failed: 1,
        total: specs.len(),
        partial_result: Box::new(brep.clone()),
        failures: vec![(
            fallback_edge_for_specs(specs, &edge_geometries),
            BlendFailure::InvalidTopology {
                reason,
                location: None,
            },
        )],
    };
    repair_stale_shell_anchors(&mut apply_result.brep)
        .map_err(|err| stale_anchor_failure(format!("repair stale shell anchors: {err}")))?;
    emit!(
        progress::FilletPhase::ApplyTopology,
        t4,
        edge_geometries.len()
    );
    check_cancel!();
    let t5 = popkern_core::time::Stopwatch::start();
    // Structured tail timing (PK-278 profile): the four passes below
    // carried ~2.9 s of the 5.9 s G2 cube before parallelization —
    // recorded as a diagnostic fact so tail regressions stay visible.
    let mut tail_timing = construct::pipeline_timing::RepairTailTiming::default();
    // Whole-model by measurement (PK-274 postmortem): apply-phase
    // repairs (edge reconcile chord rebuilds, vertex snaps) touch
    // faces arbitrarily far from the blend region, so no k-ring
    // scoping of this heal is sound — the blend-on-blend torus
    // fixture shipped a silent 7.7e-2-stale pcurve two rings out
    // under a scoped pass. The passes PK-274 kept scoped (per-corner
    // consolidation heals, boolean fill-only tails) carry their own
    // touched-face inventories; this pipeline-tail net stays global.
    let _ = apply_result.brep.recompute_all_pcurves();
    let ((), ms) = construct::pipeline_timing::timed_ms(|| apply_result.brep.fix_stale_pcurves());
    tail_timing.fix_stale_first_ms = ms;
    check_cancel!();
    let _ = construct::clamp_bounded_pcurves(&mut apply_result.brep);
    repair_stale_shell_anchors(&mut apply_result.brep).map_err(|err| {
        stale_anchor_failure(format!(
            "repair stale shell anchors after pcurve repair: {err}"
        ))
    })?;
    surface_boundary_evidence::audit_gregory_analytic_boundary_evidence(&mut apply_result.brep);
    construct::restore_exact_boundary_cap_pcurves(
        &mut apply_result.brep,
        &mut apply_result.repair_report,
    );
    apply_result
        .repair_report
        .edge_curve_endpoint_coherence_repairs +=
        apply_result.brep.separate_nonmanifold_edge_uses();
    apply_result
        .repair_report
        .edge_curve_endpoint_coherence_repairs += apply_result
        .brep
        .repair_edge_curve_endpoint_coherence(tolerance.linear().max(1.0e-8));
    construct::corner_patch_support_completion::complete_corner_patch_support_pcurves(
        &mut apply_result.brep,
        &apply_result.corner_patches,
        &mut apply_result.repair_report.apex_diagnostics,
    );
    construct::kp20_polar_support_cell_shared_trim_quality::append_final_support_cell_shared_trim_quality(
        &apply_result.brep,
        &apply_result.corner_patches,
        &mut apply_result.repair_report.apex_diagnostics,
    );
    // Final pcurve fidelity pass: intermediate repair stages can
    // (re)write low-fidelity pcurves after the first fix_stale run
    // (PK-095: straight 2-point pcurves on fillet-sweep trims, 5e-2
    // off, were the dominant tessellation seam crack). Whole-model —
    // see the PK-274 postmortem note on the first pass above.
    let ((), ms) = construct::pipeline_timing::timed_ms(|| apply_result.brep.fix_stale_pcurves());
    tail_timing.fix_stale_second_ms = ms;
    check_cancel!();
    // PK-141 strip cleanup: seam slits + foreign pockets (see
    // cleanup_strip_slits_and_pockets docs).
    if apply_result.repair_report.region_resolved_trims > 0 {
        construct::junction_remnant_absorb::cleanup_strip_slits_and_pockets(
            &mut apply_result.brep,
            &tolerance,
        );
    }
    if apply_result.repair_report.region_resolved_trims > 0 {
        let _ = construct::junction_remnant_absorb::absorb_junction_remnants(
            &mut apply_result.brep,
            &edge_geometries,
            &tolerance,
        );
    }
    // PK-159: weld duplicated mid-face parallel-fillet seams first.
    construct::parallel_seam_weld::weld_parallel_seams(brep, &mut apply_result.brep, &tolerance);
    // PK-148: settle stale chain-end corners after all pcurve authors.
    apply_result.repair_report.terminal_arc_trim =
        construct::terminal_arc_trim::settle_terminal_arc_stubs(
            brep,
            &mut apply_result.brep,
            &tolerance,
            construct::terminal_arc_trim::PlanarStubRecognition::ChartGated,
        );
    apply_result.repair_report.blend_seam_corner =
        construct::blend_seam_corner::settle_blend_seam_corners(
            brep,
            &mut apply_result.brep,
            &tolerance,
        );
    // PK-141 settle: AFTER corner settlements, region bodies only.
    if apply_result.repair_report.region_resolved_trims > 0 {
        let _ = construct::stale_vertex_settle::settle_stale_strip_vertices(
            &mut apply_result.brep,
            &tolerance,
        );
    }
    // PK-139: co-parameterize swept-strip rails whose pcurves diverge
    // from the rail curves at the same fraction (curved variable
    // rails; the imprecision-scaled stale gate cannot serve them).
    {
        let pre_faces: std::collections::HashSet<popkern_core::topology::FaceHandle> =
            brep.faces.keys().collect();
        let created: Vec<popkern_core::topology::FaceHandle> = apply_result
            .brep
            .faces
            .keys()
            .filter(|f| !pre_faces.contains(f))
            .collect();
        let (_, ms) = construct::pipeline_timing::timed_ms(|| {
            construct::strip_rail_coparameterization::coparameterize_strip_rails(
                &mut apply_result.brep,
                &created,
                &tolerance,
            )
        });
        tail_timing.coparameterize_ms = ms;
    }
    // PK-175 final bounded-UV clamp; PK-141 period-branch reconcile
    // on periodic supports; PK-177 declared-accuracy recording.
    let _ = construct::clamp_bounded_pcurves(&mut apply_result.brep);
    let _ = construct::periodic_branch_finalize::reconcile_all_periodic_pcurve_branches(
        &mut apply_result.brep,
    );
    emit!(
        progress::FilletPhase::RepairPcurves,
        t5,
        apply_result.brep.half_edges.len()
    );
    // Merge overlap warnings with topology warnings, plus §13.1 pre-pass warnings.
    let mut all_warnings = apply_result.warnings;
    let fallback_edge = specs
        .first()
        .map(|s| s.edge)
        .or_else(|| edge_geometries.first().map(|g| g.analysis.edge))
        .unwrap_or_default();
    for warn in overlap_warnings {
        all_warnings.push((fallback_edge, warn));
    }
    for warn in prepass_warnings {
        let edge_h = match &warn {
            BlendWarning::CollapsedEdgeSkipped { edge, .. } => *edge,
            _ => fallback_edge,
        };
        all_warnings.push((edge_h, warn));
    }
    continuity::emit_downgrade_warnings(&edge_geometries, options, &mut all_warnings);
    continuity::emit_corner_downgrade_warnings(
        &apply_result.corner_patches,
        options,
        fallback_edge,
        &mut all_warnings,
    );
    apply_result.repair_report.residual_scaffolding =
        residual_scaffolding::retire_unglued_construction_residual_scaffolding(
            &mut apply_result.brep,
        );
    let (_, ms) = construct::pipeline_timing::timed_ms(|| {
        construct::final_declared_accuracy::record_final_declared_accuracy(&mut apply_result.brep)
    });
    tail_timing.declared_accuracy_ms = ms;
    apply_result
        .repair_report
        .apex_diagnostics
        .push(tail_timing.diagnostic());
    if let Some(reason) = first_structural_validation_error(&apply_result.brep) {
        return Err(BlendError::PartialFailure {
            failed: 1,
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures: vec![(
                fallback_edge,
                BlendFailure::InvalidTopology {
                    reason,
                    location: None,
                },
            )],
        });
    }
    // Assign stable names when tracking is active.
    if let Some(ctx) = tracking.as_deref_mut() {
        assign_fillet_names(&mut apply_result.brep, ctx, &pre_faces);
    }
    // PK-211 displaced-pcurve + PK-294 loop-winding backstops: see
    // `fillet_pipeline_pcurve_backstop` for both defect classes.
    let winding_repairs =
        crate::fillet_pipeline_pcurve_backstop::final_pcurve_backstops(&mut apply_result.brep);
    apply_result
        .repair_report
        .edge_curve_endpoint_coherence_repairs += winding_repairs;
    // PK-244 serving-boundary flag reconcile: the orientation flag is
    // derived state (loop winding × surface geometry), so once every
    // pcurve author above has run, each measurable face's flag is set
    // to its honest value. This must be the LAST pass that runs after
    // pcurve finalization — an earlier run (post_pipeline) measured a
    // clean model whose Gregory-cap pcurves did not yet exist.
    apply_result.repair_report.orientation_flag_reconciliations +=
        apply_result.brep.reconcile_orientation_flags();
    // PK-220 serving-boundary mid-seam refresh: the commit-time cap
    // seam measurement runs through not-yet-final pcurves and
    // over-reads the pole half of each spoke; re-measure on the final
    // body and stamp `served_mid_seam_normal` so gates pin the truth.
    let _ = crate::apex_construction::exact_boundary_cap::cap_quality_postcondition::
        refresh_served_mid_seam_diagnostics(&apply_result.brep, &mut apply_result.repair_report);
    Ok(BlendResult {
        brep: apply_result.brep,
        edge_warnings: all_warnings,
        corner_patches: apply_result.corner_patches,
        repair_report: apply_result.repair_report,
    })
}
#[path = "fillet_pipeline_postconditions.rs"]
mod fillet_pipeline_postconditions;
use fillet_pipeline_postconditions::{fallback_edge_for_specs, first_structural_validation_error};
