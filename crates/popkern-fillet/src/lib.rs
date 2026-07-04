#![warn(missing_docs)]
//! popkern-fillet: Fillets and chamfers for the popkern BRep kernel.
//! Public entry points: [`fillet_edges`], [`chamfer_edges`], [`vertex_blend_3d`].
pub mod analysis;
pub mod apex_construction;
pub mod blend;
pub mod blend_diagnostics;
pub(crate) mod blend_warnings;
pub mod cancellation;
pub mod cap;
pub mod center_curve;
pub mod center_curve_helpers;
pub mod center_curve_walker;
pub mod chamfer_2d;
pub mod chamfer_2d_by_angle;
pub(crate) mod chamfer_2d_ray;
pub mod chamfer_curve;
pub(crate) mod collapsed_skip;
pub mod construct;
pub mod contact;
pub(crate) mod continuity;
pub(crate) mod corner_2d_helpers;

/// Test-only hooks for gated work-in-progress construction arms.
pub mod construct_test_hooks {
    /// PK-097 Phase A: run the PnS setback candidate (G1 cap on the
    /// setback border) and return its quality report.
    pub fn pns_setback_candidate(
        fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
        vertex_position: popkern_core::geometry::Point3,
    ) -> Result<crate::construct::pns_setback_candidate::PnsSetbackReport, String> {
        crate::construct::pns_setback_candidate::try_pns_setback_candidate(
            fillet_surfaces,
            vertex_position,
        )
        .map(|c| c.report)
    }

    /// Force COMMIT mode on (tests; production defaults to on).
    pub fn enable_setback_topology_commit() {
        crate::construct::setback_vertex_topology::SETBACK_TOPOLOGY_COMMIT
            .store(true, std::sync::atomic::Ordering::Relaxed);
        crate::construct::setback_vertex_topology::SETBACK_TOPOLOGY_COMMIT_SET
            .store(true, std::sync::atomic::Ordering::Relaxed);
    }

    /// Force COMMIT mode off (tests needing the legacy corner path).
    pub fn disable_setback_topology_commit() {
        crate::construct::setback_vertex_topology::SETBACK_TOPOLOGY_COMMIT
            .store(false, std::sync::atomic::Ordering::Relaxed);
        crate::construct::setback_vertex_topology::SETBACK_TOPOLOGY_COMMIT_SET
            .store(true, std::sync::atomic::Ordering::Relaxed);
    }

    /// Read the last WIP checkpoint report (integration tests).
    pub fn take_setback_topology_checkpoint() -> Option<String> {
        crate::construct::setback_vertex_topology::SETBACK_TOPOLOGY_LAST
            .lock()
            .unwrap_or_else(|p| p.into_inner())
            .take()
    }
}
pub mod curvature_comb;
pub(crate) mod cusp_subdivide;
pub(crate) mod degeneracy;
pub(crate) mod degeneracy_surface;
pub(crate) mod direct_apex_trace;
pub mod edges_2d;
pub mod error;
pub mod fillet_2d;
pub(crate) mod fillet_2d_plan;
pub(crate) mod fillet_2d_solvers;
pub mod fillet_surface;
pub mod fillet_surface_analytic;
#[cfg(test)]
#[path = "g2_acceptance_tests.rs"]
mod g2_acceptance_tests;
pub mod g2_cross_section;
pub(crate) mod geometry_compute;
pub(crate) mod geometry_compute_caps;
pub(crate) mod geometry_compute_cone;
pub(crate) mod geometry_compute_feasibility;
pub mod highlight_line;
pub mod history;
pub mod identify;
pub(crate) mod identify_helpers;
pub mod modify;
pub mod modify_inplace;
pub(crate) mod modify_inplace_closed_sweep;
pub(crate) mod modify_inplace_geometry;
pub(crate) mod modify_inplace_sweep;
pub(crate) mod naming;
pub mod offset;
pub mod overlap;
pub mod precheck;
pub mod preview;
pub mod progress;
pub mod quality;
pub mod reflection_line;
pub mod remediation;
pub mod sensitivity;
mod surface_boundary_evidence;
pub mod surface_validate;
pub mod trim;
pub mod types;

#[cfg(test)]
#[path = "types_tests.rs"]
mod types_tests;
pub mod unfillet;
pub mod unfillet_closed;
pub(crate) mod unfillet_closed_helpers;
pub(crate) mod unfillet_helpers;
pub(crate) mod unfillet_surgery;
pub mod validate;
pub mod vertex_blend;
pub(crate) mod vertex_blend_topology;
pub mod vertex_classifier;
pub mod walker;
pub(crate) mod walker_adaptive;
pub mod zebra;
pub use cancellation::CancellationToken;
pub use chamfer_2d::{chamfer_corner_2d, chamfer_corner_2d_by_angle, chamfer_corners_2d};
pub use edges_2d::*;
pub use error::{
    BlendError, BlendFailure, BlendWarning, Remediation, RemediationAction, RemediationContext,
};
pub use fillet_2d::{fillet_corner_2d, fillet_corners_2d};
pub use history::{
    ChamferOperationRecord, FilletOperationRecord, chamfer_recorded, fillet_recorded,
};
pub use identify::{
    BlendInfo, BlendKind, IdentificationSource, blend_neighbours, blend_radius,
    blend_support_faces, identify_all_blends, identify_blend,
};
pub use modify::modify_fillet_radius;
pub use precheck::{
    BlendQuality, ComputeQualityInput, FilletFeasibility, FilletIssue, FilletRadiusBounds,
    feasible_radius_bounds, fillet_quality, precheck_fillet,
};
pub use preview::{
    CornerPatchPreview, PreviewEdgeGeometry, PreviewGeometry, preview_chamfer,
    preview_corner_patches, preview_fillet, preview_fillet_with_continuity,
};
pub use progress::*;
pub use types::*;
pub use unfillet::{UnfilletPlan, compute_unfillet_plan, unfillet_face, unfillet_via_record};
pub use vertex_blend::{vertex_blend_3d, vertex_blend_3d_variable};

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::apex_construction::repeated_apex_partition_topology::repair_stale_shell_anchors;
use crate::construct::EdgeGeometry;

/// Fillet the specified edges of a BRep solid. See
/// [`BlendError`] for failure modes.
pub fn fillet_edges(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    fillet_edges_internal(brep, specs, options, tolerance, tracking, None, None)
}

/// Cancellable variant of [`fillet_edges`] — returns
/// [`BlendError::Cancelled`] if the token fires between phases.
#[rustfmt::skip]
pub fn fillet_edges_cancellable(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    cancel: &cancellation::CancellationToken,
) -> Result<BlendResult, BlendError> {
    fillet_edges_internal(brep, specs, options, tolerance, tracking, Some(cancel), None)
}

/// Progress-reporting variant of [`fillet_edges`]. See [`ProgressReporter`].
#[rustfmt::skip]
pub fn fillet_edges_with_progress(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    reporter: &mut dyn progress::ProgressReporter,
) -> Result<BlendResult, BlendError> {
    fillet_edges_internal(brep, specs, options, tolerance, tracking, None, Some(reporter))
}

/// Pipeline body for [`fillet_edges`] + cancellable/progress variants.
fn fillet_edges_internal(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    cancel: Option<&cancellation::CancellationToken>,
    mut reporter: Option<&mut dyn progress::ProgressReporter>,
) -> Result<BlendResult, BlendError> {
    use std::time::Instant;

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

    // PHASE 0: Validate.
    let t0 = Instant::now();
    validate::validate_fillet_specs(brep, specs)?;
    emit!(progress::FilletPhase::Validate, t0, specs.len());
    check_cancel!();

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
    let t1 = Instant::now();
    let edges: Vec<EdgeHandle> = specs.iter().map(|s| s.edge).collect();
    let analyses = analysis::analyze_edges(brep, &edges, options.propagate_tangent)
        .map_err(BlendError::NoValidEdges)?;
    emit!(progress::FilletPhase::Analyze, t1, analyses.len());
    check_cancel!();

    // PHASE 2: Compute geometry for each edge.
    let t2 = Instant::now();
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
        let edge_start = Instant::now();
        let compute_result = compute_fillet_geometry(
            &analyses[i],
            &spec.radius,
            brep,
            &tolerance,
            options.continuity,
            options.cap_kind,
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
    emit!(
        progress::FilletPhase::ComputeGeometry,
        t2,
        edge_geometries.len()
    );
    check_cancel!();

    if edge_geometries.is_empty() && !failures.is_empty() {
        return Err(BlendError::PartialFailure {
            failed: failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures,
        });
    }

    // PHASE 3: Detect overlaps (warnings only — don't block the operation).
    let t3 = Instant::now();
    let mut overlap_warnings = Vec::new();
    let (overlap_trims, _face_trim_map) = build_overlap_inputs(&edge_geometries, brep);

    let overlaps = overlap::detect_overlaps(&overlap_trims, &tolerance);
    for ov in &overlaps {
        overlap_warnings.push(BlendWarning::FilletOverlap {
            fillet_a: ov.fillet_a,
            fillet_b: ov.fillet_b,
        });
    }
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
    let t4 = Instant::now();
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
    )
    .map_err(|edge_failures| BlendError::PartialFailure {
        failed: edge_failures.len(),
        total: specs.len(),
        partial_result: Box::new(brep.clone()),
        failures: edge_failures,
    })?;
    repair_stale_shell_anchors(&mut apply_result.brep).map_err(|err| {
        BlendError::PartialFailure {
            failed: 1,
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures: vec![(
                fallback_edge_for_specs(specs, &edge_geometries),
                BlendFailure::InvalidTopology {
                    reason: format!("repair stale shell anchors: {err}"),
                    location: None,
                },
            )],
        }
    })?;
    emit!(
        progress::FilletPhase::ApplyTopology,
        t4,
        edge_geometries.len()
    );
    check_cancel!();

    let t5 = Instant::now();
    let _ = apply_result.brep.recompute_all_pcurves();
    apply_result.brep.fix_stale_pcurves();
    let _ = construct::clamp_bounded_pcurves_post_phase5(&mut apply_result.brep);
    repair_stale_shell_anchors(&mut apply_result.brep).map_err(|err| {
        BlendError::PartialFailure {
            failed: 1,
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures: vec![(
                fallback_edge_for_specs(specs, &edge_geometries),
                BlendFailure::InvalidTopology {
                    reason: format!("repair stale shell anchors after pcurve repair: {err}"),
                    location: None,
                },
            )],
        }
    })?;
    surface_boundary_evidence::audit_gregory_analytic_boundary_evidence(&mut apply_result.brep);
    construct::restore_exact_boundary_cap_pcurves_post_phase5(
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
    // Final pcurve fidelity pass: the intermediate repair stages can
    // (re)write low-fidelity pcurves after the first fix_stale run
    // (PK-095: straight 2-point pcurves on fillet-sweep trims, 5e-2
    // off, were the dominant tessellation seam crack). Deviation-gated:
    // only pcurves whose surface trace leaves the edge curve beyond
    // the accumulated imprecision are cleared and recomputed.
    apply_result.brep.fix_stale_pcurves();
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
    if let Some(ctx) = tracking {
        assign_fillet_names(&mut apply_result.brep, ctx, &pre_faces);
    }

    Ok(BlendResult {
        brep: apply_result.brep,
        edge_warnings: all_warnings,
        corner_patches: apply_result.corner_patches,
        repair_report: apply_result.repair_report,
    })
}

fn fallback_edge_for_specs(specs: &[FilletSpec], edge_geometries: &[EdgeGeometry]) -> EdgeHandle {
    specs
        .first()
        .map(|s| s.edge)
        .or_else(|| edge_geometries.first().map(|g| g.analysis.edge))
        .unwrap_or_default()
}

fn first_structural_validation_error(brep: &BRep) -> Option<String> {
    [
        brep.validate_half_edge_pairing().err(),
        brep.validate_half_edge_cycles().err(),
        brep.validate_vertex_consistency().err(),
        brep.validate_manifold_edges().err(),
        brep.validate_loop_closure().err(),
    ]
    .into_iter()
    .flatten()
    .find(|err| {
        matches!(
            err,
            popkern_core::error::ValidationError::BrokenCycle { .. }
                | popkern_core::error::ValidationError::BrokenTwin { .. }
                | popkern_core::error::ValidationError::InconsistentVertex { .. }
                | popkern_core::error::ValidationError::NonManifoldEdge { .. }
                | popkern_core::error::ValidationError::OpenLoop { .. }
        )
    })
    .map(|err| format!("post-blend structural validation failed: {err:?}"))
}

#[path = "lib_chamfer.rs"]
mod lib_chamfer;
pub use lib_chamfer::chamfer_edges;

/// Fillet the specified edges and bundle the result with post-op
/// [`BlendQuality`] metrics (`R_min`, `patch_count`, `surface_area`,
/// composite score). Same parameter validation and error taxonomy
/// as [`fillet_edges`]; quality is computed from the new-face delta
/// via [`precheck::fillet_quality`].
#[path = "lib_convenience.rs"]
mod lib_convenience;
pub use lib_convenience::{chamfer_edge, fillet_edge, fillet_edges_with_quality};

use geometry_compute::{build_overlap_inputs, compute_chamfer_geometry, compute_fillet_geometry};
use naming::assign_fillet_names;
#[cfg(test)]
#[path = "corner_warning_tests.rs"]
mod corner_warning_tests;
#[cfg(test)]
#[path = "lib_tests.rs"]
mod tests;
