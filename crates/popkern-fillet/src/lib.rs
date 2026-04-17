#![warn(missing_docs)]
//! popkern-fillet: Fillets and chamfers for the popkern BRep kernel.
//! Public entry points: [`fillet_edges`], [`chamfer_edges`], [`vertex_blend_3d`].
pub mod analysis;
pub mod blend;
pub mod blend_diagnostics;
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
pub mod curvature_comb;
pub(crate) mod cusp_subdivide;
pub(crate) mod degeneracy;
pub(crate) mod degeneracy_surface;
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
pub mod remediation;
pub mod sensitivity;
pub mod surface_validate;
pub mod trim;
pub mod types;
pub mod unfillet;
pub mod unfillet_closed;
pub(crate) mod unfillet_closed_helpers;
pub(crate) mod unfillet_helpers;
pub(crate) mod unfillet_surgery;
pub mod validate;
pub mod vertex_blend;
pub mod vertex_classifier;
pub mod walker;
pub(crate) mod walker_adaptive;
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
pub use preview::{PreviewEdgeGeometry, PreviewGeometry, preview_chamfer, preview_fillet};
pub use progress::*;
pub use types::*;
pub use unfillet::{UnfilletPlan, compute_unfillet_plan, unfillet_face, unfillet_via_record};
pub use vertex_blend::{vertex_blend_3d, vertex_blend_3d_variable};

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

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
    let mut apply_result =
        construct::apply_blend_geometry(brep, &edge_geometries, &tolerance, &options.vertex_mode)
            .map_err(|edge_failures| BlendError::PartialFailure {
            failed: edge_failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures: edge_failures,
        })?;
    emit!(
        progress::FilletPhase::ApplyTopology,
        t4,
        edge_geometries.len()
    );
    check_cancel!();

    // PHASE 5: Ensure pcurve coverage and correctness.
    let t5 = Instant::now();
    let _ = apply_result.brep.recompute_all_pcurves();
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

    // Assign stable names when tracking is active.
    if let Some(ctx) = tracking {
        assign_fillet_names(&mut apply_result.brep, ctx, &pre_faces);
    }

    Ok(BlendResult {
        brep: apply_result.brep,
        edge_warnings: all_warnings,
        corner_patches: apply_result.corner_patches,
    })
}

/// Chamfer the specified edges of a BRep solid. See [`BlendError`] for failure modes.
pub fn chamfer_edges(
    brep: &BRep,
    specs: &[ChamferSpec],
    options: &ChamferOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    // §10.1 — reject invalid parameters at the API boundary.
    validate::validate_chamfer_specs(brep, specs)?;

    // Capture pre-operation faces for naming.
    let pre_faces: std::collections::HashSet<FaceHandle> = brep.faces.keys().collect();

    // 1. Analyze edges.
    let edges: Vec<EdgeHandle> = specs.iter().map(|s| s.edge).collect();
    let analyses = analysis::analyze_edges(brep, &edges, options.propagate_tangent)
        .map_err(BlendError::NoValidEdges)?;

    // 2. Compute geometry for each edge.
    let mut edge_geometries = Vec::new();
    let mut failures = Vec::new();

    for (i, spec) in specs.iter().enumerate() {
        if i >= analyses.len() {
            continue;
        }
        match compute_chamfer_geometry(&analyses[i], &spec.kind, brep, &tolerance) {
            Ok(geom) => edge_geometries.push(geom),
            Err(failure) => failures.push((spec.edge, failure)),
        }
    }

    if edge_geometries.is_empty() && !failures.is_empty() {
        return Err(BlendError::PartialFailure {
            failed: failures.len(),
            total: specs.len(),
            partial_result: Box::new(brep.clone()),
            failures,
        });
    }

    // 3. Apply topology. Chamfers always use mitred vertex mode.
    let mut apply_result =
        construct::apply_blend_geometry(brep, &edge_geometries, &tolerance, &VertexMode::Mitred)
            .map_err(|edge_failures| BlendError::PartialFailure {
                failed: edge_failures.len(),
                total: specs.len(),
                partial_result: Box::new(brep.clone()),
                failures: edge_failures,
            })?;

    // Assign stable names when tracking is active.
    if let Some(ctx) = tracking {
        assign_fillet_names(&mut apply_result.brep, ctx, &pre_faces);
    }

    Ok(BlendResult {
        brep: apply_result.brep,
        edge_warnings: apply_result.warnings,
        corner_patches: apply_result.corner_patches,
    })
}

/// Fillet the specified edges and bundle the result with post-op
/// [`BlendQuality`] metrics (`R_min`, `patch_count`, `surface_area`,
/// composite score). Same parameter validation and error taxonomy
/// as [`fillet_edges`]; quality is computed from the new-face delta
/// via [`precheck::fillet_quality`].
pub fn fillet_edges_with_quality(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<(BlendResult, BlendQuality), BlendError> {
    let pre_faces: std::collections::HashSet<FaceHandle> = brep.faces.keys().collect();
    let result = fillet_edges(brep, specs, options, tolerance, tracking)?;
    // Use the spec's largest requested radius as the scoring
    // denominator — for constant radius this is the request, for
    // variable specs it's the peak.
    let requested_radius = specs
        .iter()
        .map(|s| geometry_compute_feasibility::largest_requested_radius(&s.radius))
        .fold(0.0_f64, f64::max);
    let quality = precheck::fillet_quality(&precheck::ComputeQualityInput {
        post_brep: &result.brep,
        pre_faces: &pre_faces,
        requested_radius,
    });
    Ok((result, quality))
}

/// Convenience: fillet a single edge with a constant radius.
/// Equivalent to [`fillet_edges`] with one [`FilletSpec`] + default options.
/// Returns a [`BlendError`] if the edge cannot be filleted.
pub fn fillet_edge(
    brep: &BRep,
    edge: EdgeHandle,
    radius: f64,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    fillet_edges(
        brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        }],
        &FilletOptions::default(),
        tolerance,
        tracking,
    )
}

/// Convenience: chamfer a single edge with equal distances.
/// Equivalent to [`chamfer_edges`] with one [`ChamferSpec`] + default options.
/// Returns a [`BlendError`] if the edge cannot be chamfered.
pub fn chamfer_edge(
    brep: &BRep,
    edge: EdgeHandle,
    distance: f64,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    chamfer_edges(
        brep,
        &[ChamferSpec {
            edge,
            kind: ChamferKind::Equal(distance),
        }],
        &ChamferOptions::default(),
        tolerance,
        tracking,
    )
}

use geometry_compute::{build_overlap_inputs, compute_chamfer_geometry, compute_fillet_geometry};
use naming::assign_fillet_names;
#[cfg(test)]
#[path = "corner_warning_tests.rs"]
mod corner_warning_tests;
#[cfg(test)]
#[path = "lib_tests.rs"]
mod tests;
