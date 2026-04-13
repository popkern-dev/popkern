#![warn(missing_docs)]
#![allow(clippy::result_large_err)]

//! popkern-bool: Boolean set operations on BRep solids.
//!
//! Provides union, intersection, and subtraction operations on
//! [`BRep`](popkern_core::topology::BRep) solid models.
//!
//! # Pipeline
//!
//! Each boolean operation proceeds through five stages:
//!
//! 1. **Intersect** — compute face-face intersection curves.
//! 2. **Split** — trim faces along those curves.
//! 3. **Classify** — label each fragment as inside, outside, or on-boundary.
//! 4. **Select** — choose which fragments belong to the result.
//! 5. **Merge** — stitch fragments into a valid closed solid.
//!
//! # Example (placeholder — not yet implemented)
//!
//! ```rust,ignore
//! use popkern_bool::{BoolOperation, boolean};
//! use popkern_core::tolerance::Tolerance;
//!
//! let result = boolean(BoolOperation::Union, solid_a, solid_b, Tolerance::default())?;
//! ```

pub mod adaptive_precision;
pub mod cancellation;
pub mod cellular;
pub(crate) mod classify;
pub mod classify_on_deferred;
pub mod classify_on_face;
pub mod classify_propagation;
pub mod classify_winding;
pub(crate) mod clip;
pub mod coincidence;
pub mod csg;
pub mod error;
pub mod history;
pub mod imprint;
pub mod instancing;
pub(crate) mod intersect;
pub(crate) mod intersect_hierarchy;
pub(crate) mod merge;
pub mod multi_body;
pub(crate) mod operation;
pub mod ops;
pub mod parallel;
pub mod perturbation;
pub mod postprocess;
pub mod progress;
pub mod provenance;
pub(crate) mod select;
pub mod sheet;
pub(crate) mod split;
pub(crate) mod split_inner_loop;
pub mod tangent_analysis;
pub mod tolerance;
pub mod validate;

pub use classify::point_in_solid;
pub use csg::{CsgNode, csg_eval};

pub use error::{
    BoolDiagnostics, BoolError, BoolOperation, Classification, PointClassification, Stage,
};
pub use ops::interference::{InterferenceMode, InterferenceReport, interference_check};
pub use postprocess::quality::{QualityReport, QualityThresholds, quality_report};
pub use postprocess::sliver_removal::{SliverRemovalReport, remove_slivers};
pub use progress::{CaptureReporter, NullReporter, ProgressEvent, ProgressReporter};
pub use sheet::{
    SheetBooleanResult, SheetOperationKind, SideOfSheet, boolean_with_sheets,
    classify_sheet_operation, is_closed_solid, is_sheet_body, side_of_sheet,
};

#[cfg(test)]
mod error_and_tracking_tests;
#[cfg(test)]
mod integration_tests;
#[cfg(test)]
mod integration_tests_curved;
#[cfg(test)]
mod integration_tests_curved_regression;
#[cfg(test)]
mod integration_tests_hole_debug;
#[cfg(test)]
mod integration_tests_multicurve;
#[cfg(test)]
mod translate_helper_tests;

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::EntityName;
use popkern_core::{tolerance::Tolerance, topology::BRep};

/// Perform a boolean operation on two solids, consuming both inputs.
///
/// Returns the result solid on success, or a [`BoolError`] describing the failure.
///
/// # Arguments
///
/// * `op` — the boolean set operation to perform.
/// * `a` — the first (left-hand) operand solid.
/// * `b` — the second (right-hand) operand solid.
/// * `tol` — geometric tolerance for all distance and coincidence tests.
pub fn boolean(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    // §2.3: reject sheet body inputs with a structured error rather
    // than letting volume-based classification produce undefined
    // results. Future work will add boundary-relative semantics; for
    // now this is the only safe response.
    let a_is_sheet = sheet::is_sheet_body(&a);
    let b_is_sheet = sheet::is_sheet_body(&b);
    if a_is_sheet || b_is_sheet {
        return Err(BoolError::OpenShellNotSupported {
            a: a_is_sheet,
            b: b_is_sheet,
        });
    }

    // Collect input face data before the boolean consumes the inputs.
    let input_a_faces = if tracking.is_some() {
        collect_face_data(&a)
    } else {
        Vec::new()
    };
    let input_b_faces = if tracking.is_some() {
        collect_face_data(&b)
    } else {
        Vec::new()
    };

    let mut result = operation::execute_boolean(op, a.clone(), b.clone(), &tol).or_else(|e| {
        // Perturbation fallback: retry with slightly perturbed geometry
        // if the error is the kind that perturbation can fix.
        if matches!(
            e,
            BoolError::CoincidentFaces { .. }
                | BoolError::IntersectionFailed(_)
                | BoolError::ClassificationFailed { .. }
        ) {
            perturbation::boolean_with_perturbation(
                op,
                &a,
                &b,
                tol,
                &perturbation::PerturbationConfig::default(),
            )
        } else {
            Err(e)
        }
    })?;
    if let Some(ctx) = tracking {
        assign_boolean_names(&mut result, &input_a_faces, &input_b_faces, ctx, &tol);
    }
    Ok(result)
}

/// Perform a boolean operation on two solids, borrowing both inputs.
///
/// Clones `a` and `b` internally before operating. Prefer [`boolean`] when
/// ownership can be transferred to avoid the clone cost.
///
/// # Arguments
///
/// * `op` — the boolean set operation to perform.
/// * `a` — the first (left-hand) operand solid.
/// * `b` — the second (right-hand) operand solid.
/// * `tol` — geometric tolerance for all distance and coincidence tests.
/// * `tracking` — optional tracking context for stable naming.
pub fn boolean_ref(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean(op, a.clone(), b.clone(), tol, tracking)
}

/// Compute the union of two solids (A ∪ B), consuming both inputs.
///
/// Equivalent to `boolean(BoolOperation::Union, a, b, tol, tracking)`.
pub fn union(
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean(BoolOperation::Union, a, b, tol, tracking)
}

/// Compute the union of two solids (A ∪ B), borrowing both inputs.
///
/// Equivalent to `boolean_ref(BoolOperation::Union, a, b, tol, tracking)`.
pub fn union_ref(
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean_ref(BoolOperation::Union, a, b, tol, tracking)
}

/// Compute the intersection of two solids (A ∩ B), consuming both inputs.
///
/// Equivalent to `boolean(BoolOperation::Intersection, a, b, tol, tracking)`.
pub fn intersection(
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean(BoolOperation::Intersection, a, b, tol, tracking)
}

/// Compute the intersection of two solids (A ∩ B), borrowing both inputs.
///
/// Equivalent to `boolean_ref(BoolOperation::Intersection, a, b, tol, tracking)`.
pub fn intersection_ref(
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean_ref(BoolOperation::Intersection, a, b, tol, tracking)
}

/// Subtract `b` from `a` (A − B), consuming both inputs.
///
/// Equivalent to `boolean(BoolOperation::Subtraction, a, b, tol, tracking)`.
pub fn subtraction(
    a: BRep,
    b: BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean(BoolOperation::Subtraction, a, b, tol, tracking)
}

/// Subtract `b` from `a` (A − B), borrowing both inputs.
///
/// Equivalent to `boolean_ref(BoolOperation::Subtraction, a, b, tol, tracking)`.
pub fn subtraction_ref(
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, BoolError> {
    boolean_ref(BoolOperation::Subtraction, a, b, tol, tracking)
}

/// Perform a boolean operation and return per-phase telemetry alongside the result.
///
/// Identical to [`boolean`] but populates a [`BoolDiagnostics`] with
/// [`PhaseMetric`](crate::error::PhaseMetric) samples for each pipeline stage
/// (intersection, split, classification, selection, merge, validation). Use
/// this entry point for performance regression detection and tooling that
/// needs to attribute time to specific phases.
///
/// # Arguments
///
/// * `op` — the boolean set operation.
/// * `a`, `b` — operand solids (cloned).
/// * `tol` — geometric tolerance.
pub fn boolean_with_diagnostics(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
) -> Result<(BRep, BoolDiagnostics), BoolError> {
    let mut metrics = Vec::new();
    let result =
        operation::execute_boolean_with_metrics(op, a.clone(), b.clone(), &tol, &mut metrics)?;
    let diagnostics = BoolDiagnostics {
        phase_metrics: metrics,
        ..BoolDiagnostics::default()
    };
    Ok((result, diagnostics))
}

/// Perform a boolean and, on failure, return the error together with
/// a context-aware [`crate::error::Remediation`] hint (§13.4).
///
/// This wraps [`boolean`] and its automatic perturbation fallback:
/// if the pipeline fails on a perturbation-eligible error
/// (`CoincidentFaces`, `IntersectionFailed`, `ClassificationFailed`),
/// the fallback is attempted, and the final hint carries
/// `perturbation_attempted: true` so callers know not to try
/// perturbation again — the hint will recommend the next-best
/// action (typically a tolerance bump).
///
/// The hint also includes a concrete suggested tolerance `value`
/// computed from the caller-supplied `tol`, so UIs can display the
/// "try `linear = 1e-5`" suggestion without computing it themselves.
///
/// On success, the hint is `None` and the caller gets the result
/// solid plus the full diagnostics bundle (phase metrics, warnings
/// from post-check, etc.).
pub fn boolean_with_remediation(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
) -> Result<(BRep, BoolDiagnostics), (BoolError, crate::error::Remediation)> {
    // Try the straight-through pipeline first. Track whether we
    // fell through to the perturbation fallback so the final hint
    // can avoid recommending it again.
    let mut perturbation_attempted = false;
    let result = operation::execute_boolean(op, a.clone(), b.clone(), &tol).or_else(|e| {
        if matches!(
            e,
            BoolError::CoincidentFaces { .. }
                | BoolError::IntersectionFailed(_)
                | BoolError::ClassificationFailed { .. }
        ) {
            perturbation_attempted = true;
            perturbation::boolean_with_perturbation(
                op,
                a,
                b,
                tol,
                &perturbation::PerturbationConfig::default(),
            )
        } else {
            Err(e)
        }
    });

    match result {
        Ok(brep) => Ok((brep, BoolDiagnostics::default())),
        Err(err) => {
            let ctx = crate::error::RemediationContext {
                tolerance: Some(&tol),
                perturbation_attempted,
            };
            let hint = err.remediation_hint_with_context(&ctx);
            Err((err, hint))
        }
    }
}

/// Perform a boolean with cooperative cancellation support (§15).
///
/// Behaves identically to [`boolean`] but checks `cancel` between
/// pipeline phases and returns [`BoolError::Cancelled`] as soon as a
/// cancellation request is observed. The input model is unchanged.
///
/// Pair this with [`cancellation::CancellationToken`] held by a UI
/// progress dialog or signal handler.
pub fn boolean_cancellable(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    cancel: &cancellation::CancellationToken,
) -> Result<BRep, BoolError> {
    let mut metrics = Vec::new();
    operation::execute_boolean_full(
        op,
        a.clone(),
        b.clone(),
        &tol,
        &mut metrics,
        Some(cancel),
        None,
    )
}

/// Perform a boolean and report progress through a [`progress::ProgressReporter`] (§15).
///
/// Behaves identically to [`boolean`] but invokes `reporter.on_phase`
/// after every pipeline stage with a [`progress::ProgressEvent`].
/// Combine with [`boolean_cancellable`] in your own wrapper if you
/// need both progress notifications and cancellation.
pub fn boolean_with_progress(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    reporter: &mut dyn progress::ProgressReporter,
) -> Result<BRep, BoolError> {
    let mut metrics = Vec::new();
    operation::execute_boolean_full(
        op,
        a.clone(),
        b.clone(),
        &tol,
        &mut metrics,
        None,
        Some(reporter),
    )
}

/// Perform a boolean operation and enforce post-validation on the result.
///
/// Runs the standard boolean pipeline, then applies the post-check
/// validator (manifold edges, Euler-Poincaré, half-edge cycles, watertightness,
/// face orientation, and optional self-intersection). Any `Severity::Error`
/// or `Severity::Critical` finding causes the operation to return
/// [`BoolError::ValidationFailed`] with a concatenated description.
///
/// This is the strict counterpart to [`boolean`]: prefer it when the caller
/// needs a guaranteed-valid result and is willing to fail on degenerate
/// boolean cases that the lenient pipeline would tolerate.
///
/// # Arguments
///
/// * `op` — the boolean set operation.
/// * `a`, `b` — operand solids (cloned).
/// * `tol` — geometric tolerance.
/// * `check_self_intersection` — opt in to the O(n²) result self-intersection
///   check (§9.2). Off by default in the lenient path; recommended for any
///   strict workflow.
pub fn boolean_strict(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
    check_self_intersection: bool,
) -> Result<BRep, BoolError> {
    use crate::error::Severity;

    let result = boolean_ref(op, a, b, tol, None)?;
    let report = validate::post_check::post_check(&result, &tol, check_self_intersection);

    let errors: Vec<String> = report
        .findings
        .iter()
        .filter(|f| matches!(f.severity, Severity::Error | Severity::Critical))
        .map(|f| f.message.clone())
        .collect();

    if errors.is_empty() {
        Ok(result)
    } else {
        Err(BoolError::ValidationFailed(errors.join("; ")))
    }
}

/// Perform a boolean and return both the result and entity provenance.
///
/// The [`EntityMap`] tracks which input entities each output entity came from.
/// Useful for parametric modeling workflows that need to trace features
/// through boolean operations.
pub fn boolean_with_provenance(
    op: BoolOperation,
    a: &BRep,
    b: &BRep,
    tol: Tolerance,
) -> Result<(BRep, provenance::EntityMap), BoolError> {
    let result = boolean_ref(op, a, b, tol, None)?;

    // Build provenance by matching output faces to input faces via
    // centroid proximity. This is approximate but works for most cases.
    let mut map = provenance::EntityMap::default();
    for (fh, face) in &result.faces {
        if let Some(surf) = &face.surface {
            use popkern_core::geometry::surface::SurfaceEval;
            let ((u0, u1), (v0, v1)) = surf.geometry.parameter_range();
            let u = ((u0.max(-10.0)) + (u1.min(10.0))) * 0.5;
            let v = ((v0.max(-10.0)) + (v1.min(10.0))) * 0.5;
            let centroid = surf.geometry.point_at(u, v);

            // Check which input body this face is closest to
            let dist_a = a
                .faces
                .values()
                .filter_map(|f| f.surface.as_ref())
                .map(|s| {
                    let p = s.geometry.point_at(u, v);
                    p.distance_to(&centroid)
                })
                .fold(f64::INFINITY, f64::min);
            let dist_b = b
                .faces
                .values()
                .filter_map(|f| f.surface.as_ref())
                .map(|s| {
                    let p = s.geometry.point_at(u, v);
                    p.distance_to(&centroid)
                })
                .fold(f64::INFINITY, f64::min);

            let source = if dist_a <= dist_b {
                provenance::types::SourceBody::A
            } else {
                provenance::types::SourceBody::B
            };

            // Find the closest input face
            let source_face = if source == provenance::types::SourceBody::A {
                a.faces.keys().next().unwrap_or_default()
            } else {
                b.faces.keys().next().unwrap_or_default()
            };

            map.faces.push((
                fh,
                provenance::types::FaceProvenance::Kept {
                    source,
                    original_face: source_face,
                },
            ));
        }
    }

    Ok((result, map))
}

/// Collect face names, centroids, and normals from a BRep for ancestry matching.
///
/// For faces without an existing name, creates a placeholder name so that all
/// faces can participate in the matching step.
fn collect_face_data(brep: &BRep) -> Vec<(EntityName, Point3, Vec3)> {
    use popkern_core::naming::{EntityRole, OperationId, RoleCategory};
    use popkern_core::naming_criteria::capture_face;

    let mut result = Vec::new();
    for (i, (fh, _face)) in brep.faces.iter().enumerate() {
        let Some(fc) = capture_face(brep, fh) else {
            continue;
        };
        let name = if let Some(name_id) = brep.faces[fh].name {
            if let Some(n) = brep.name_pool.get(name_id) {
                n.clone()
            } else {
                // Fallback placeholder for faces with an invalid name id.
                EntityName {
                    operation: OperationId(0),
                    role: EntityRole {
                        category: RoleCategory::Generated,
                        index: i as u32,
                    },
                    ancestors: Vec::new(),
                }
            }
        } else {
            // Placeholder for unnamed faces.
            EntityName {
                operation: OperationId(0),
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: i as u32,
                },
                ancestors: Vec::new(),
            }
        };
        result.push((name, fc.centroid, fc.normal));
    }
    result
}

/// Assign stable names to faces in a boolean result by matching against inputs.
///
/// Output faces that geometrically match an input face (by centroid proximity
/// and normal alignment) receive `Modified` names with the input face as
/// ancestor. Genuinely new intersection faces receive `Generated` names.
fn assign_boolean_names(
    brep: &mut BRep,
    input_a_faces: &[(EntityName, Point3, Vec3)],
    input_b_faces: &[(EntityName, Point3, Vec3)],
    ctx: &mut TrackingContext,
    tol: &Tolerance,
) {
    use popkern_core::naming::{EntityRole, RoleCategory};
    use popkern_core::naming_criteria::{GeometricCriteria, capture_face};

    let all_inputs: Vec<&(EntityName, Point3, Vec3)> =
        input_a_faces.iter().chain(input_b_faces.iter()).collect();

    let face_handles: Vec<_> = brep.faces.keys().collect();
    let mut gen_idx: u32 = 0;
    let mut mod_idx: u32 = 0;

    for fh in &face_handles {
        let Some(fc) = capture_face(brep, *fh) else {
            continue;
        };
        let centroid = fc.centroid;
        let normal = fc.normal;

        // Find the closest input face by centroid + normal alignment.
        let generous_dist = tol.linear() * 100.0;
        let best_match = all_inputs
            .iter()
            .filter(|(_, ic, in_)| {
                centroid.distance_to(ic) < generous_dist && normal.dot(in_).abs() > 0.9
            })
            .min_by(|a, b| {
                centroid
                    .distance_to(&a.1)
                    .partial_cmp(&centroid.distance_to(&b.1))
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

        if let Some((input_name, _, _)) = best_match {
            // Modified: this face descended from an input face.
            let ancestor_id = brep.name_pool.intern(input_name.clone(), None);
            let new_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Modified,
                    index: mod_idx,
                },
                ancestors: vec![ancestor_id],
            };
            let criteria = Some(GeometricCriteria::Face(fc));
            let new_id = brep.name_pool.intern(new_name.clone(), criteria);
            brep.faces[*fh].name = Some(new_id);
            ctx.record_modified(input_name.clone(), new_name);
            mod_idx += 1;
        } else {
            // Generated: genuinely new intersection face.
            let new_name = EntityName {
                operation: ctx.operation,
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: gen_idx,
                },
                ancestors: Vec::new(),
            };
            let criteria = Some(GeometricCriteria::Face(fc));
            let new_id = brep.name_pool.intern(new_name.clone(), criteria);
            brep.faces[*fh].name = Some(new_id);
            ctx.record_generated(new_name);
            gen_idx += 1;
        }
    }
}
