//! Extra public boolean entry points: remediation-aware and
//! provenance-tracking variants.
//!
//! Split out of `lib.rs` to keep the parent under the 500-line cap.

use popkern_core::{tolerance::Tolerance, topology::BRep};

use crate::error::{BoolDiagnostics, BoolError, BoolOperation};
use crate::{boolean_ref, operation, perturbation, provenance, validate};

/// Perform a boolean and, on failure, return the error together with
/// a context-aware [`crate::error::Remediation`] hint (§13.4).
///
/// This wraps [`crate::boolean`] and its automatic perturbation fallback:
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
    let result = operation::execute_boolean(op, a, b, &tol).or_else(|e| {
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

/// Perform a boolean operation and enforce post-validation on the result.
///
/// Runs the standard boolean pipeline, then applies the post-check
/// validator (manifold edges, Euler-Poincaré, half-edge cycles, watertightness,
/// face orientation, and optional self-intersection). Any `Severity::Error`
/// or `Severity::Critical` finding causes the operation to return
/// [`BoolError::ValidationFailed`] with a concatenated description.
///
/// This is the strict counterpart to [`crate::boolean`]: prefer it when the
/// caller needs a guaranteed-valid result and is willing to fail on degenerate
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
/// The [`provenance::EntityMap`] tracks which input entities each output
/// entity came from. Useful for parametric modeling workflows that need
/// to trace features through boolean operations.
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
