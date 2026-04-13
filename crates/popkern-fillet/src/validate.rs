//! Parameter validation at the public API boundary.
//!
//! Every entry point in [`crate`] funnels its caller-supplied
//! arguments through these validators before touching any geometry.
//! Failures return [`BlendError::InvalidParameter`] with a specific
//! parameter name, the offending value rendered as a string, and a
//! human-readable reason. The input BRep is guaranteed unchanged
//! when validation fails — validators are pure and side-effect-free.
//!
//! §10.1 "Parameter validation" requires the engine to reject
//! non-finite, NaN, and out-of-range numbers at the API boundary
//! rather than let them fail deep in the geometry pipeline with
//! cryptic errors. This module is the single place that enforces
//! that contract.

use popkern_core::topology::BRep;

use crate::error::BlendError;
use crate::types::{ChamferKind, ChamferSpec, FilletRadius, FilletSpec};

/// Validate a batch of fillet specs against a BRep.
///
/// Checks, in order:
///
/// 1. Spec list is non-empty.
/// 2. Every `edge` handle exists in `brep`.
/// 3. Every `radius` is a finite, strictly-positive value.
/// 4. `FilletRadius::Variable` start/end values are both positive.
/// 5. `FilletRadius::ByPoints` has at least one point; all points
///    have positive radii and strictly-monotonic-increasing
///    parameter values in `[0, 1]`.
///
/// Returns `Ok(())` on success or a specific
/// [`BlendError::InvalidParameter`] on the first failure. Validators
/// stop at the first problem so callers can iterate: fix, retry.
pub fn validate_fillet_specs(brep: &BRep, specs: &[FilletSpec]) -> Result<(), BlendError> {
    if specs.is_empty() {
        return Err(BlendError::InvalidParameter {
            parameter: "specs".to_string(),
            value: "[]".to_string(),
            reason: "at least one fillet spec is required".to_string(),
        });
    }
    for (i, spec) in specs.iter().enumerate() {
        validate_edge_handle(brep, spec.edge, &format!("specs[{i}].edge"))?;
        validate_fillet_radius(&spec.radius, &format!("specs[{i}].radius"))?;
    }
    Ok(())
}

/// Validate a batch of chamfer specs against a BRep.
pub fn validate_chamfer_specs(brep: &BRep, specs: &[ChamferSpec]) -> Result<(), BlendError> {
    if specs.is_empty() {
        return Err(BlendError::InvalidParameter {
            parameter: "specs".to_string(),
            value: "[]".to_string(),
            reason: "at least one chamfer spec is required".to_string(),
        });
    }
    for (i, spec) in specs.iter().enumerate() {
        validate_edge_handle(brep, spec.edge, &format!("specs[{i}].edge"))?;
        validate_chamfer_kind(brep, &spec.kind, &format!("specs[{i}].kind"))?;
    }
    Ok(())
}

/// Validate a single fillet radius spec.
pub fn validate_fillet_radius(radius: &FilletRadius, param: &str) -> Result<(), BlendError> {
    match radius {
        FilletRadius::Constant(r) => validate_positive_finite(*r, param),
        FilletRadius::Variable { start, end } => {
            validate_positive_finite(*start, &format!("{param}.start"))?;
            validate_positive_finite(*end, &format!("{param}.end"))?;
            Ok(())
        }
        FilletRadius::ByPoints(pts) => {
            if pts.is_empty() {
                return Err(BlendError::InvalidParameter {
                    parameter: param.to_string(),
                    value: "[]".to_string(),
                    reason: "ByPoints requires at least one (t, r) pair".to_string(),
                });
            }
            let mut prev_t: Option<f64> = None;
            for (i, (t, r)) in pts.iter().enumerate() {
                let pname = format!("{param}[{i}]");
                if !t.is_finite() {
                    return Err(BlendError::InvalidParameter {
                        parameter: format!("{pname}.t"),
                        value: format!("{t}"),
                        reason: "parameter must be finite".to_string(),
                    });
                }
                if !(0.0..=1.0).contains(t) {
                    return Err(BlendError::InvalidParameter {
                        parameter: format!("{pname}.t"),
                        value: format!("{t}"),
                        reason: "parameter must be in [0, 1]".to_string(),
                    });
                }
                if let Some(prev) = prev_t
                    && *t <= prev
                {
                    return Err(BlendError::InvalidParameter {
                        parameter: format!("{pname}.t"),
                        value: format!("{t}"),
                        reason: format!(
                            "parameters must be strictly increasing (previous was {prev})"
                        ),
                    });
                }
                prev_t = Some(*t);
                validate_positive_finite(*r, &format!("{pname}.r"))?;
            }
            Ok(())
        }
    }
}

/// Validate a single chamfer spec variant.
pub fn validate_chamfer_kind(
    brep: &BRep,
    kind: &ChamferKind,
    param: &str,
) -> Result<(), BlendError> {
    match kind {
        ChamferKind::Equal(d) => validate_positive_finite(*d, &format!("{param}.distance")),
        ChamferKind::TwoDistance {
            face_1,
            dist_1,
            dist_2,
        } => {
            validate_face_handle(brep, *face_1, &format!("{param}.face_1"))?;
            validate_positive_finite(*dist_1, &format!("{param}.dist_1"))?;
            validate_positive_finite(*dist_2, &format!("{param}.dist_2"))?;
            Ok(())
        }
        ChamferKind::DistanceAngle {
            face,
            distance,
            angle,
        } => {
            validate_face_handle(brep, *face, &format!("{param}.face"))?;
            validate_positive_finite(*distance, &format!("{param}.distance"))?;
            if !angle.is_finite() {
                return Err(BlendError::InvalidParameter {
                    parameter: format!("{param}.angle"),
                    value: format!("{angle}"),
                    reason: "angle must be finite".to_string(),
                });
            }
            if *angle <= 0.0 || *angle >= std::f64::consts::PI {
                return Err(BlendError::InvalidParameter {
                    parameter: format!("{param}.angle"),
                    value: format!("{angle}"),
                    reason: "angle must be in (0, π) radians — exclusive".to_string(),
                });
            }
            Ok(())
        }
    }
}

fn validate_positive_finite(v: f64, param: &str) -> Result<(), BlendError> {
    if v.is_nan() {
        return Err(BlendError::InvalidParameter {
            parameter: param.to_string(),
            value: "NaN".to_string(),
            reason: "value must be a finite number, not NaN".to_string(),
        });
    }
    if !v.is_finite() {
        return Err(BlendError::InvalidParameter {
            parameter: param.to_string(),
            value: format!("{v}"),
            reason: "value must be finite (not ±∞)".to_string(),
        });
    }
    if v <= 0.0 {
        return Err(BlendError::InvalidParameter {
            parameter: param.to_string(),
            value: format!("{v}"),
            reason: "value must be strictly positive (> 0)".to_string(),
        });
    }
    Ok(())
}

fn validate_edge_handle(
    brep: &BRep,
    edge: popkern_core::topology::EdgeHandle,
    param: &str,
) -> Result<(), BlendError> {
    if !brep.edges.contains_key(edge) {
        return Err(BlendError::InvalidParameter {
            parameter: param.to_string(),
            value: format!("{edge:?}"),
            reason: "edge handle does not exist in the BRep".to_string(),
        });
    }
    Ok(())
}

fn validate_face_handle(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
    param: &str,
) -> Result<(), BlendError> {
    if !brep.faces.contains_key(face) {
        return Err(BlendError::InvalidParameter {
            parameter: param.to_string(),
            value: format!("{face:?}"),
            reason: "face handle does not exist in the BRep".to_string(),
        });
    }
    Ok(())
}

#[cfg(test)]
#[path = "validate_tests.rs"]
mod tests;
