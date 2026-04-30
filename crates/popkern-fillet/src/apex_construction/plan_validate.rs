//! Phase 6c-A2.1.b — input validation for the apex planner.
//!
//! Checks the planner's inputs are well-formed before any sphere-fit
//! math runs.  Each variant of [`ApexPlanError`] surfaced here is
//! user-actionable.

use popkern_core::topology::{BRep, VertexHandle};

use crate::types::{FilletRadius, FilletSpec};

use super::error::ApexPlanError;

/// Validate planner inputs.  Returns `Ok(())` when every input is
/// well-formed; otherwise the most specific [`ApexPlanError`].
///
/// Checks:
/// 1. ≥ 3 fillet specs (sphere-fit cap centre needs ≥ 3 face planes).
/// 2. Apex vertex exists in the BRep.
/// 3. Each fillet's edge exists in the BRep.
/// 4. Each fillet's edge is incident to the apex vertex.
/// 5. Each fillet's radius is positive, finite, and well-formed.
pub(crate) fn validate_inputs(
    brep: &BRep,
    apex: VertexHandle,
    fillet_specs: &[FilletSpec],
) -> Result<(), ApexPlanError> {
    // 1. Enough fillets.
    if fillet_specs.len() < 3 {
        return Err(ApexPlanError::InsufficientFillets {
            apex,
            n: fillet_specs.len(),
        });
    }
    // 2. Apex vertex exists.
    if brep.vertices.get(apex).is_none() {
        return Err(ApexPlanError::InvalidInput {
            field: "apex",
            value: format!("vertex {apex:?} not present in BRep"),
        });
    }
    // 3, 4, 5. Per-fillet checks.
    for (idx, spec) in fillet_specs.iter().enumerate() {
        // 3. Edge exists.
        let Some(edge) = brep.edges.get(spec.edge) else {
            return Err(ApexPlanError::InvalidInput {
                field: "fillet_specs[].edge",
                value: format!("fillet {idx}: edge {:?} not present in BRep", spec.edge),
            });
        };
        // 4. Edge incident to apex.
        let primary = edge.half_edge;
        let Some(primary_he) = brep.half_edges.get(primary) else {
            return Err(ApexPlanError::InvalidInput {
                field: "fillet_specs[].edge",
                value: format!(
                    "fillet {idx}: edge {:?}'s primary half-edge missing",
                    spec.edge
                ),
            });
        };
        let twin = primary_he.twin;
        let Some(twin_he) = brep.half_edges.get(twin) else {
            return Err(ApexPlanError::InvalidInput {
                field: "fillet_specs[].edge",
                value: format!(
                    "fillet {idx}: edge {:?}'s twin half-edge missing",
                    spec.edge
                ),
            });
        };
        if primary_he.origin != apex && twin_he.origin != apex {
            return Err(ApexPlanError::InvalidInput {
                field: "fillet_specs[].edge",
                value: format!(
                    "fillet {idx}: edge {:?} not incident to apex {apex:?}",
                    spec.edge
                ),
            });
        }
        // 5. Radius valid.
        validate_radius(idx, &spec.radius)?;
    }
    Ok(())
}

/// Validate a single fillet's radius profile.  Returns `Ok(())`
/// when all radii are positive, finite, non-NaN, non-Inf.
fn validate_radius(idx: usize, radius: &FilletRadius) -> Result<(), ApexPlanError> {
    let check = |r: f64, label: &str| -> Result<(), ApexPlanError> {
        if !r.is_finite() {
            return Err(ApexPlanError::InvalidInput {
                field: "fillet_specs[].radius",
                value: format!("fillet {idx}: {label} radius {r} is not finite (NaN or Inf)"),
            });
        }
        if r <= 0.0 {
            return Err(ApexPlanError::InvalidInput {
                field: "fillet_specs[].radius",
                value: format!("fillet {idx}: {label} radius {r} is not strictly positive"),
            });
        }
        Ok(())
    };
    match radius {
        FilletRadius::Constant(r) => check(*r, "constant")?,
        FilletRadius::Variable { start, end } => {
            check(*start, "start")?;
            check(*end, "end")?;
        }
        FilletRadius::ByPoints(points) => {
            if points.is_empty() {
                return Err(ApexPlanError::InvalidInput {
                    field: "fillet_specs[].radius",
                    value: format!("fillet {idx}: ByPoints profile has no points"),
                });
            }
            for (t, r) in points {
                if !t.is_finite() || !(0.0..=1.0).contains(t) {
                    return Err(ApexPlanError::InvalidInput {
                        field: "fillet_specs[].radius",
                        value: format!("fillet {idx}: ByPoints t={t} not in [0, 1] or not finite"),
                    });
                }
                check(*r, "ByPoints")?;
            }
        }
    }
    Ok(())
}
