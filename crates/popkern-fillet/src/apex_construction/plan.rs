//! Phase 6c-A2 apex planner.
//!
//! `plan_apex_resolution` is the **pure** entry point that produces
//! an [`ApexPlan`] from a BRep + apex vertex + fillet specs.
//! It does NOT mutate the BRep — callers can preview, inspect
//! telemetry, or cancel.
//!
//! Algorithm (per ADR §6.2):
//! 1. Validate inputs (radii positive/finite, ≥3 fillets at apex).
//! 2. Compute boundary segments (via `compute_blend_boundary`).
//! 3. Per fillet: compute T_a[k], T_b[k] via `sphere_fit_termination`.
//! 4. Decide cap kind: Sphere if all radii equal AND bridges collapse,
//!    else Gregory.
//! 5. Build [`TargetTopology`] (apex_vertices, apex_edges, cap_face,
//!    face_boundary_updates).
//! 6. Build [`PlanTelemetry`] (cap_kind_reason, diagnostic_values).
//! 7. Return [`ApexPlan`].
//!
//! Failure paths fall back through the [`ConstructionRung`] ladder:
//! Direct → Setback → Open.  Each fallback emits a
//! [`BlendWarning::ContinuityDowngraded`] (legacy pattern).
//!
//! # Status
//!
//! Stub implementation (A2.1.a) — returns
//! `ApexPlanError::InvalidInput("planner not yet implemented")`.
//! A2.1.b through A2.1.f flesh out each step.

use popkern_core::topology::{BRep, VertexHandle};

use crate::types::{FilletRadius, FilletSpec};

use super::cap_kind::CapKind;
use super::error::ApexPlanError;
use super::target_topology::TargetTopology;
use super::telemetry::PlanTelemetry;

/// A plan for resolving an apex vertex's fillet topology.
///
/// Produced by [`plan_apex_resolution`]; consumed by
/// [`construct_apex_topology`](super::super::apex_construction::construct_apex_topology)
/// (when A2.2 lands).  Introspectable for UI preview.
///
/// See ADR `docs/decisions/2026-04-30-phase-6c-a2-direct-apex-construction.md` §3.1.
#[derive(Debug, Clone)]
pub struct ApexPlan {
    /// The apex vertex this plan resolves.
    pub apex_vertex: VertexHandle,
    /// Cap kind decision (Sphere or Gregory) with computed geometry.
    pub cap_kind: CapKind,
    /// Final topology spec — vertices, edges, faces with their geometry.
    pub target_topology: TargetTopology,
    /// Decision telemetry — explains why this plan was chosen, which
    /// rung was used, and exposes diagnostic numerics.  UI-friendly.
    pub telemetry: PlanTelemetry,
}

impl ApexPlan {
    /// Number of fillets meeting at the apex.
    pub fn n_fillets(&self) -> usize {
        self.cap_kind.n_fillets()
    }
}

/// Compute an [`ApexPlan`] for resolving an apex vertex's fillet
/// topology.  Pure — does not mutate the BRep.
///
/// # Inputs
///
/// * `brep` — the BRep being filleted (read-only).
/// * `apex` — the vertex where multiple fillets meet.
/// * `fillet_specs` — the fillet specs (one per incident fillet edge),
///   in the same order as the cap-loop walk would visit them.
/// * `tolerance` — geometric tolerance for cap-kind dispatch and
///   coincidence checks.
///
/// # Output
///
/// `Ok(ApexPlan)` on success; `Err(ApexPlanError)` describing the
/// failure mode for user-actionable feedback.
///
/// # Errors
///
/// See [`ApexPlanError`] for the full failure taxonomy.
///
/// # Status
///
/// **Stub**: A2.1.a scaffolds the type API.  Actual implementation
/// lands in A2.1.b through A2.1.f.
pub fn plan_apex_resolution(
    brep: &BRep,
    apex: VertexHandle,
    fillet_specs: &[FilletSpec],
    _tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<ApexPlan, ApexPlanError> {
    // Step 1: validate inputs (A2.1.b).
    validate_inputs(brep, apex, fillet_specs)?;
    // A2.1.c through A2.1.f flesh out the rest.
    Err(ApexPlanError::InvalidInput {
        field: "planner",
        value: "stub — A2.1.c through A2.1.f not yet implemented".into(),
    })
}

/// Validate planner inputs.  Returns `Ok(())` when every input is
/// well-formed; otherwise the most specific [`ApexPlanError`].
///
/// Checks:
/// 1. ≥ 3 fillet specs (sphere-fit cap centre needs ≥ 3 face planes).
/// 2. Apex vertex exists in the BRep.
/// 3. Each fillet's edge exists in the BRep.
/// 4. Each fillet's edge is incident to the apex vertex.
/// 5. Each fillet's radius is positive, finite, and well-formed.
fn validate_inputs(
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
