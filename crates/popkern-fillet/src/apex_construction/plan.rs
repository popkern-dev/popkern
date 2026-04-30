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

use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::blend::boundary::{BlendBoundary, BoundarySegment, compute_blend_boundary};
use crate::blend::sphere_fit_boundary::{
    face_outward_normal_at, fillet_nominal_radius, select_per_fillet_centre,
};
use crate::fillet_surface::FilletSurfaceResult;
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
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<ApexPlan, ApexPlanError> {
    // Step 1: validate inputs (A2.1.b).
    validate_inputs(brep, apex, fillet_specs)?;
    // Step 2: compute boundary segments (A2.1.c).
    let boundary = compute_apex_boundary(brep, apex, fillet_surfaces, tolerance)?;
    // Step 3: compute per-fillet T_a, T_b, FilletEnd curve (A2.1.c).
    let _terminations =
        compute_per_fillet_terminations(brep, apex, &boundary, fillet_surfaces, tolerance)?;
    // A2.1.d through A2.1.f flesh out the rest.
    Err(ApexPlanError::InvalidInput {
        field: "planner",
        value: "stub — A2.1.d (cap kind), A2.1.e (target topology), A2.1.f (telemetry) not yet implemented".into(),
    })
}

/// One fillet's sphere-fit termination data.  Internal to the
/// planner; aggregated into the [`ApexPlan`] in A2.1.d/e.
///
/// Per ADR §3.1: each fillet has its own pair of termination
/// points (T_a on face_A side, T_b on face_B side) plus the
/// FilletEnd 3-D arc connecting them.
#[derive(Debug, Clone)]
#[allow(dead_code)] // fields consumed by A2.1.d/e (cap kind decision + topology assembly)
pub(crate) struct PerFilletTermination {
    /// Index into `fillet_surfaces` (this fillet's surface).
    pub fillet_idx: usize,
    /// T_a — termination on the face_A (CCW-previous) side.
    pub endpoint_a: Point3,
    /// T_b — termination on the face_B (CCW-next) side.
    pub endpoint_b: Point3,
    /// Cap sphere centre for this fillet (per Phase 6c-A's
    /// per-fillet sphere-fit; uniform for symmetric corners).
    pub cap_centre: Point3,
    /// Cap sphere radius (= fillet's nominal radius).
    pub cap_radius: f64,
    /// FilletEnd arc 3-D curve, going T_a → T_b along the
    /// fillet's cylinder.
    pub curve_3d: NurbsCurve,
    /// Outward normal of face_A at the apex.
    pub face_a_normal: Vec3,
    /// Outward normal of face_B at the apex.
    pub face_b_normal: Vec3,
    /// Face_A handle (CCW-previous face shared with fillet k-1).
    pub face_a: FaceHandle,
    /// Face_B handle (CCW-next face shared with fillet k+1).
    pub face_b: FaceHandle,
}

/// Compute the apex blend boundary using the existing
/// `compute_blend_boundary` infrastructure.  Translates the
/// returned `Result<BlendBoundary, String>` into our typed error.
fn compute_apex_boundary(
    brep: &BRep,
    apex: VertexHandle,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<BlendBoundary, ApexPlanError> {
    let local_indices: Vec<usize> = (0..fillet_surfaces.len()).collect();
    compute_blend_boundary(apex, &local_indices, fillet_surfaces, brep, tolerance).map_err(|e| {
        ApexPlanError::InvalidInput {
            field: "boundary",
            value: format!("compute_blend_boundary failed at apex {apex:?}: {e}"),
        }
    })
}

/// Compute per-fillet sphere-fit terminations T_a[k], T_b[k] plus
/// the FilletEnd arc curve3d.  Reuses `select_per_fillet_centre`
/// and `sphere_fit_termination` from the existing substitution
/// machinery — the **math** is correct in the legacy path; only
/// the substitute-then-patch retrofit is being replaced.
///
/// Per Phase 6c-A continuation §4.2: each fillet uses its own cap
/// radius (= nominal radius) and its own cap centre, yielding
/// distinct termination points for variable-radius corners.
///
/// Returns one `PerFilletTermination` per fillet, in cap-loop walk
/// order (matching `boundary.segments[2k]` — the FilletEnd
/// segments).
fn compute_per_fillet_terminations(
    brep: &BRep,
    apex: VertexHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<Vec<PerFilletTermination>, ApexPlanError> {
    use popkern_core::geometry::apex_blend::sphere_fit::sphere_fit_termination;

    let n_seg = boundary.segments.len();
    if n_seg < 4 || !n_seg.is_multiple_of(2) {
        return Err(ApexPlanError::InvalidInput {
            field: "boundary.segments",
            value: format!(
                "expected alternating FilletEnd/FaceEdge segments (n_seg even, ≥4); got n_seg={n_seg}"
            ),
        });
    }
    let n_fil = n_seg / 2;

    // Validate alternation pattern.  Index 2k is FilletEnd; 2k+1 is FaceEdge.
    for (i, seg) in boundary.segments.iter().enumerate() {
        match (i.is_multiple_of(2), seg) {
            (true, BoundarySegment::FilletEnd { .. }) => {}
            (false, BoundarySegment::FaceEdge { .. }) => {}
            _ => {
                return Err(ApexPlanError::InvalidInput {
                    field: "boundary.segments",
                    value: format!("segment {i} doesn't follow FilletEnd/FaceEdge alternation"),
                });
            }
        }
    }

    // Step A: collect face_A/face_B normals per fillet, plus the
    // unique adjacent faces (for cap-centre solver).
    let apex_pos =
        brep.vertices
            .get(apex)
            .map(|v| v.position)
            .ok_or_else(|| ApexPlanError::InvalidInput {
                field: "apex",
                value: format!("vertex {apex:?} disappeared between validation and planning"),
            })?;

    let mut face_a_normals: Vec<Vec3> = Vec::with_capacity(n_fil);
    let mut face_b_normals: Vec<Vec3> = Vec::with_capacity(n_fil);
    let mut face_a_handles: Vec<FaceHandle> = Vec::with_capacity(n_fil);
    let mut face_b_handles: Vec<FaceHandle> = Vec::with_capacity(n_fil);
    let mut unique_face_normals: Vec<(FaceHandle, Vec3)> = Vec::new();

    for k in 0..n_fil {
        // FaceEdge BEFORE FilletEnd[k] is at boundary segment 2k-1
        // (or n_seg-1 for k=0 — wraparound).
        let prev_idx = if k == 0 { n_seg - 1 } else { 2 * k - 1 };
        let next_idx = 2 * k + 1;
        let face_a = match &boundary.segments[prev_idx] {
            BoundarySegment::FaceEdge { face, .. } => *face,
            _ => unreachable!("alternation already validated"),
        };
        let face_b = match &boundary.segments[next_idx] {
            BoundarySegment::FaceEdge { face, .. } => *face,
            _ => unreachable!("alternation already validated"),
        };
        let n_a = face_outward_normal_at(brep, face_a, &apex_pos).ok_or_else(|| {
            ApexPlanError::InvalidInput {
                field: "boundary.face_a",
                value: format!(
                    "could not evaluate outward normal at apex for face {face_a:?} (fillet {k})"
                ),
            }
        })?;
        let n_b = face_outward_normal_at(brep, face_b, &apex_pos).ok_or_else(|| {
            ApexPlanError::InvalidInput {
                field: "boundary.face_b",
                value: format!(
                    "could not evaluate outward normal at apex for face {face_b:?} (fillet {k})"
                ),
            }
        })?;
        face_a_normals.push(n_a);
        face_b_normals.push(n_b);
        face_a_handles.push(face_a);
        face_b_handles.push(face_b);
        for (fh, n) in [(face_a, n_a), (face_b, n_b)] {
            if !unique_face_normals.iter().any(|(f, _)| *f == fh) {
                unique_face_normals.push((fh, n));
            }
        }
    }

    // Need ≥ 3 unique face planes through the apex to solve the
    // cap-centre system (degenerate-case enumeration §4.4).
    if unique_face_normals.len() < 3 {
        return Err(ApexPlanError::InsufficientFillets {
            apex,
            n: unique_face_normals.len(),
        });
    }
    let face_normals: Vec<Vec3> = unique_face_normals.iter().map(|(_, n)| *n).collect();

    // Step B: per-fillet sphere-fit termination.
    let mut out: Vec<PerFilletTermination> = Vec::with_capacity(n_fil);
    for k in 0..n_fil {
        let fi = match &boundary.segments[2 * k] {
            BoundarySegment::FilletEnd { fillet_index, .. } => *fillet_index,
            _ => unreachable!(),
        };
        let fillet = fillet_surfaces
            .get(fi)
            .ok_or_else(|| ApexPlanError::InvalidInput {
                field: "fillet_surfaces",
                value: format!("FilletEnd[{k}].fillet_index={fi} out of range"),
            })?;
        let r_k = fillet_nominal_radius(&fillet.surface).ok_or_else(|| {
            // Non-cylinder, non-torus fillet — A2 v1 doesn't support these.
            let kind = match &fillet.surface {
                popkern_core::geometry::surface::SurfaceGeometry::Plane { .. } => "Plane",
                popkern_core::geometry::surface::SurfaceGeometry::Cone { .. } => "Cone",
                popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. } => "Sphere",
                popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_) => "Nurbs",
                popkern_core::geometry::surface::SurfaceGeometry::Gregory(_) => "Gregory",
                _ => "Unknown",
            };
            ApexPlanError::UnsupportedFilletKind {
                fillet_idx: k,
                kind,
            }
        })?;
        // Cap-centre system is singular when the face planes are
        // (nearly) coplanar.  Condition number not yet computed by
        // select_per_fillet_centre; report large sentinel value.
        let cap_centre = select_per_fillet_centre(apex_pos, &face_normals, r_k, &fillet.surface)
            .ok_or(ApexPlanError::SingularCapCentre {
                apex,
                condition_number: f64::INFINITY,
            })?;
        let result = sphere_fit_termination(
            apex_pos,
            cap_centre,
            r_k,
            &fillet.surface,
            face_a_normals[k],
            face_b_normals[k],
            tolerance,
        )
        .map_err(|_| ApexPlanError::SphereCylinderNoIntersection { fillet_idx: k })?;
        out.push(PerFilletTermination {
            fillet_idx: fi,
            endpoint_a: result.endpoint_a,
            endpoint_b: result.endpoint_b,
            cap_centre,
            cap_radius: r_k,
            curve_3d: result.curve_3d,
            face_a_normal: face_a_normals[k],
            face_b_normal: face_b_normals[k],
            face_a: face_a_handles[k],
            face_b: face_b_handles[k],
        });
    }
    Ok(out)
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
