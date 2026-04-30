//! Phase 6c-A2.1.c — per-fillet sphere-fit termination computation.
//!
//! Produces a `Vec<PerFilletTermination>` from the BlendBoundary
//! and fillet surfaces, by walking each fillet's boundary segment
//! and computing T_a/T_b via the existing `sphere_fit_termination`
//! math.
//!
//! Two passes are exposed:
//! * **Strict**: returns the first per-fillet error.  Used by the
//!   Direct rung.
//! * **Lenient**: skips per-fillet failures, returning whichever
//!   fillets succeeded.  Used by the Setback rung
//!   ([`plan_rung_fallback`](super::plan_rung_fallback)).

use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::blend::boundary::{BlendBoundary, BoundarySegment};
use crate::blend::sphere_fit_centre::{
    face_outward_normal_at, fillet_nominal_radius, select_per_fillet_centre,
};
use crate::fillet_surface::FilletSurfaceResult;

use super::error::ApexPlanError;

/// One fillet's sphere-fit termination data.  Internal to the
/// planner; aggregated into the [`ApexPlan`](super::plan::ApexPlan)
/// in A2.1.d/e.
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
pub(crate) fn compute_apex_boundary(
    brep: &BRep,
    apex: VertexHandle,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<BlendBoundary, ApexPlanError> {
    use crate::blend::boundary::compute_blend_boundary;
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
/// machinery.
///
/// Per Phase 6c-A continuation §4.2: each fillet uses its own cap
/// radius (= nominal radius) and its own cap centre, yielding
/// distinct termination points for variable-radius corners.
///
/// Returns one `PerFilletTermination` per fillet, in cap-loop walk
/// order (matching `boundary.segments[2k]` — the FilletEnd
/// segments).
///
/// **Strict**: returns the first per-fillet error.  Direct rung uses
/// this.  See [`compute_per_fillet_terminations_lenient`] for the
/// Setback-rung variant that tolerates individual failures.
pub(crate) fn compute_per_fillet_terminations(
    brep: &BRep,
    apex: VertexHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<Vec<PerFilletTermination>, ApexPlanError> {
    compute_per_fillet_terminations_inner(brep, apex, boundary, fillet_surfaces, tolerance, true)
}

/// Like [`compute_per_fillet_terminations`] but skips fillets whose
/// individual sphere-fit fails instead of returning the first error.
/// Used by the Setback rung.
pub(crate) fn compute_per_fillet_terminations_lenient(
    brep: &BRep,
    apex: VertexHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<Vec<PerFilletTermination>, ApexPlanError> {
    compute_per_fillet_terminations_inner(brep, apex, boundary, fillet_surfaces, tolerance, false)
}

/// Shared implementation used by strict and lenient passes.  The
/// `strict` flag controls per-fillet error handling: strict returns
/// on the first failure; lenient skips and continues.
///
/// Boundary-shape validation, face-normal collection, and unique-face
/// counting are common to both passes — those errors always
/// propagate.
fn compute_per_fillet_terminations_inner(
    brep: &BRep,
    apex: VertexHandle,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
    strict: bool,
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

    // Validate alternation pattern — both strict and lenient
    // require this structural correctness.
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

    let apex_pos =
        brep.vertices
            .get(apex)
            .map(|v| v.position)
            .ok_or_else(|| ApexPlanError::InvalidInput {
                field: "apex",
                value: format!("vertex {apex:?} disappeared between validation and planning"),
            })?;

    // Collect face_A/face_B normals per fillet, plus the unique
    // adjacent faces (for cap-centre solver).  These structural
    // checks fail-fast in both strict and lenient.
    let mut face_a_normals: Vec<Vec3> = Vec::with_capacity(n_fil);
    let mut face_b_normals: Vec<Vec3> = Vec::with_capacity(n_fil);
    let mut face_a_handles: Vec<FaceHandle> = Vec::with_capacity(n_fil);
    let mut face_b_handles: Vec<FaceHandle> = Vec::with_capacity(n_fil);
    let mut unique_face_normals: Vec<(FaceHandle, Vec3)> = Vec::new();
    for k in 0..n_fil {
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

    if unique_face_normals.len() < 3 {
        return Err(ApexPlanError::InsufficientFillets {
            apex,
            n: unique_face_normals.len(),
        });
    }
    let face_normals: Vec<Vec3> = unique_face_normals.iter().map(|(_, n)| *n).collect();

    // Per-fillet sphere-fit termination loop.  This is where
    // strict/lenient differ.
    let mut out: Vec<PerFilletTermination> = Vec::with_capacity(n_fil);
    for k in 0..n_fil {
        let fi = match &boundary.segments[2 * k] {
            BoundarySegment::FilletEnd { fillet_index, .. } => *fillet_index,
            _ => unreachable!(),
        };
        // Strict: bail on missing fillet_surfaces entry.  Lenient: skip.
        let fillet = match fillet_surfaces.get(fi) {
            Some(f) => f,
            None if !strict => continue,
            None => {
                return Err(ApexPlanError::InvalidInput {
                    field: "fillet_surfaces",
                    value: format!("FilletEnd[{k}].fillet_index={fi} out of range"),
                });
            }
        };
        let r_k = match fillet_nominal_radius(&fillet.surface) {
            Some(r) => r,
            None if !strict => continue,
            None => {
                let kind = match &fillet.surface {
                    popkern_core::geometry::surface::SurfaceGeometry::Plane { .. } => "Plane",
                    popkern_core::geometry::surface::SurfaceGeometry::Cone { .. } => "Cone",
                    popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. } => "Sphere",
                    popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_) => "Nurbs",
                    popkern_core::geometry::surface::SurfaceGeometry::Gregory(_) => "Gregory",
                    _ => "Unknown",
                };
                return Err(ApexPlanError::UnsupportedFilletKind {
                    fillet_idx: k,
                    kind,
                });
            }
        };
        let cap_centre =
            match select_per_fillet_centre(apex_pos, &face_normals, r_k, &fillet.surface) {
                Some(c) => c,
                None if !strict => continue,
                None => {
                    return Err(ApexPlanError::SingularCapCentre {
                        apex,
                        condition_number: f64::INFINITY,
                    });
                }
            };
        let result = match sphere_fit_termination(
            apex_pos,
            cap_centre,
            r_k,
            &fillet.surface,
            face_a_normals[k],
            face_b_normals[k],
            tolerance,
        ) {
            Ok(r) => r,
            Err(_) if !strict => continue,
            Err(_) => {
                return Err(ApexPlanError::SphereCylinderNoIntersection { fillet_idx: k });
            }
        };
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
