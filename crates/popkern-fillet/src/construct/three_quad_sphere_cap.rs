//! Sphere-cap selection and seam preflight helpers for three-quad apply.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

use super::three_quad_face_matching::primary_cylinder_identity_of_face;
use crate::blend::boundary::{BlendBoundary, BoundarySegment};
use crate::fillet_surface::FilletSurfaceResult;

#[derive(Debug, Clone, PartialEq)]
pub(super) struct SphereCapBoundaryDegeneracy {
    pub boundary_index: usize,
    pub segment_kind: &'static str,
    pub support_face: Option<FaceHandle>,
    pub code: &'static str,
    pub reason: String,
    pub residual: Option<f64>,
}

impl SphereCapBoundaryDegeneracy {
    pub(super) fn describe(&self) -> String {
        format!(
            "boundary[{}] {}: {}",
            self.boundary_index, self.segment_kind, self.reason
        )
    }

    pub(super) fn to_diagnostic(&self) -> crate::types::ApexConstructionDiagnostic {
        crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
            crate::types::ExactBoundaryCapDiagnostic {
                stage: "native_gregory_sphere_cap_preflight".to_string(),
                code: self.code.to_string(),
                detail_code: None,
                sector_name: format!(
                    "sphere_cap_boundary[{}].{}",
                    self.boundary_index, self.segment_kind
                ),
                edge: None,
                support_half_edge: None,
                support_face: self.support_face,
                residual: self.residual,
                outer_side: None,
                sample_parameter: None,
                tolerance: None,
                shared_fact_ids: Vec::new(),
            },
        )
    }
}

/// Selected sphere-cap rung state for the surface-assignment pass.
pub(super) struct SphereCapRungSelection {
    /// Per-sub-face sphere data (always computed when sphere-fit
    /// substitution succeeded, even if the rung is not selected).
    pub cap_sphere_per_subface: Vec<Option<(Point3, f64)>>,
    /// Whether the per-fillet sphere-cap rung replaces the Gregory
    /// sub-patches.
    pub rung: bool,
}

/// Decide whether the per-fillet sphere-cap rung applies: sphere-fit
/// substitution must have succeeded, the rung must be motivated
/// (self-intersecting Gregory intermediate, severe radius spread, or
/// the explicit env override), and every planned seam must be
/// non-degenerate on the sphere and its neighbour.
pub(super) fn select_sphere_cap_rung(
    brep: &BRep,
    boundary: &BlendBoundary,
    sub_faces: &[FaceHandle],
    fillet_surfaces: &[FilletSurfaceResult],
    intermediate_self_intersects: bool,
    sphere_cap_boundary_preflight: bool,
) -> SphereCapRungSelection {
    let sphere_fit_available = boundary.sphere_fit_substituted
        && !boundary.per_fillet_caps.is_empty()
        && boundary.per_fillet_caps.len() == sub_faces.len();
    let cap_sphere_per_subface: Vec<Option<(Point3, f64)>> = if sphere_fit_available {
        let mut matched = cap_spheres_for_subfaces(brep, sub_faces, &boundary.per_fillet_caps);
        for (k, entry) in matched.iter_mut().enumerate() {
            if entry.is_none() {
                *entry = Some(boundary.per_fillet_caps[k]);
            }
        }
        matched
    } else {
        vec![None; sub_faces.len()]
    };
    let requested = sphere_fit_available
        && (intermediate_self_intersects
            || severe_radius_spread(&boundary.per_fillet_caps)
            || std::env::var("POPKERN_SPHERE_CAP_RUNG").is_ok_and(|v| v == "1"));
    let rung = requested
        && sphere_cap_boundary_preflight
        && sphere_cap_boundary_seams_are_nondegenerate(
            brep,
            boundary,
            fillet_surfaces,
            &cap_sphere_per_subface,
        );
    SphereCapRungSelection {
        cap_sphere_per_subface,
        rung,
    }
}

/// Diagnostic for a sphere-cap seam authoring failure that the rung
/// preflight did not predict. The decomposition recommits with the
/// Gregory sub-patches; this record explains why the analytic rung
/// did not ship.
pub(super) fn sphere_cap_seam_authoring_failure_diagnostic(
    reason: &str,
) -> crate::types::ApexConstructionDiagnostic {
    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
        crate::types::ExactBoundaryCapDiagnostic {
            stage: "native_gregory_sphere_cap_rung".to_string(),
            code: "sphere_cap_seam_authoring_failed_gregory_fallback".to_string(),
            detail_code: Some(reason.to_string()),
            sector_name: "sphere_cap_rung".to_string(),
            edge: None,
            support_half_edge: None,
            support_face: None,
            residual: None,
            outer_side: None,
            sample_parameter: None,
            tolerance: None,
            shared_fact_ids: Vec::new(),
        },
    )
}

pub(super) fn severe_radius_spread(per_fillet_caps: &[(Point3, f64)]) -> bool {
    let mut min_radius = f64::INFINITY;
    let mut max_radius = 0.0_f64;
    for &(_, radius) in per_fillet_caps {
        if radius <= 0.0 {
            return false;
        }
        min_radius = min_radius.min(radius);
        max_radius = max_radius.max(radius);
    }
    min_radius.is_finite() && max_radius / min_radius >= 4.0
}

pub(super) fn cap_spheres_for_subfaces(
    brep: &BRep,
    sub_faces: &[FaceHandle],
    per_fillet_caps: &[(Point3, f64)],
) -> Vec<Option<(Point3, f64)>> {
    let mut used = vec![false; per_fillet_caps.len()];
    let mut assigned = vec![None; sub_faces.len()];

    for (face_index, &face) in sub_faces.iter().enumerate() {
        let Some((cyl_origin, cyl_axis, cyl_radius)) =
            primary_cylinder_identity_of_face(brep, face)
        else {
            continue;
        };
        let Some(axis_unit) = cyl_axis.normalized() else {
            continue;
        };
        let mut best: Option<(usize, f64)> = None;
        for (cap_index, &(centre, radius)) in per_fillet_caps.iter().enumerate() {
            if used[cap_index] {
                continue;
            }
            let radius_error = (radius - cyl_radius).abs();
            if radius_error > cyl_radius.max(radius) * 0.05 + 1.0e-9 {
                continue;
            }
            let v = centre - cyl_origin;
            let along = v.dot(&axis_unit);
            let projected = cyl_origin + axis_unit * along;
            let axis_error = centre.distance_to(&projected);
            let score = axis_error + radius_error;
            if best
                .map(|(_, best_score)| score < best_score)
                .unwrap_or(true)
            {
                best = Some((cap_index, score));
            }
        }
        if let Some((cap_index, _)) = best {
            used[cap_index] = true;
            assigned[face_index] = Some(per_fillet_caps[cap_index]);
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace] sphere-cap assignment face={face:?} -> cap[{cap_index}] r={:.3}",
                    per_fillet_caps[cap_index].1
                );
            }
        } else if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace] sphere-cap assignment face={face:?} unmatched cyl_radius={cyl_radius:.3}"
            );
        }
    }

    assigned
}

pub(super) fn sphere_cap_boundary_degeneracy_in_boundary_order(
    brep: &BRep,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Option<SphereCapBoundaryDegeneracy> {
    if boundary.segments.len() != 2 * boundary.per_fillet_caps.len() {
        return Some(SphereCapBoundaryDegeneracy {
            boundary_index: 0,
            segment_kind: "Boundary",
            support_face: None,
            code: "sphere_cap_boundary_count_mismatch",
            reason: format!(
                "boundary segment count {} does not match per-fillet cap count {}",
                boundary.segments.len(),
                boundary.per_fillet_caps.len()
            ),
            residual: None,
        });
    }
    for (k, &(center, radius)) in boundary.per_fillet_caps.iter().enumerate() {
        for segment in [&boundary.segments[2 * k], &boundary.segments[2 * k + 1]] {
            if let Some(mut degeneracy) =
                sphere_cap_segment_degeneracy(brep, fillet_surfaces, center, radius, segment)
            {
                degeneracy.boundary_index = k;
                if crate::direct_apex_trace::enabled() {
                    eprintln!(
                        "[direct-apex-trace] sphere-cap preflight declined boundary[{k}] segment={}: {reason}",
                        boundary_segment_kind(segment),
                        reason = degeneracy.reason,
                    );
                }
                return Some(degeneracy);
            }
        }
    }
    None
}

pub(super) fn sphere_cap_boundary_seams_are_nondegenerate(
    brep: &BRep,
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    cap_sphere_per_subface: &[Option<(Point3, f64)>],
) -> bool {
    if boundary.segments.len() != 2 * cap_sphere_per_subface.len() {
        return false;
    }
    for (k, cap) in cap_sphere_per_subface.iter().enumerate() {
        let Some((center, radius)) = cap else {
            return false;
        };
        for segment in [&boundary.segments[2 * k], &boundary.segments[2 * k + 1]] {
            if let Some(degeneracy) =
                sphere_cap_segment_degeneracy(brep, fillet_surfaces, *center, *radius, segment)
            {
                if crate::direct_apex_trace::enabled() {
                    eprintln!(
                        "[direct-apex-trace] sphere-cap rung declined subface[{k}] segment={}: {reason}",
                        boundary_segment_kind(segment),
                        reason = degeneracy.reason,
                    );
                }
                return false;
            }
        }
    }
    true
}

fn sphere_cap_segment_degeneracy(
    brep: &BRep,
    fillet_surfaces: &[FilletSurfaceResult],
    sphere_center: Point3,
    sphere_radius: f64,
    segment: &BoundarySegment,
) -> Option<SphereCapBoundaryDegeneracy> {
    match segment {
        BoundarySegment::FaceEdge { face, .. } => {
            // A planned support face that is missing or surfaceless is a
            // degeneracy, not a pass: the seam authoring must evaluate the
            // neighbour surface, so "cannot check" here means "cannot
            // author". (Mid-partition constructive-builder states carry
            // surfaceless support faces by design — the rung declines and
            // the Gregory sub-patches commit instead.)
            let Some(surface) = brep.faces.get(*face).and_then(|f| f.surface.as_ref()) else {
                return Some(SphereCapBoundaryDegeneracy {
                    boundary_index: 0,
                    segment_kind: "FaceEdge",
                    support_face: Some(*face),
                    code: "sphere_cap_support_surface_missing",
                    reason: format!(
                        "planned support face {face:?} is missing or has no surface; \
                         the sphere-plane seam cannot be authored against it"
                    ),
                    residual: None,
                });
            };
            if let SurfaceGeometry::Plane { origin, normal, .. } = &surface.geometry {
                let normal = normal.normalized()?;
                let plane_offset = (sphere_center - *origin).dot(&normal).abs();
                let circle_radius_sq = sphere_radius * sphere_radius - plane_offset * plane_offset;
                let tol = sphere_radius.max(1.0) * 1.0e-8;
                if circle_radius_sq <= tol * tol {
                    return Some(SphereCapBoundaryDegeneracy {
                        boundary_index: 0,
                        segment_kind: "FaceEdge",
                        support_face: Some(*face),
                        code: "sphere_plane_boundary_tangent_or_empty",
                        reason: format!(
                            "sphere-plane boundary is tangent or empty \
                             (sphere radius {sphere_radius:.6e}, plane offset {plane_offset:.6e})"
                        ),
                        residual: Some((sphere_radius - plane_offset).abs()),
                    });
                }
            }
        }
        BoundarySegment::FilletEnd { fillet_index, .. } => {
            let Some(fillet) = fillet_surfaces.get(*fillet_index) else {
                return Some(SphereCapBoundaryDegeneracy {
                    boundary_index: 0,
                    segment_kind: "FilletEnd",
                    support_face: None,
                    code: "sphere_cap_fillet_index_missing",
                    reason: format!("fillet index {fillet_index} missing"),
                    residual: None,
                });
            };
            if let SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } = &fillet.surface
            {
                let Some(axis_unit) = axis.normalized() else {
                    return Some(SphereCapBoundaryDegeneracy {
                        boundary_index: 0,
                        segment_kind: "FilletEnd",
                        support_face: None,
                        code: "sphere_cap_fillet_cylinder_axis_degenerate",
                        reason: "fillet cylinder axis is degenerate".to_string(),
                        residual: None,
                    });
                };
                let to_center = sphere_center - *origin;
                let along = to_center.dot(&axis_unit);
                let projected = *origin + axis_unit * along;
                let axis_error = sphere_center.distance_to(&projected);
                let radius_error = (sphere_radius - *radius).abs();
                let tol = sphere_radius.max(*radius).max(1.0) * 1.0e-5;
                if axis_error > tol || radius_error > tol {
                    return Some(SphereCapBoundaryDegeneracy {
                        boundary_index: 0,
                        segment_kind: "FilletEnd",
                        support_face: None,
                        code: "sphere_cylinder_boundary_not_tangent",
                        reason: format!(
                            "sphere-cylinder tangent circle is not well-defined \
                             (axis error {axis_error:.3e}, radius error {radius_error:.3e})"
                        ),
                        residual: Some(axis_error.max(radius_error)),
                    });
                }
            }
        }
    }
    None
}

fn boundary_segment_kind(segment: &BoundarySegment) -> &'static str {
    match segment {
        BoundarySegment::FilletEnd { .. } => "FilletEnd",
        BoundarySegment::FaceEdge { .. } => "FaceEdge",
    }
}
