//! Sphere-cap selection and seam preflight helpers for three-quad apply.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

use super::three_quad_face_matching::primary_cylinder_identity_of_face;
use crate::blend::boundary::{BlendBoundary, BoundarySegment};
use crate::fillet_surface::FilletSurfaceResult;

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
) -> Option<String> {
    if boundary.segments.len() != 2 * boundary.per_fillet_caps.len() {
        return Some(format!(
            "boundary segment count {} does not match per-fillet cap count {}",
            boundary.segments.len(),
            boundary.per_fillet_caps.len()
        ));
    }
    for (k, &(center, radius)) in boundary.per_fillet_caps.iter().enumerate() {
        for segment in [&boundary.segments[2 * k], &boundary.segments[2 * k + 1]] {
            if let Some(reason) =
                sphere_cap_segment_degeneracy(brep, fillet_surfaces, center, radius, segment)
            {
                if crate::direct_apex_trace::enabled() {
                    eprintln!(
                        "[direct-apex-trace] sphere-cap preflight declined boundary[{k}] segment={}: {reason}",
                        boundary_segment_kind(segment)
                    );
                }
                return Some(format!(
                    "boundary[{k}] {}: {reason}",
                    boundary_segment_kind(segment)
                ));
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
            if let Some(reason) =
                sphere_cap_segment_degeneracy(brep, fillet_surfaces, *center, *radius, segment)
            {
                if crate::direct_apex_trace::enabled() {
                    eprintln!(
                        "[direct-apex-trace] sphere-cap rung declined subface[{k}] segment={}: {reason}",
                        boundary_segment_kind(segment)
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
) -> Option<String> {
    match segment {
        BoundarySegment::FaceEdge { face, .. } => {
            let surface = brep.faces.get(*face)?.surface.as_ref()?;
            if let SurfaceGeometry::Plane { origin, normal, .. } = &surface.geometry {
                let normal = normal.normalized()?;
                let plane_offset = (sphere_center - *origin).dot(&normal).abs();
                let circle_radius_sq = sphere_radius * sphere_radius - plane_offset * plane_offset;
                let tol = sphere_radius.max(1.0) * 1.0e-8;
                if circle_radius_sq <= tol * tol {
                    return Some(format!(
                        "sphere-plane boundary is tangent or empty \
                         (sphere radius {sphere_radius:.6e}, plane offset {plane_offset:.6e})"
                    ));
                }
            }
        }
        BoundarySegment::FilletEnd { fillet_index, .. } => {
            let Some(fillet) = fillet_surfaces.get(*fillet_index) else {
                return Some(format!("fillet index {fillet_index} missing"));
            };
            if let SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } = &fillet.surface
            {
                let Some(axis_unit) = axis.normalized() else {
                    return Some("fillet cylinder axis is degenerate".to_string());
                };
                let to_center = sphere_center - *origin;
                let along = to_center.dot(&axis_unit);
                let projected = *origin + axis_unit * along;
                let axis_error = sphere_center.distance_to(&projected);
                let radius_error = (sphere_radius - *radius).abs();
                let tol = sphere_radius.max(*radius).max(1.0) * 1.0e-5;
                if axis_error > tol || radius_error > tol {
                    return Some(format!(
                        "sphere-cylinder tangent circle is not well-defined \
                         (axis error {axis_error:.3e}, radius error {radius_error:.3e})"
                    ));
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
