//! Semantic cap identity for sphere-fit edge wiring.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

use crate::blend::boundary::BlendBoundary;

pub(super) fn semantic_cap_indices_for_subfaces(
    brep: &BRep,
    boundary: &BlendBoundary,
    sub_faces: &[FaceHandle],
) -> Vec<Option<usize>> {
    let mut assigned = vec![None; sub_faces.len()];
    if boundary.per_fillet_caps.len() != sub_faces.len() {
        return assigned;
    }

    let mut used = vec![false; boundary.per_fillet_caps.len()];
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
        for (cap_index, &(centre, radius)) in boundary.per_fillet_caps.iter().enumerate() {
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
            assigned[face_index] = Some(cap_index);
        }
    }
    assign_remaining_cap_indices(&mut assigned, &used);
    assigned
}

fn assign_remaining_cap_indices(assigned: &mut [Option<usize>], used: &[bool]) {
    let remaining_faces: Vec<usize> = assigned
        .iter()
        .enumerate()
        .filter_map(|(face_index, cap_index)| cap_index.is_none().then_some(face_index))
        .collect();
    let remaining_caps: Vec<usize> = used
        .iter()
        .enumerate()
        .filter_map(|(cap_index, is_used)| (!is_used).then_some(cap_index))
        .collect();
    if remaining_faces.len() == remaining_caps.len() {
        for (face_index, cap_index) in remaining_faces.into_iter().zip(remaining_caps) {
            assigned[face_index] = Some(cap_index);
        }
    }
}

fn primary_cylinder_identity_of_face(
    brep: &BRep,
    face: FaceHandle,
) -> Option<(Point3, popkern_core::geometry::Vec3, f64)> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;

    type CylKey = (i64, i64, i64, i64, i64, i64);
    type CylBucket = (Point3, popkern_core::geometry::Vec3, f64, usize);
    let mut counts: std::collections::BTreeMap<CylKey, CylBucket> =
        std::collections::BTreeMap::new();
    let mut cursor = start;
    let mut safety = brep.half_edges.len() + 1;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        let he = brep.half_edges.get(cursor)?;
        let next = he.next;
        let twin_loop = brep.half_edges.get(he.twin)?.loop_;
        let twin_face = brep.loops.get(twin_loop)?.face;
        if let Some(SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        }) = brep
            .faces
            .get(twin_face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry)
        {
            let key = (
                (origin.x * 1e6).round() as i64,
                (origin.y * 1e6).round() as i64,
                (origin.z * 1e6).round() as i64,
                (axis.x * 1e6).round() as i64,
                (axis.y * 1e6).round() as i64,
                (axis.z * 1e6).round() as i64,
            );
            counts.entry(key).or_insert((*origin, *axis, *radius, 0)).3 += 1;
        }
        cursor = next;
        if cursor == start {
            break;
        }
    }

    counts
        .values()
        .copied()
        .max_by_key(|(_, _, _, count)| *count)
        .map(|(origin, axis, radius, _)| (origin, axis, radius))
}
