use super::*;

pub(super) struct GeneratorTrimSupportOwnerQuery<'a> {
    pub(super) cylinder_face: FaceHandle,
    pub(super) cap_face: FaceHandle,
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
    pub(super) label: &'static str,
    pub(super) fillet_index: usize,
    pub(super) allowed_support_faces: &'a std::collections::BTreeSet<FaceHandle>,
}

pub(super) fn support_face_for_generator_trim_side(
    brep: &BRep,
    query: GeneratorTrimSupportOwnerQuery<'_>,
) -> Result<FaceHandle, ApexBuildError> {
    let GeneratorTrimSupportOwnerQuery {
        cylinder_face,
        cap_face,
        start,
        end,
        label,
        fillet_index,
        allowed_support_faces,
    } = query;
    let start_point = brep.vertices[start].position;
    let end_point = brep.vertices[end].position;
    let tol = brep.tolerance.linear().max(1.0e-9);
    let candidates = brep
        .faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cylinder_face || face == cap_face {
                return None;
            }
            if !allowed_support_faces.contains(&face) {
                return None;
            }
            let SurfaceGeometry::Plane { origin, normal, .. } = face_data
                .surface
                .as_ref()
                .map(|surface| &surface.geometry)?
            else {
                return None;
            };
            let start_distance = (start_point - *origin).dot(normal).abs();
            let end_distance = (end_point - *origin).dot(normal).abs();
            if start_distance > tol || end_distance > tol {
                return None;
            }
            let start_on_boundary = vertex_is_on_face_outer_loop(brep, face, start)
                || containing_outer_loop_edge_by_curve(brep, face, start_point).is_some();
            let end_on_boundary = vertex_is_on_face_outer_loop(brep, face, end)
                || containing_outer_loop_edge_by_curve(brep, face, end_point).is_some();
            (start_on_boundary || end_on_boundary).then_some(face)
        })
        .collect::<Vec<_>>();
    if candidates.len() == 1 {
        return Ok(candidates[0]);
    }
    if candidates.is_empty() {
        let planar_candidates = generator_planar_support_owner_candidates(
            brep,
            cylinder_face,
            cap_face,
            start_point,
            end_point,
            allowed_support_faces,
        );
        if planar_candidates.len() == 1 {
            return Ok(planar_candidates[0]);
        }
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] cannot choose a unique \
             support owner for generator trim {label} {start:?}->{end:?}: boundary support \
             candidates {:?} from allowed faces {:?}. {}. Generator trim support ownership must \
             come from the real support face adjacent to the target cylinder, not from a coplanar \
             duplicate.",
            candidates,
            allowed_support_faces,
            generator_support_owner_candidate_diagnostic(
                brep,
                start_point,
                end_point,
                allowed_support_faces,
            )
        ),
    })
}

fn generator_planar_support_owner_candidates(
    brep: &BRep,
    cylinder_face: FaceHandle,
    cap_face: FaceHandle,
    start_point: Point3,
    end_point: Point3,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Vec<FaceHandle> {
    let tol = brep.tolerance.linear().max(1.0e-9);
    brep.faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cylinder_face || face == cap_face || !allowed_support_faces.contains(&face) {
                return None;
            }
            let SurfaceGeometry::Plane { origin, normal, .. } = face_data
                .surface
                .as_ref()
                .map(|surface| &surface.geometry)?
            else {
                return None;
            };
            let start_distance = (start_point - *origin).dot(normal).abs();
            let end_distance = (end_point - *origin).dot(normal).abs();
            (start_distance <= tol && end_distance <= tol).then_some(face)
        })
        .collect()
}

fn generator_support_owner_candidate_diagnostic(
    brep: &BRep,
    start_point: Point3,
    end_point: Point3,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> String {
    let tol = brep.tolerance.linear().max(1.0e-9);
    let reports = allowed_support_faces
        .iter()
        .map(|face| {
            let Some(face_data) = brep.faces.get(*face) else {
                return format!("{face:?}: stale face");
            };
            let Some(SurfaceGeometry::Plane { origin, normal, .. }) =
                face_data.surface.as_ref().map(|surface| &surface.geometry)
            else {
                return format!("{face:?}: non-plane");
            };
            let start_distance = (start_point - *origin).dot(normal).abs();
            let end_distance = (end_point - *origin).dot(normal).abs();
            let start_edge = containing_outer_loop_edge_by_curve(brep, *face, start_point);
            let end_edge = containing_outer_loop_edge_by_curve(brep, *face, end_point);
            format!(
                "{face:?}: distance={start_distance:.3e}/{end_distance:.3e} tol={tol:.3e} \
                 start_edge={start_edge:?} end_edge={end_edge:?}"
            )
        })
        .collect::<Vec<_>>()
        .join("; ");
    format!("allowed support diagnostics [{reports}]")
}

pub(super) fn unique_far_side_support_face_for_trim_side(
    brep: &BRep,
    cylinder_face: FaceHandle,
    half_edge: HalfEdgeHandle,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Result<FaceHandle, ApexBuildError> {
    let start = brep
        .half_edges
        .get(half_edge)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot find far-side support owner for stale \
                 trim half-edge {half_edge:?}"
            ),
        })?;
    let end =
        half_edge_end_vertex(brep, half_edge).ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot find far-side support owner for trim \
             half-edge {half_edge:?}: stale twin"
            ),
        })?;
    let start_point = brep.vertices[start].position;
    let end_point = brep.vertices[end].position;
    let tol = brep.tolerance.linear().max(1.0e-9);
    let mut candidates = brep
        .faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cylinder_face {
                return None;
            }
            if !allowed_support_faces.contains(&face) {
                return None;
            }
            let SurfaceGeometry::Plane { origin, normal, .. } = face_data
                .surface
                .as_ref()
                .map(|surface| &surface.geometry)?
            else {
                return None;
            };
            let start_distance = (start_point - *origin).dot(normal).abs();
            let end_distance = (end_point - *origin).dot(normal).abs();
            (start_distance <= tol && end_distance <= tol).then_some(face)
        })
        .collect::<Vec<_>>();
    candidates.sort();
    candidates.dedup();
    if candidates.len() == 1 {
        return Ok(candidates[0]);
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition cannot choose a unique far-side support face for \
             generated trim {half_edge:?} {start:?}->{end:?}: adjacent-support candidates {:?} \
             from allowed faces {:?}. Far-side support ownership must stay within the two \
             FaceEdge support owners adjacent to the FilletEnd.",
            candidates, allowed_support_faces
        ),
    })
}
