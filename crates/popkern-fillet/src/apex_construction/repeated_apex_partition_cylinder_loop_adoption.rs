use super::*;

pub(super) fn adopt_fillet_end_endpoints_into_target_cylinders(
    brep: &mut BRep,
    fillet_ends: &[FilletEndOwnership],
) -> Result<(), ApexBuildError> {
    for fillet_end in fillet_ends {
        let start = brep
            .half_edges
            .get(fillet_end.cap_half_edge)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} is stale",
                    fillet_end.fillet_index, fillet_end.cap_half_edge
                ),
            })?;
        let end = half_edge_end_vertex(brep, fillet_end.cap_half_edge).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] cap half-edge {:?} has \
                     stale twin",
                    fillet_end.fillet_index, fillet_end.cap_half_edge
                ),
            }
        })?;
        adopt_vertex_into_target_cylinder_loop(
            brep,
            fillet_end.target_cylinder_face,
            start,
            fillet_end.fillet_index,
        )?;
        adopt_vertex_into_target_cylinder_loop(
            brep,
            fillet_end.target_cylinder_face,
            end,
            fillet_end.fillet_index,
        )?;
    }
    Ok(())
}

fn adopt_vertex_into_target_cylinder_loop(
    brep: &mut BRep,
    target_cylinder_face: FaceHandle,
    vertex: VertexHandle,
    fillet_index: usize,
) -> Result<(), ApexBuildError> {
    if vertex_is_on_face_outer_loop(brep, target_cylinder_face, vertex) {
        return Ok(());
    }
    let point = brep
        .vertices
        .get(vertex)
        .map(|vertex_data| vertex_data.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] endpoint vertex \
                 {vertex:?} is stale"
            ),
        })?;
    if let Some(coincident) = coincident_vertex_on_face_outer_loop(
        brep,
        target_cylinder_face,
        point,
        brep.tolerance.linear(),
    ) {
        if coincident != vertex {
            brep.weld_vertex_pair(vertex, coincident).map_err(|err| {
                ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] failed to \
                         adopt coincident cylinder-loop vertex {coincident:?} into {vertex:?}: {err:?}"
                    ),
                }
            })?;
        }
        return Ok(());
    }

    let edge = if let Some(edge) =
        containing_outer_loop_edge_by_curve(brep, target_cylinder_face, point)
    {
        edge
    } else if let Some(edge) =
        author_missing_on_cylinder_line_trim_for_endpoint(brep, target_cylinder_face, point)?
    {
        edge
    } else {
        let surface_status = if point_lies_on_cylinder_face(brep, target_cylinder_face, point) {
            "point lies on the target cylinder surface but no current trim edge owns it"
        } else {
            "point is not on the target cylinder surface"
        };
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] endpoint {vertex:?} \
                 at ({:.6},{:.6},{:.6}) does not lie on any curve-owned outer-loop edge of \
                 target cylinder face {target_cylinder_face:?}; {surface_status}",
                point.x, point.y, point.z
            ),
        });
    };
    let tolerance = brep.tolerance;
    let split = brep
        .split_edge_at_point(edge, point, &tolerance)
        .map_err(|err| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] failed to split \
                 target cylinder edge {edge:?} at endpoint {vertex:?}: {err:?}"
            ),
        })?;
    brep.weld_vertex_pair(vertex, split.vertex)
        .map_err(|err| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] failed to adopt \
                 split vertex {:?} into endpoint {vertex:?}: {err:?}",
                split.vertex
            ),
        })?;
    Ok(())
}
