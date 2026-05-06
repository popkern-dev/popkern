use super::*;

pub(super) fn materialize_missing_generator_anchors(
    brep: &mut BRep,
    cap_face: FaceHandle,
    fillet_ends: &[FilletEndOwnership],
    face_edges: &[FaceEdgeOwnership],
) -> Result<(), ApexBuildError> {
    for fillet_end in fillet_ends {
        let allowed_support_faces = target_cylinder_support_faces(
            brep,
            fillet_end.target_cylinder_face,
            cap_face,
            face_edges,
        );
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
        let start_point = brep
            .vertices
            .get(start)
            .map(|vertex| vertex.position)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] start endpoint {start:?} \
                         is stale while materializing generator anchors",
                    fillet_end.fillet_index
                ),
            })?;
        let end_point = brep
            .vertices
            .get(end)
            .map(|vertex| vertex.position)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition FilletEnd[{}] end endpoint {end:?} is \
                         stale while materializing generator anchors",
                    fillet_end.fillet_index
                ),
            })?;
        let start_has_anchor =
            generator_anchor_on_cylinder_loop(brep, fillet_end.target_cylinder_face, start_point)
                .is_some();
        let end_has_anchor =
            generator_anchor_on_cylinder_loop(brep, fillet_end.target_cylinder_face, end_point)
                .is_some();
        if start_has_anchor && end_has_anchor {
            continue;
        }
        let (start_anchor, end_anchor) = support_topology_generator_anchor_pair(
            brep,
            fillet_end.target_cylinder_face,
            start_point,
            end_point,
            &allowed_support_faces,
        )
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{}] cannot materialize support-aware \
                 generator anchors on target cylinder face {:?}: start {}; end {}",
                fillet_end.fillet_index,
                fillet_end.target_cylinder_face,
                generator_anchor_blocker(brep, fillet_end.target_cylinder_face, start_point),
                generator_anchor_blocker(brep, fillet_end.target_cylinder_face, end_point)
            ),
        })?;
        for (endpoint, has_anchor, anchor_point) in [
            (start, start_has_anchor, start_anchor),
            (end, end_has_anchor, end_anchor),
        ] {
            if has_anchor {
                continue;
            }
            let anchor = materialize_support_boundary_generator_anchor(
                brep,
                fillet_end.target_cylinder_face,
                endpoint,
                &anchor_point,
                &allowed_support_faces,
                fillet_end.fillet_index,
            )?;
            if std::env::var_os("POPKERN_A2_TRACE").is_some() {
                let support_boundary_hits = support_plane_faces_containing_point(
                    brep,
                    fillet_end.target_cylinder_face,
                    anchor_point.point,
                )
                .into_iter()
                .map(|support_face| {
                    let boundary_edge =
                        containing_outer_loop_edge_by_curve(brep, support_face, anchor_point.point);
                    format!("{support_face:?}: boundary_edge={boundary_edge:?}")
                })
                .collect::<Vec<_>>()
                .join(", ");
                eprintln!(
                    "[a2-trace §14.9/cylinder-anchor] materialized support-aware FilletEnd[{}] \
                     anchor {anchor:?} at ({:.6},{:.6},{:.6}) for endpoint {endpoint:?} on {:?}; \
                     support hits [{}]",
                    fillet_end.fillet_index,
                    anchor_point.point.x,
                    anchor_point.point.y,
                    anchor_point.point.z,
                    fillet_end.target_cylinder_face,
                    support_boundary_hits
                );
            }
        }
    }
    Ok(())
}
