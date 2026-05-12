use super::*;

pub(super) fn materialized_far_side_trim_sides(
    brep: &mut BRep,
    cylinder_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Result<Vec<SupportOwnedTrimSide>, ApexBuildError> {
    if let Ok(support_face) = unique_support_face_for_trim_endpoints(
        brep,
        cylinder_face,
        start,
        end,
        label,
        fillet_index,
        allowed_support_faces,
    ) {
        let (half_edge, actual_support_face, support_side_owned, support_source) =
            materialized_generator_trim_side(
                brep,
                GeneratorTrimMaterializationRequest {
                    cylinder_face,
                    support_face,
                    start,
                    end,
                    label,
                    fillet_index,
                    side_support_bridge_loops: &[],
                },
            )?;
        return Ok(vec![SupportOwnedTrimSide {
            cylinder_face,
            support_face: actual_support_face,
            fillet_index,
            half_edge,
            support_side_owned,
            support_source,
            requires_support_ownership: true,
        }]);
    }

    let (edge, half_edge, twin) = insert_unwired_edge_pair(brep, start, end);
    match exact_cylinder_trim_curve_between(brep, cylinder_face, start, end, label, fillet_index) {
        Ok(curve) => {
            match unique_far_side_support_face_for_trim_side(
                brep,
                cylinder_face,
                half_edge,
                allowed_support_faces,
            ) {
                Ok(support_face) => {
                    brep.edges[edge].curve = Some(curve);
                    return Ok(vec![SupportOwnedTrimSide {
                        cylinder_face,
                        support_face,
                        fillet_index,
                        half_edge,
                        support_side_owned: false,
                        support_source: None,
                        requires_support_ownership: true,
                    }]);
                }
                Err(_) => {
                    brep.half_edges.remove(half_edge);
                    if brep.half_edges.contains_key(twin) {
                        brep.half_edges.remove(twin);
                    }
                    brep.edges.remove(edge);
                }
            }
        }
        Err(_) => {
            brep.half_edges.remove(half_edge);
            if brep.half_edges.contains_key(twin) {
                brep.half_edges.remove(twin);
            }
            brep.edges.remove(edge);
        }
    }

    let start_point = brep.vertices[start].position;
    let end_point = brep.vertices[end].position;
    let start_supports = support_plane_faces_containing_point(brep, cylinder_face, start_point)
        .into_iter()
        .filter(|support_face| allowed_support_faces.contains(support_face))
        .collect::<Vec<_>>();
    let end_supports = support_plane_faces_containing_point(brep, cylinder_face, end_point)
        .into_iter()
        .filter(|support_face| allowed_support_faces.contains(support_face))
        .collect::<Vec<_>>();
    let mut candidates = Vec::new();
    for start_support in &start_supports {
        for end_support in &end_supports {
            for joint_point in
                plane_plane_cylinder_joint_points(brep, cylinder_face, *start_support, *end_support)
            {
                let score =
                    joint_point.distance_to(&start_point) + joint_point.distance_to(&end_point);
                candidates.push((*start_support, *end_support, joint_point, score));
            }
        }
    }
    candidates.sort_by(|a, b| a.3.partial_cmp(&b.3).unwrap_or(std::cmp::Ordering::Equal));
    if let Some((start_support, end_support, joint_point, _)) = candidates.into_iter().next() {
        let joint = brep.vertices.insert(Vertex {
            position: joint_point,
            outgoing: HalfEdgeHandle::default(),
            name: None,
            tolerance: Some(brep.tolerance.linear()),
        });
        let (edge_a, half_edge_a, _twin_a) = insert_unwired_edge_pair(brep, start, joint);
        let (edge_b, half_edge_b, _twin_b) = insert_unwired_edge_pair(brep, joint, end);
        let curve_a = exact_support_cylinder_intersection_curve_between(
            brep,
            cylinder_face,
            start,
            joint,
            start_point,
            joint_point,
            label,
            fillet_index,
        );
        let curve_b = exact_support_cylinder_intersection_curve_between(
            brep,
            cylinder_face,
            joint,
            end,
            joint_point,
            end_point,
            label,
            fillet_index,
        );
        match (curve_a, curve_b) {
            (Ok(curve_a), Ok(curve_b)) => {
                brep.edges[edge_a].curve = Some(curve_a);
                brep.edges[edge_b].curve = Some(curve_b);
                if crate::a2_trace::enabled() {
                    eprintln!(
                        "[a2-trace §14.9/cylinder-loop-replace] split far-side FilletEnd[{}] \
                         through joint {joint:?} using supports {:?}->{:?}",
                        fillet_index, start_support, end_support
                    );
                }
                return Ok(vec![
                    SupportOwnedTrimSide {
                        cylinder_face,
                        support_face: start_support,
                        fillet_index,
                        half_edge: half_edge_a,
                        support_side_owned: false,
                        support_source: None,
                        requires_support_ownership: true,
                    },
                    SupportOwnedTrimSide {
                        cylinder_face,
                        support_face: end_support,
                        fillet_index,
                        half_edge: half_edge_b,
                        support_side_owned: false,
                        support_source: None,
                        requires_support_ownership: true,
                    },
                ]);
            }
            (left, right) => {
                let left_status = match left {
                    Ok(_) => "ok".to_string(),
                    Err(err) => format!("{err:?}"),
                };
                let right_status = match right {
                    Ok(_) => "ok".to_string(),
                    Err(err) => format!("{err:?}"),
                };
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] found \
                         plane-plane-cylinder joint {joint:?} for {label}, but could not author \
                         both exact support/cylinder sub-spans: left {start:?}->{joint:?} \
                         via {start_support:?}: {left_status}; right {joint:?}->{end:?} via \
                         {end_support:?}: {right_status}"
                    ),
                });
            }
        }
    }
    if let Ok(curve) =
        exact_cylinder_arc_trim_curve_between(brep, cylinder_face, start, end, label, fillet_index)
    {
        let (edge, half_edge) = insert_boundary_half_edge(brep, start);
        brep.edges[edge].curve = Some(curve);
        if crate::a2_trace::enabled() {
            eprintln!(
                "[a2-trace §14.9/cylinder-loop-replace] FilletEnd[{fillet_index}] authored \
                 cylinder-only far-side boundary arc {start:?}->{end:?}; support transition \
                 has no plane-plane-cylinder joint"
            );
        }
        return Ok(vec![SupportOwnedTrimSide {
            cylinder_face,
            support_face: start_supports
                .first()
                .copied()
                .or_else(|| end_supports.first().copied())
                .unwrap_or_default(),
            fillet_index,
            half_edge,
            support_side_owned: false,
            support_source: None,
            requires_support_ownership: false,
        }]);
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] cannot split exact \
             far-side cylinder trim {start:?}->{end:?}: no plane-plane-cylinder transition \
             point from start supports {:?} to end supports {:?}",
            start_supports, end_supports
        ),
    })
}

fn unique_support_face_for_trim_endpoints(
    brep: &BRep,
    cylinder_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Result<FaceHandle, ApexBuildError> {
    let start_point = brep
        .vertices
        .get(start)
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot choose support \
                 owner for {label}: stale start vertex {start:?}"
            ),
        })?;
    let end_point = brep
        .vertices
        .get(end)
        .map(|vertex| vertex.position)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition FilletEnd[{fillet_index}] cannot choose support \
                 owner for {label}: stale end vertex {end:?}"
            ),
        })?;
    let tol = brep.tolerance.linear().max(1.0e-9);
    let mut candidates = brep
        .faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cylinder_face || !allowed_support_faces.contains(&face) {
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
            "planned repeated-apex partition FilletEnd[{fillet_index}] cannot choose unique \
             support owner for {label} {start:?}->{end:?}: candidates {:?} from allowed faces \
             {:?}",
            candidates, allowed_support_faces
        ),
    })
}
