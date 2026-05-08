use super::*;

pub(super) fn materialized_generator_trim_side(
    brep: &mut BRep,
    cylinder_face: FaceHandle,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    label: &'static str,
    fillet_index: usize,
) -> Result<
    (
        HalfEdgeHandle,
        FaceHandle,
        bool,
        Option<SupportOwnedTrimSource>,
    ),
    ApexBuildError,
> {
    let curve =
        exact_cylinder_trim_curve_between(brep, cylinder_face, start, end, label, fillet_index)?;
    match materialize_support_face_mef_trim(
        brep,
        support_face,
        start,
        end,
        curve.clone(),
        label,
        fillet_index,
    ) {
        Ok((half_edge, actual_support_face, support_source)) => {
            return Ok((half_edge, actual_support_face, true, Some(support_source)));
        }
        Err(err) => {
            if crate::a2_trace::enabled() {
                eprintln!(
                    "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] cannot author \
                     {label} {start:?}->{end:?} through support-face mef on {support_face:?}: \
                     {err}"
                );
            }
        }
    }
    let (edge, half_edge, _twin) = insert_unwired_edge_pair(brep, start, end);
    brep.edges[edge].curve = Some(curve);
    Ok((half_edge, support_face, false, None))
}

fn materialize_support_face_mef_trim(
    brep: &mut BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
    curve: CurveOnEdge,
    label: &'static str,
    fillet_index: usize,
) -> Result<(HalfEdgeHandle, FaceHandle, SupportOwnedTrimSource), ApexBuildError> {
    let rollback = brep.clone();
    let result = (|| {
        let (actual_support_face, support_start, support_end) =
            support_residual_face_split_half_edges(brep, support_face, start, end).ok_or_else(
                || ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                         {label} {start:?}->{end:?} by support-face mef: no residual support \
                         face has both endpoints on one outer loop. {}",
                        support_family_endpoint_diagnostic(brep, support_face, start, end)
                    ),
                },
            )?;
        let result =
            brep.mef(support_start, support_end)
                .map_err(|err| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "planned repeated-apex partition FilletEnd[{fillet_index}] cannot author \
                     {label} {start:?}->{end:?} by support-face mef on {actual_support_face:?}: \
                     {err:?}"
                    ),
                })?;
        brep.edges[result.edge].curve = Some(curve);
        let support_source = SupportOwnedTrimSource {
            loop_: brep.half_edges[result.half_edge].loop_,
            prev: brep.half_edges[result.half_edge].prev,
            next: brep.half_edges[result.half_edge].next,
        };
        if crate::a2_trace::enabled() {
            eprintln!(
                "[a2-trace §14.9/trim-subcomplex] FilletEnd[{fillet_index}] authored {label} \
                 {start:?}->{end:?} by support-face mef on {actual_support_face:?}; cylinder \
                 side {:?}, support side {:?}",
                result.half_edge, result.half_edge_twin
            );
        }
        Ok((result.half_edge, actual_support_face, support_source))
    })();
    if result.is_err() {
        *brep = rollback;
    }
    result
}

fn support_family_endpoint_diagnostic(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> String {
    let mut candidate_faces = Vec::new();
    candidate_faces.push(support_face);
    candidate_faces.extend(coplanar_residual_support_faces(brep, support_face));
    let reports = candidate_faces
        .into_iter()
        .map(|face| {
            let start_on_loop = vertex_is_on_face_outer_loop(brep, face, start);
            let end_on_loop = vertex_is_on_face_outer_loop(brep, face, end);
            let start_edge = brep.vertices.get(start).and_then(|vertex| {
                containing_outer_loop_edge_by_curve(brep, face, vertex.position)
            });
            let end_edge = brep.vertices.get(end).and_then(|vertex| {
                containing_outer_loop_edge_by_curve(brep, face, vertex.position)
            });
            format!(
                "{face:?}: start_on_loop={start_on_loop} end_on_loop={end_on_loop} \
                 start_edge={start_edge:?} end_edge={end_edge:?}"
            )
        })
        .collect::<Vec<_>>()
        .join("; ");
    format!("support family endpoint status [{reports}]")
}

fn support_residual_face_split_half_edges(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<(FaceHandle, HalfEdgeHandle, HalfEdgeHandle)> {
    let mut candidate_faces = Vec::new();
    candidate_faces.push(support_face);
    candidate_faces.extend(coplanar_residual_support_faces(brep, support_face));
    candidate_faces.into_iter().find_map(|candidate_face| {
        support_face_split_half_edges(brep, candidate_face, start, end)
            .map(|(support_start, support_end)| (candidate_face, support_start, support_end))
    })
}

fn support_face_split_half_edges(
    brep: &BRep,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<(HalfEdgeHandle, HalfEdgeHandle)> {
    let face = brep.faces.get(support_face)?;
    let loop_handle = face.outer_loop;
    let loop_data = brep.loops.get(loop_handle)?;
    let mut start_edges = Vec::new();
    let mut end_edges = Vec::new();
    let mut cursor = loop_data.half_edge;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.loop_ != loop_handle {
            return None;
        }
        if half_edge.origin == start {
            start_edges.push(cursor);
        }
        if half_edge.origin == end {
            end_edges.push(cursor);
        }
        cursor = half_edge.next;
        if cursor == loop_data.half_edge {
            break;
        }
    }
    match (start_edges.as_slice(), end_edges.as_slice()) {
        ([support_start], [support_end]) if support_start != support_end => {
            Some((*support_start, *support_end))
        }
        _ => None,
    }
}

fn coplanar_residual_support_faces(brep: &BRep, support_face: FaceHandle) -> Vec<FaceHandle> {
    let Some(source_face) = brep.faces.get(support_face) else {
        return Vec::new();
    };
    let Some(SurfaceGeometry::Plane {
        origin: source_origin,
        normal: source_normal,
        ..
    }) = source_face
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)
    else {
        return Vec::new();
    };
    let tol = brep.tolerance.linear().max(1.0e-9);
    brep.faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == support_face || face_data.shell != source_face.shell {
                return None;
            }
            let SurfaceGeometry::Plane { origin, normal, .. } = face_data
                .surface
                .as_ref()
                .map(|surface| &surface.geometry)?
            else {
                return None;
            };
            let normals_parallel = source_normal.cross(normal).length() <= tol;
            let same_plane = (*origin - *source_origin).dot(source_normal).abs() <= tol;
            (normals_parallel && same_plane).then_some(face)
        })
        .collect()
}
