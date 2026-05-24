use super::*;

pub(super) fn selected_far_side_trim_materialization_preflight(
    brep: &BRep,
    cap_face: FaceHandle,
    face_edges: &[FaceEdgeOwnership],
    rebuild_plans: &[SemanticCylinderTrimRebuildPlan],
    topology_plans: &[CylinderTrimRewriteTopologyPlan],
) -> (bool, Vec<crate::types::ApexConstructionDiagnostic>) {
    let mut candidate = brep.clone();
    match pre_materialize_far_side_trim_sides(
        &mut candidate,
        cap_face,
        face_edges,
        rebuild_plans,
        topology_plans,
        true,
    ) {
        Ok(_) => (true, Vec::new()),
        Err(error) => {
            let diagnostics = materialization_preflight_diagnostics(
                brep,
                cap_face,
                face_edges,
                rebuild_plans,
                topology_plans,
                &error.to_string(),
            );
            (false, diagnostics)
        }
    }
}

fn materialization_preflight_diagnostics(
    brep: &BRep,
    cap_face: FaceHandle,
    face_edges: &[FaceEdgeOwnership],
    rebuild_plans: &[SemanticCylinderTrimRebuildPlan],
    topology_plans: &[CylinderTrimRewriteTopologyPlan],
    error: &str,
) -> Vec<crate::types::ApexConstructionDiagnostic> {
    rebuild_plans
        .iter()
        .zip(topology_plans)
        .filter_map(|(rebuild_plan, topology_plan)| {
            rebuild_plan.fragment_for_strategy(topology_plan.loop_candidate.far_side_strategy)?;
            let far_side = topology_plan
                .loop_candidate
                .segments
                .iter()
                .find(|segment| segment.label.contains("far-side"))?;
            let allowed_support_faces = target_cylinder_support_faces(
                brep,
                rebuild_plan.cylinder_face,
                cap_face,
                face_edges,
            );
            let start_support_faces =
                support_faces_containing_vertex(brep, far_side.start, &allowed_support_faces);
            let end_support_faces =
                support_faces_containing_vertex(brep, far_side.end, &allowed_support_faces);
            let reason = if error.contains("no plane-plane-cylinder transition point") {
                crate::types::CylinderTrimMaterializationDeclineReason::NoPlanePlaneCylinderTransition
            } else {
                crate::types::CylinderTrimMaterializationDeclineReason::Other
            };
            let strategy = match topology_plan.loop_candidate.far_side_strategy {
                FarSideTrimStrategy::Forward => crate::types::CylinderTrimFarSideStrategy::Forward,
                FarSideTrimStrategy::Reverse => crate::types::CylinderTrimFarSideStrategy::Reverse,
            };
            Some(
                crate::types::ApexConstructionDiagnostic::CylinderTrimMaterialization(
                    crate::types::CylinderTrimMaterializationDiagnostic {
                        fillet_index: topology_plan.fillet_index,
                        cylinder_face: topology_plan.cylinder_face,
                        start: far_side.start,
                        end: far_side.end,
                        strategy,
                        start_support_faces,
                        end_support_faces,
                        reason,
                        exact_authorable: false,
                    },
                ),
            )
        })
        .collect()
}

fn support_faces_containing_vertex(
    brep: &BRep,
    vertex: VertexHandle,
    allowed_support_faces: &std::collections::BTreeSet<FaceHandle>,
) -> Vec<FaceHandle> {
    let Some(point) = brep.vertices.get(vertex).map(|vertex| vertex.position) else {
        return Vec::new();
    };
    let tolerance = brep.tolerance.linear().max(1.0e-9);
    let mut faces = allowed_support_faces
        .iter()
        .copied()
        .filter(|face| {
            let Some(surface) = brep.faces.get(*face).and_then(|face| face.surface.as_ref()) else {
                return false;
            };
            let SurfaceGeometry::Plane { origin, normal, .. } = surface.geometry else {
                return false;
            };
            (point - origin).dot(&normal).abs() <= tolerance
        })
        .collect::<Vec<_>>();
    faces.sort();
    faces.dedup();
    faces
}

pub(super) fn plan_semantic_cylinder_trim_rebuilds(
    brep: &BRep,
    fillet_ends: &[FilletEndOwnership],
) -> Result<Vec<SemanticCylinderTrimRebuildPlan>, ApexBuildError> {
    fillet_ends
        .iter()
        .map(|fillet_end| SemanticCylinderTrimRebuildPlan::build(brep, fillet_end))
        .collect()
}
