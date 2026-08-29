use super::*;

pub(super) fn validate_replacement_loop_support_side_ownership(
    brep: &BRep,
    cylinder_face: FaceHandle,
    fillet_index: usize,
    cycle: &[HalfEdgeHandle],
) -> Result<(), ApexBuildError> {
    let mut missing = Vec::new();
    for half_edge in cycle {
        let Some(twin) = brep.half_edges.get(*half_edge).map(|data| data.twin) else {
            missing.push(format!("{half_edge:?}: stale cylinder half-edge"));
            continue;
        };
        let Some(twin_data) = brep.half_edges.get(twin) else {
            missing.push(format!("{half_edge:?}: stale twin {twin:?}"));
            continue;
        };
        let Some(twin_loop) = brep.loops.get(twin_data.loop_) else {
            missing.push(format!(
                "{half_edge:?}/twin {twin:?}: no support-side loop owns the new trim side; {}",
                replacement_trim_support_split_diagnostic(brep, cylinder_face, *half_edge)
            ));
            continue;
        };
        if twin_loop.face == cylinder_face {
            missing.push(format!(
                "{half_edge:?}/twin {twin:?}: twin is still owned by target cylinder face \
                 {cylinder_face:?}"
            ));
        }
    }
    if missing.is_empty() {
        return Ok(());
    }
    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "planned repeated-apex partition FilletEnd[{fillet_index}] cannot accept \
             materialized target-cylinder loop yet: support-side ownership is incomplete for \
             replacement trim side(s): {}. The next construction slice must split or create the \
             adjacent support loops/faces that own these exact cylinder trims before exposing the \
             rewritten cylinder face.",
            missing.join("; ")
        ),
    })
}

fn replacement_trim_support_split_diagnostic(
    brep: &BRep,
    cylinder_face: FaceHandle,
    half_edge: HalfEdgeHandle,
) -> String {
    let Some(start) = brep.half_edges.get(half_edge).map(|data| data.origin) else {
        return "stale trim start prevents support split planning".to_string();
    };
    let Some(end) = half_edge_end_vertex(brep, half_edge) else {
        return "stale trim end prevents support split planning".to_string();
    };
    let Some(start_point) = brep.vertices.get(start).map(|vertex| vertex.position) else {
        return format!("stale trim start vertex {start:?} prevents support split planning");
    };
    let Some(end_point) = brep.vertices.get(end).map(|vertex| vertex.position) else {
        return format!("stale trim end vertex {end:?} prevents support split planning");
    };
    let tol = brep.tolerance.linear().max(1.0e-9);
    let mut support_faces = Vec::new();
    for (face, face_data) in brep.faces.iter() {
        if face == cylinder_face {
            continue;
        }
        let Some(SurfaceGeometry::Plane { origin, normal, .. }) =
            face_data.surface.as_ref().map(|surface| &surface.geometry)
        else {
            continue;
        };
        let start_distance = (start_point - *origin).dot(normal).abs();
        let end_distance = (end_point - *origin).dot(normal).abs();
        if start_distance > tol || end_distance > tol {
            continue;
        }
        let start_on_boundary = vertex_is_on_face_outer_loop(brep, face, start);
        let end_on_boundary = vertex_is_on_face_outer_loop(brep, face, end);
        let start_edge = containing_outer_loop_edge_by_curve(brep, face, start_point);
        let end_edge = containing_outer_loop_edge_by_curve(brep, face, end_point);
        support_faces.push(format!(
            "{face:?}: start_boundary={start_on_boundary}, end_boundary={end_on_boundary}, \
             start_edge={start_edge:?}, end_edge={end_edge:?}"
        ));
    }
    if support_faces.is_empty() {
        return format!(
            "no real support plane contains trim endpoints {start:?}->{end:?}; construction must \
             refuse rather than inventing a support owner"
        );
    }
    format!(
        "candidate support split face(s) [{}]; anchors that are neither existing boundary vertices \
         nor splittable boundary-edge points require an interior support-face subdivision before \
         the cylinder loop can be exposed",
        support_faces.join("; ")
    )
}
