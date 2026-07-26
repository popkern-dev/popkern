use super::*;

pub(super) fn retired_trim_runs_by_adjacent_face(
    brep: &BRep,
    cylinder_face: FaceHandle,
    path: &[HalfEdgeHandle],
) -> Vec<RetiredCylinderTrimRun> {
    let mut runs = Vec::new();
    let mut current_face = None;
    let mut current_path = Vec::new();
    for half_edge in path {
        let Some(adjacent_face) = adjacent_face_across_half_edge(brep, cylinder_face, *half_edge)
        else {
            continue;
        };
        if current_face.is_some_and(|face| face != adjacent_face)
            && !current_path.is_empty()
            && let Some(face) = current_face
        {
            runs.push(RetiredCylinderTrimRun::from_path(
                brep,
                cylinder_face,
                face,
                std::mem::take(&mut current_path),
            ));
        }
        current_face = Some(adjacent_face);
        current_path.push(*half_edge);
    }
    if let Some(adjacent_face) = current_face
        && !current_path.is_empty()
    {
        runs.push(RetiredCylinderTrimRun::from_path(
            brep,
            cylinder_face,
            adjacent_face,
            current_path,
        ));
    }
    runs
}

impl RetiredCylinderTrimRun {
    fn from_path(
        brep: &BRep,
        cylinder_face: FaceHandle,
        adjacent_face: FaceHandle,
        path: Vec<HalfEdgeHandle>,
    ) -> Self {
        let start = path
            .first()
            .and_then(|half_edge| brep.half_edges.get(*half_edge).map(|data| data.origin))
            .unwrap_or_default();
        let end = path
            .last()
            .and_then(|half_edge| half_edge_end_vertex(brep, *half_edge))
            .unwrap_or_default();
        let vertices = path
            .iter()
            .filter_map(|half_edge| brep.half_edges.get(*half_edge).map(|data| data.origin))
            .collect();
        let off_surface_vertices = path_off_cylinder_vertices(brep, cylinder_face, &path);
        Self {
            adjacent_face,
            path,
            start,
            end,
            vertices,
            off_surface_vertices,
        }
    }
}

pub(super) fn adjacent_face_across_half_edge(
    brep: &BRep,
    cylinder_face: FaceHandle,
    half_edge: HalfEdgeHandle,
) -> Option<FaceHandle> {
    let twin = brep.half_edges.get(half_edge)?.twin;
    let twin_loop = brep.half_edges.get(twin)?.loop_;
    let twin_face = brep.loops.get(twin_loop)?.face;
    (twin_face != cylinder_face).then_some(twin_face)
}

pub(super) fn face_surface_kind(brep: &BRep, face: FaceHandle) -> &'static str {
    match brep
        .faces
        .get(face)
        .and_then(|face_data| face_data.surface.as_ref())
        .map(|surface| &surface.geometry)
    {
        Some(SurfaceGeometry::Plane { .. }) => "Plane",
        Some(SurfaceGeometry::Cylinder { .. }) => "Cylinder",
        Some(SurfaceGeometry::Cone { .. }) => "Cone",
        Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
        Some(SurfaceGeometry::Torus { .. }) => "Torus",
        Some(SurfaceGeometry::Nurbs(_)) => "Nurbs",
        Some(SurfaceGeometry::Gregory(_)) => "Gregory",
        Some(SurfaceGeometry::PolyhedralSpline(_)) => "PolyhedralSpline",
        None => "None",
    }
}

pub(super) fn replacement_run_exact_curve_kind(
    brep: &BRep,
    cylinder_face: FaceHandle,
    run: &RetiredCylinderTrimRun,
) -> Option<&'static str> {
    let support_surface = brep
        .faces
        .get(run.adjacent_face)?
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)?;
    let SurfaceGeometry::Plane {
        origin: plane_origin,
        normal: plane_normal,
        ..
    } = support_surface
    else {
        return None;
    };
    let start_point = brep.vertices.get(run.start)?.position;
    let end_point = brep.vertices.get(run.end)?.position;
    let start_plane_distance = (start_point - *plane_origin).dot(plane_normal).abs();
    let end_plane_distance = (end_point - *plane_origin).dot(plane_normal).abs();
    let tol = brep.tolerance.linear().max(1.0e-9);
    if start_plane_distance > tol || end_plane_distance > tol {
        return None;
    }
    if !point_lies_on_cylinder_face(brep, cylinder_face, start_point)
        || !point_lies_on_cylinder_face(brep, cylinder_face, end_point)
    {
        return None;
    }
    let cylinder = CylinderTrimSurface::from_face(brep, cylinder_face, 0).ok()?;
    if cylinder
        .validate_generator_segment(brep, 0, run.start, run.end)
        .is_ok()
    {
        return Some("Line");
    }
    let SurfaceGeometry::Cylinder { axis, .. } = brep
        .faces
        .get(cylinder_face)?
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)?
    else {
        return None;
    };
    if plane_normal.dot(axis).abs() > 1.0e-9 {
        return Some("Ellipse");
    }
    None
}

pub(super) fn fragment_replacement_authoring_status(
    brep: &BRep,
    cylinder_face: FaceHandle,
    fragment: &RetiredCylinderTrimFragment,
) -> FragmentAuthoringStatus {
    let mut labels = Vec::new();
    let mut authorable = true;
    for run in &fragment.adjacent_face_runs {
        if !run.requires_replacement() {
            labels.push(format!(
                "{:?} {:?}->{:?}: reuse",
                run.adjacent_face, run.start, run.end
            ));
            continue;
        }
        if let Some(kind) = replacement_run_exact_curve_kind(brep, cylinder_face, run) {
            labels.push(format!(
                "{:?} {:?}->{:?}: exact {kind}",
                run.adjacent_face, run.start, run.end
            ));
            continue;
        }
        if let Some(plan) = replacement_run_reownership_plan(brep, cylinder_face, run) {
            labels.push(format!(
                "{:?} {:?}->{:?}: split/re-own {}",
                run.adjacent_face,
                run.start,
                run.end,
                plan.diagnostic_label()
            ));
            continue;
        }
        authorable = false;
        labels.push(format!(
            "{:?} {:?}->{:?}: blocked ({})",
            run.adjacent_face,
            run.start,
            run.end,
            replacement_run_curve_blocker(brep, cylinder_face, run)
        ));
    }
    FragmentAuthoringStatus {
        authorable,
        diagnostic: labels.join("; "),
    }
}

pub(super) fn replacement_run_reownership_plan(
    brep: &BRep,
    cylinder_face: FaceHandle,
    run: &RetiredCylinderTrimRun,
) -> Option<ReplacementRunReownershipPlan> {
    if run.path.len() <= 1 {
        return None;
    }
    let mut spans = Vec::new();
    for half_edge in &run.path {
        let start = brep.half_edges.get(*half_edge)?.origin;
        let end = half_edge_end_vertex(brep, *half_edge)?;
        let candidate_faces = replacement_span_support_candidates(brep, cylinder_face, start, end);
        if candidate_faces.is_empty() {
            return None;
        }
        spans.push(ReplacementRunSpanOwnership {
            start,
            end,
            candidate_faces,
        });
    }
    Some(ReplacementRunReownershipPlan { spans })
}

pub(super) fn replacement_run_reownership_blocker(
    brep: &BRep,
    cylinder_face: FaceHandle,
    run: &RetiredCylinderTrimRun,
) -> String {
    if run.path.len() <= 1 {
        return "run has one span, so split/re-own decomposition does not apply".to_string();
    }
    let span_reports = run
        .path
        .iter()
        .filter_map(|half_edge| {
            let start = brep.half_edges.get(*half_edge)?.origin;
            let end = half_edge_end_vertex(brep, *half_edge)?;
            let candidates = replacement_span_support_candidates(brep, cylinder_face, start, end);
            let candidate_label = if candidates.is_empty() {
                "no exact support-plane candidates".to_string()
            } else {
                candidates
                    .iter()
                    .map(|candidate| format!("{:?}/{}", candidate.face, candidate.curve_kind))
                    .collect::<Vec<_>>()
                    .join("|")
            };
            Some(format!("{start:?}->{end:?}: {candidate_label}"))
        })
        .collect::<Vec<_>>()
        .join(", ");
    format!("reownership spans [{span_reports}]")
}

fn replacement_span_support_candidates(
    brep: &BRep,
    cylinder_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Vec<ReplacementRunSupportCandidate> {
    brep.faces
        .iter()
        .filter_map(|(face, face_data)| {
            if face == cylinder_face {
                return None;
            }
            let SurfaceGeometry::Plane { .. } = face_data.surface.as_ref()?.geometry else {
                return None;
            };
            let curve_kind =
                replacement_span_exact_curve_kind(brep, cylinder_face, face, start, end)?;
            Some(ReplacementRunSupportCandidate { face, curve_kind })
        })
        .collect()
}

fn replacement_span_exact_curve_kind(
    brep: &BRep,
    cylinder_face: FaceHandle,
    support_face: FaceHandle,
    start: VertexHandle,
    end: VertexHandle,
) -> Option<&'static str> {
    let support_surface = brep
        .faces
        .get(support_face)?
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)?;
    let SurfaceGeometry::Plane {
        origin: plane_origin,
        normal: plane_normal,
        ..
    } = support_surface
    else {
        return None;
    };
    let start_point = brep.vertices.get(start)?.position;
    let end_point = brep.vertices.get(end)?.position;
    let tol = brep.tolerance.linear().max(1.0e-9);
    if (start_point - *plane_origin).dot(plane_normal).abs() > tol
        || (end_point - *plane_origin).dot(plane_normal).abs() > tol
    {
        return None;
    }
    if !point_lies_on_cylinder_face(brep, cylinder_face, start_point)
        || !point_lies_on_cylinder_face(brep, cylinder_face, end_point)
    {
        return None;
    }
    let cylinder = CylinderTrimSurface::from_face(brep, cylinder_face, 0).ok()?;
    if cylinder
        .validate_generator_segment(brep, 0, start, end)
        .is_ok()
    {
        return Some("Line");
    }
    let SurfaceGeometry::Cylinder { axis, .. } = brep
        .faces
        .get(cylinder_face)?
        .surface
        .as_ref()
        .map(|surface| &surface.geometry)?
    else {
        return None;
    };
    if plane_normal.dot(axis).abs() > 1.0e-9 {
        return Some("Ellipse");
    }
    None
}

pub(super) fn replacement_run_curve_blocker(
    brep: &BRep,
    cylinder_face: FaceHandle,
    run: &RetiredCylinderTrimRun,
) -> String {
    let Some(start_point) = brep.vertices.get(run.start).map(|vertex| vertex.position) else {
        return format!("start vertex {:?} is stale", run.start);
    };
    let Some(end_point) = brep.vertices.get(run.end).map(|vertex| vertex.position) else {
        return format!("end vertex {:?} is stale", run.end);
    };
    let Some(SurfaceGeometry::Plane {
        origin: plane_origin,
        normal: plane_normal,
        ..
    }) = brep
        .faces
        .get(run.adjacent_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return "support surface is not a plane".to_string();
    };
    let Some(SurfaceGeometry::Cylinder {
        origin: cylinder_origin,
        axis,
        ..
    }) = brep
        .faces
        .get(cylinder_face)
        .and_then(|face| face.surface.as_ref())
        .map(|surface| &surface.geometry)
    else {
        return "target face is not a cylinder".to_string();
    };
    let plane_start = (start_point - *plane_origin).dot(plane_normal).abs();
    let plane_end = (end_point - *plane_origin).dot(plane_normal).abs();
    let plane_axis_dot = plane_normal.dot(axis).abs();
    let radial_alignment = match (
        cylinder_radial_vector(start_point, *cylinder_origin, *axis)
            .and_then(|radial| radial.normalized()),
        cylinder_radial_vector(end_point, *cylinder_origin, *axis)
            .and_then(|radial| radial.normalized()),
    ) {
        (Some(start_radial), Some(end_radial)) => start_radial.dot(&end_radial),
        _ => f64::NAN,
    };
    format!(
        "plane distances start={plane_start:.3e}, end={plane_end:.3e}, \
         |plane_normal.cylinder_axis|={plane_axis_dot:.3e}, radial_alignment={radial_alignment:.6}"
    )
}
