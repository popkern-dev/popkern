//! Apex vertex authoring for constructive apex building.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve_inversion::invert_curve_at_point;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::fillet_surface::FilletSurfaceResult;

use super::error::ApexBuildError;
use super::plan::ApexPlan;

pub(super) fn author_sphere_apex_vertices(
    brep: &mut BRep,
    plan: &ApexPlan,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<Vec<VertexHandle>, ApexBuildError> {
    let tol = tolerance.linear();
    let trace = crate::a2_trace::enabled();
    let mut out: Vec<VertexHandle> = Vec::with_capacity(plan.target_topology.apex_vertices.len());
    let mut used: std::collections::HashSet<VertexHandle> = std::collections::HashSet::new();
    for (i, target) in plan.target_topology.apex_vertices.iter().enumerate() {
        match find_vertex_within(brep, *target, tol) {
            Some(vh) if !used.contains(&vh) => {
                if trace {
                    eprintln!("[a2-trace §14.4/A] sphere apex_vertices[{i}] tight-match {vh:?}");
                }
                used.insert(vh);
                out.push(vh);
            }
            Some(vh) => {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "author_sphere_apex_vertices: apex_vertices[{i}] {target:?} tight-matched \
                         already-used vertex {vh:?}"
                    ),
                });
            }
            None => {
                let nearest = nearest_vertex_distance(brep, *target);
                let split = split_sphere_target_vertex(
                    brep,
                    plan,
                    fillet_surfaces,
                    *target,
                    tolerance,
                )
                .map_err(|reason| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "author_sphere_apex_vertices: apex_vertices[{i}] {target:?} not within \
                                 {tol:e} m of any existing vertex (nearest {nearest:.6e} m), \
                                 and curve-inversion split failed: {reason}"
                    ),
                })?;
                if trace {
                    eprintln!(
                        "[a2-trace §14.4/A] sphere apex_vertices[{i}] -> CREATED {:?} via \
                         split_edge_at_point(edge={:?})",
                        split.vertex, split.edge_a
                    );
                }
                used.insert(split.vertex);
                out.push(split.vertex);
            }
        }
    }
    Ok(out)
}

fn split_sphere_target_vertex(
    brep: &mut BRep,
    plan: &ApexPlan,
    fillet_surfaces: &[FilletSurfaceResult],
    target: Point3,
    tolerance: &Tolerance,
) -> Result<popkern_core::topology::SplitEdgeResult, String> {
    let candidates = collect_sphere_target_candidate_edges(brep, fillet_surfaces)?;
    let mut matches: Vec<EdgeHandle> = Vec::new();
    for edge in candidates {
        let Some(edge_data) = brep.edges.get(edge) else {
            continue;
        };
        let Some(curve) = edge_data.curve.as_ref() else {
            continue;
        };
        let Ok(t_curve) = invert_curve_at_point(curve, target, tolerance) else {
            continue;
        };
        let point_on_curve = curve.geometry.point_at(t_curve);
        if point_on_curve.distance_to(&target) <= tolerance.linear() {
            matches.push(edge);
        }
    }
    matches.sort_by_key(|edge| format!("{edge:?}"));
    matches.dedup();
    match matches.as_slice() {
        [] => Err(format!(
            "no fillet cylinder outer-loop edge contains target for apex {:?}",
            plan.apex_vertex
        )),
        [edge] => brep
            .split_edge_at_point(*edge, target, tolerance)
            .map_err(|e| format!("split_edge_at_point({edge:?}): {e}")),
        many => Err(format!(
            "{} candidate fillet cylinder edges contain target: {many:?}",
            many.len()
        )),
    }
}

fn collect_sphere_target_candidate_edges(
    brep: &BRep,
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<Vec<EdgeHandle>, String> {
    let fillet_faces = build_fillet_face_map(brep, fillet_surfaces);
    if fillet_faces.len() != fillet_surfaces.len() {
        return Err(format!(
            "matched {} of {} fillet cylinder faces",
            fillet_faces.len(),
            fillet_surfaces.len()
        ));
    }
    let mut edges: Vec<EdgeHandle> = Vec::new();
    for face in fillet_faces.values().copied() {
        collect_face_outer_loop_edges(brep, face, &mut edges)?;
    }
    edges.sort_by_key(|edge| format!("{edge:?}"));
    edges.dedup();
    Ok(edges)
}

fn collect_face_outer_loop_edges(
    brep: &BRep,
    face: FaceHandle,
    out: &mut Vec<EdgeHandle>,
) -> Result<(), String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} missing"))?;
    let loop_data = brep
        .loops
        .get(face_data.outer_loop)
        .ok_or_else(|| format!("outer loop {:?} missing", face_data.outer_loop))?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;
    loop {
        let Some(he) = brep.half_edges.get(cursor) else {
            return Err(format!(
                "half-edge {cursor:?} missing in face {face:?} outer loop"
            ));
        };
        out.push(he.edge);
        cursor = he.next;
        iters += 1;
        if cursor == start {
            return Ok(());
        }
        if iters > max_iters {
            return Err(format!("face {face:?} outer loop did not close"));
        }
    }
}

fn build_fillet_face_map(
    brep: &BRep,
    fillet_surfaces: &[FilletSurfaceResult],
) -> std::collections::HashMap<usize, FaceHandle> {
    let mut out: std::collections::HashMap<usize, FaceHandle> = std::collections::HashMap::new();
    for (fi, fillet) in fillet_surfaces.iter().enumerate() {
        let SurfaceGeometry::Cylinder {
            origin: f_origin,
            axis: f_axis,
            radius: f_radius,
        } = &fillet.surface
        else {
            continue;
        };
        let Some(f_axis_unit) = f_axis.normalized() else {
            continue;
        };
        for (fh, face) in &brep.faces {
            let Some(sof) = face.surface.as_ref() else {
                continue;
            };
            let SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } = &sof.geometry
            else {
                continue;
            };
            if (*radius - *f_radius).abs() > f_radius.abs() * 1e-4 + 1e-9 {
                continue;
            }
            let Some(face_axis_unit) = axis.normalized() else {
                continue;
            };
            if face_axis_unit.dot(&f_axis_unit).abs() < 1.0 - 1e-3 {
                continue;
            }
            let d = popkern_core::geometry::Vec3::new(
                origin.x - f_origin.x,
                origin.y - f_origin.y,
                origin.z - f_origin.z,
            );
            let along = d.dot(&f_axis_unit);
            let perp = popkern_core::geometry::Vec3::new(
                d.x - along * f_axis_unit.x,
                d.y - along * f_axis_unit.y,
                d.z - along * f_axis_unit.z,
            );
            if perp.length() > 1e-3 {
                continue;
            }
            out.insert(fi, fh);
            break;
        }
    }
    out
}

fn nearest_vertex_distance(brep: &BRep, point: Point3) -> f64 {
    let mut best = f64::INFINITY;
    for (_, vertex) in &brep.vertices {
        let d = vertex.position.distance_to(&point);
        if d < best {
            best = d;
        }
    }
    best
}

fn find_vertex_within(brep: &BRep, point: Point3, tol: f64) -> Option<VertexHandle> {
    let mut best: Option<(VertexHandle, f64)> = None;
    let tol_sq = tol * tol;
    for (vh, vertex) in &brep.vertices {
        let diff = vertex.position - point;
        let d2 = diff.dot(&diff);
        if d2 <= tol_sq {
            match best {
                Some((_, b)) if b <= d2 => {}
                _ => best = Some((vh, d2)),
            }
        }
    }
    best.map(|(vh, _)| vh)
}

pub(super) fn author_gregory_apex_vertices(
    brep: &mut BRep,
    plan: &ApexPlan,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<Vec<VertexHandle>, ApexBuildError> {
    use crate::apex_construction::plan_terminations::compute_apex_boundary;
    use crate::blend::boundary::BoundarySegment;

    let n = plan.cap_kind.n_fillets();
    let n_targets = plan.target_topology.apex_vertices.len();
    if n_targets != 2 * n {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "author_gregory_apex_vertices: expected 2N={} apex_vertices for Gregory, got {n_targets}",
                2 * n
            ),
        });
    }

    let boundary = compute_apex_boundary(brep, plan.apex_vertex, fillet_surfaces, tolerance)
        .map_err(|e| ApexBuildError::PlanMismatch {
            reason: format!("author_gregory_apex_vertices: compute_apex_boundary: {e}"),
        })?;

    if boundary.segments.len() != 2 * n {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "author_gregory_apex_vertices: expected {} boundary segments, got {}",
                2 * n,
                boundary.segments.len()
            ),
        });
    }

    let mut face_a: Vec<FaceHandle> = vec![Default::default(); n];
    let mut face_b: Vec<FaceHandle> = vec![Default::default(); n];
    for k in 0..n {
        let face_b_k = match &boundary.segments[2 * k + 1] {
            BoundarySegment::FaceEdge { face, .. } => *face,
            BoundarySegment::FilletEnd { .. } => {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "author_gregory_apex_vertices: boundary segment[{}] expected FaceEdge, got FilletEnd",
                        2 * k + 1
                    ),
                });
            }
        };
        face_b[k] = face_b_k;
        face_a[(k + 1) % n] = face_b_k;
    }

    let trace = crate::a2_trace::enabled();
    let mut out: Vec<Option<VertexHandle>> = vec![None; n_targets];
    let mut t_vertices: std::collections::HashSet<VertexHandle> = std::collections::HashSet::new();
    let mut schedule: Vec<(usize, FaceHandle)> = Vec::with_capacity(n_targets);
    for (k, &fb) in face_b.iter().enumerate() {
        schedule.push((2 * k + 1, fb));
    }
    for (k, &fa) in face_a.iter().enumerate() {
        schedule.push((2 * k, fa));
    }
    for (sched_idx, &(i, host_face)) in schedule.iter().enumerate() {
        let target = plan.target_topology.apex_vertices[i];

        let he_v =
            find_v_cube_edge_half_edge_on_face(brep, plan.apex_vertex, host_face, &t_vertices)
                .ok_or_else(|| ApexBuildError::PlanMismatch {
                    reason: format!(
                        "author_gregory_apex_vertices: no V-incident cube-edge half-edge on host face {host_face:?} for apex_vertices[{i}] (sched_idx={sched_idx})"
                    ),
                })?;

        let mev = brep
            .mev(he_v, target)
            .map_err(|e| ApexBuildError::EulerOpFailed {
                op: "mev",
                step: i,
                cause: format!("{e}"),
            })?;
        if trace {
            eprintln!(
                "[a2-trace §14.8/A] sched[{sched_idx}] mev V->{:?} at {target:?} on {host_face:?} (i={i}, k={}, side={})",
                mev.vertex,
                i / 2,
                if i % 2 == 0 { "face_A" } else { "face_B" }
            );
        }
        t_vertices.insert(mev.vertex);
        out[i] = Some(mev.vertex);
    }
    let out: Vec<VertexHandle> = out
        .into_iter()
        .enumerate()
        .map(|(i, v)| {
            v.ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!("author_gregory_apex_vertices: apex_vertices[{i}] not assigned"),
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    if trace {
        eprintln!(
            "[a2-trace §14.8/A] author_gregory_apex_vertices complete: {} T-vertices authored",
            out.len()
        );
    }

    Ok(out)
}

fn find_v_cube_edge_half_edge_on_face(
    brep: &BRep,
    v: VertexHandle,
    face: FaceHandle,
    t_vertices: &std::collections::HashSet<VertexHandle>,
) -> Option<HalfEdgeHandle> {
    let outer_loop = brep.faces.get(face)?.outer_loop;
    let start = brep.loops.get(outer_loop)?.half_edge;
    let mut cursor = start;
    let cap = brep.half_edges.len() + 1;
    let mut iters = 0;
    loop {
        let he = brep.half_edges.get(cursor)?;
        if he.origin == v {
            let dest = brep.half_edges.get(he.next)?.origin;
            if !t_vertices.contains(&dest) {
                return Some(cursor);
            }
        }
        cursor = he.next;
        iters += 1;
        if cursor == start || iters > cap {
            return None;
        }
    }
}
