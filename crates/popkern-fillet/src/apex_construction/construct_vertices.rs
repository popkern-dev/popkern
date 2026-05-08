//! Apex vertex resolution for the legacy A2 mutating builder.

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::fillet_surface::FilletSurfaceResult;

use super::error::ApexBuildError;
use super::plan::ApexPlan;
use super::target_topology::TargetEdge;

pub(super) fn resolve_apex_vertices(
    brep: &mut BRep,
    plan: &ApexPlan,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<Vec<VertexHandle>, ApexBuildError> {
    let tol = tolerance.linear();
    let trace = std::env::var_os("POPKERN_A2_TRACE").is_some();
    let max_radius = max_fillet_radius(fillet_surfaces).unwrap_or(0.1);
    let wide_tol = max_radius * 0.7;
    let n_targets = plan.target_topology.apex_vertices.len();
    let mut out: Vec<Option<VertexHandle>> = vec![None; n_targets];
    let mut used: std::collections::HashSet<VertexHandle> = std::collections::HashSet::new();

    for (i, &point) in plan.target_topology.apex_vertices.iter().enumerate() {
        if let Some(vh) = find_vertex_within(brep, point, tol)
            && !used.contains(&vh)
        {
            if trace {
                eprintln!(
                    "[a2-trace] apex_vertex[{i}] at {point:?} -> FOUND existing {vh:?} \
                     (within {tol:.3e})"
                );
            }
            out[i] = Some(vh);
            used.insert(vh);
        }
    }

    let apex_vertex_fillets =
        compute_apex_vertex_fillet_map(&plan.target_topology.apex_edges, n_targets);
    let fillet_face_map = build_fillet_face_map(brep, fillet_surfaces);
    type Triple = (usize, VertexHandle, f64);
    let mut triples: Vec<Triple> = Vec::new();
    let v_pos = brep
        .vertices
        .get(plan.apex_vertex)
        .map(|v| v.position)
        .unwrap_or(Point3::new(0.0, 0.0, 0.0));
    let v_radius_sq = (max_radius * 1.5).powi(2);
    let wide_tol_sq = wide_tol * wide_tol;
    for (i, &target) in plan.target_topology.apex_vertices.iter().enumerate() {
        if out[i].is_some() {
            continue;
        }
        let expected_fillet = apex_vertex_fillets[i];
        let candidate_set: Box<dyn Iterator<Item = VertexHandle>> = match expected_fillet
            .and_then(|fi| fillet_face_map.get(&fi).copied())
            .and_then(|fh| collect_face_outer_loop_vertices(brep, fh))
        {
            Some(verts) => Box::new(verts.into_iter()),
            None => Box::new(brep.vertices.keys()),
        };
        for vh in candidate_set {
            let Some(vertex) = brep.vertices.get(vh) else {
                continue;
            };
            let p = vertex.position;
            let d_v_sq = {
                let dx = p.x - v_pos.x;
                let dy = p.y - v_pos.y;
                let dz = p.z - v_pos.z;
                dx * dx + dy * dy + dz * dz
            };
            if d_v_sq > v_radius_sq {
                continue;
            }
            let d_target_sq = {
                let dx = p.x - target.x;
                let dy = p.y - target.y;
                let dz = p.z - target.z;
                dx * dx + dy * dy + dz * dz
            };
            if d_target_sq > wide_tol_sq {
                continue;
            }
            triples.push((i, vh, d_target_sq));
        }
    }
    triples.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
    for (i, vh, dist_sq) in triples {
        if out[i].is_some() || used.contains(&vh) {
            continue;
        }
        let dist = dist_sq.sqrt();
        if trace {
            eprintln!(
                "[a2-trace] apex_vertex[{i}] at {:?} -> SNAPPED to existing {vh:?} at {:?} \
                 (dist {dist:.4} m, wide_tol {wide_tol:.4})",
                plan.target_topology.apex_vertices[i], brep.vertices[vh].position
            );
        }
        out[i] = Some(vh);
        used.insert(vh);
    }

    for (i, &point) in plan.target_topology.apex_vertices.iter().enumerate() {
        if out[i].is_some() {
            continue;
        }
        let Some((edge, t)) = crate::construct::helpers::find_edge_and_param_at(brep, point) else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "apex_vertex[{i}] at {point:?}: no existing BRep vertex within \
                     {tol:.3e} m AND no edge whose chord passes through {point:?} \
                     (A2.2.b's create branch needs an edge to subdivide)"
                ),
            });
        };
        let pre_faces = if trace {
            edge_adjacent_faces(brep, edge)
        } else {
            (None, None)
        };
        let split = brep
            .split_edge(edge, t)
            .map_err(|e| ApexBuildError::PlanMismatch {
                reason: format!(
                    "apex_vertex[{i}] at {point:?}: split_edge on edge {edge:?} at \
                     t={t:.6e} failed: {e}"
                ),
            })?;
        if trace {
            let loops = vertex_loop_face_handles(brep, split.vertex);
            eprintln!(
                "[a2-trace] apex_vertex[{i}] at {point:?} -> CREATED {:?} via \
                 split_edge(edge={edge:?}, t={t:.6e}); pre-split faces=({:?}, {:?}); \
                 post-split outer-loop face memberships={loops:?}",
                split.vertex, pre_faces.0, pre_faces.1
            );
        }
        used.insert(split.vertex);
        out[i] = Some(split.vertex);
    }
    out.into_iter()
        .enumerate()
        .map(|(i, v)| {
            v.ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!("apex_vertex[{i}] unresolved after snap + fallback passes"),
            })
        })
        .collect()
}

fn find_vertex_within(brep: &BRep, point: Point3, tol: f64) -> Option<VertexHandle> {
    let mut best: Option<(VertexHandle, f64)> = None;
    for (vh, vertex) in &brep.vertices {
        let d2 = {
            let dx = vertex.position.x - point.x;
            let dy = vertex.position.y - point.y;
            let dz = vertex.position.z - point.z;
            dx * dx + dy * dy + dz * dz
        };
        if d2 <= tol * tol {
            match best {
                Some((_, b)) if b <= d2 => {}
                _ => best = Some((vh, d2)),
            }
        }
    }
    best.map(|(vh, _)| vh)
}

fn build_fillet_face_map(
    brep: &BRep,
    fillet_surfaces: &[FilletSurfaceResult],
) -> std::collections::HashMap<usize, FaceHandle> {
    use popkern_core::geometry::surface::SurfaceGeometry;
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

fn compute_apex_vertex_fillet_map(
    apex_edges: &[TargetEdge],
    n_targets: usize,
) -> Vec<Option<usize>> {
    use super::target_topology::TargetEdgeRole;
    let mut out: Vec<Option<usize>> = vec![None; n_targets];
    for edge in apex_edges {
        if let TargetEdgeRole::CapVsCyl { fillet_idx } = edge.role {
            if edge.origin_vertex_idx < n_targets {
                out[edge.origin_vertex_idx] = Some(fillet_idx);
            }
            if edge.dest_vertex_idx < n_targets {
                out[edge.dest_vertex_idx] = Some(fillet_idx);
            }
        }
    }
    out
}

fn collect_face_outer_loop_vertices(brep: &BRep, face: FaceHandle) -> Option<Vec<VertexHandle>> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    let mut out: Vec<VertexHandle> = Vec::new();
    let mut cursor = start;
    let max_iters = brep.half_edges.len() + 1;
    let mut iters = 0;
    while let Some(he) = brep.half_edges.get(cursor) {
        out.push(he.origin);
        cursor = he.next;
        iters += 1;
        if cursor == start || iters >= max_iters {
            break;
        }
    }
    Some(out)
}

fn max_fillet_radius(fillet_surfaces: &[FilletSurfaceResult]) -> Option<f64> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let mut max_r = 0.0_f64;
    let mut any = false;
    for fs in fillet_surfaces {
        let r = match &fs.surface {
            SurfaceGeometry::Cylinder { radius, .. } => *radius,
            SurfaceGeometry::Torus { minor_radius, .. } => *minor_radius,
            _ => continue,
        };
        if r.is_finite() && r > 0.0 {
            any = true;
            if r > max_r {
                max_r = r;
            }
        }
    }
    if any { Some(max_r) } else { None }
}

fn edge_adjacent_faces(brep: &BRep, edge: EdgeHandle) -> (Option<FaceHandle>, Option<FaceHandle>) {
    let Some(edata) = brep.edges.get(edge) else {
        return (None, None);
    };
    let primary = edata.half_edge;
    let Some(primary_he) = brep.half_edges.get(primary) else {
        return (None, None);
    };
    let face_a = brep.loops.get(primary_he.loop_).map(|l| l.face);
    let face_b = brep
        .half_edges
        .get(primary_he.twin)
        .and_then(|t| brep.loops.get(t.loop_))
        .map(|l| l.face);
    (face_a, face_b)
}

pub(super) fn vertex_loop_face_handles(brep: &BRep, vertex: VertexHandle) -> Vec<FaceHandle> {
    let mut out: Vec<FaceHandle> = Vec::new();
    for he in brep.half_edges.values() {
        if he.origin != vertex {
            continue;
        }
        if let Some(loop_data) = brep.loops.get(he.loop_)
            && !out.contains(&loop_data.face)
        {
            out.push(loop_data.face);
        }
    }
    out
}
