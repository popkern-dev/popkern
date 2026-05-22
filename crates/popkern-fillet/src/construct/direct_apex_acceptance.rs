//! Native Gregory acceptance checks for DirectApex dispatch.

use super::*;
use popkern_core::geometry::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;

pub(super) fn native_gregory_path_a_error(brep: &BRep) -> Option<String> {
    use popkern_core::error::ValidationError;

    brep.validate_all().into_iter().find_map(|err| match err {
        ValidationError::BrokenCycle { .. } | ValidationError::ShellGap { .. } => {
            Some(format!("{err:?}"))
        }
        ValidationError::OpenLoop { loop_, .. } => {
            if crate::direct_apex_trace::enabled() {
                trace_loop_for_acceptance_failure(brep, loop_);
            }
            Some(format!("{err:?}"))
        }
        _ => None,
    })
}

fn trace_loop_for_acceptance_failure(brep: &BRep, loop_: popkern_core::topology::LoopHandle) {
    let Some(loop_data) = brep.loops.get(loop_) else {
        return;
    };
    let face = loop_data.face;
    let start = loop_data.half_edge;
    let mut cursor = start;
    eprintln!("[direct-apex-trace] acceptance failure loop {loop_:?} face={face:?}");
    for index in 0..=brep.half_edges.len() {
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        let origin = brep.vertices[he.origin].position;
        let edge_curve = brep.edges.get(he.edge).and_then(|edge| edge.curve.as_ref());
        let curve_summary = edge_curve
            .map(|curve| {
                use popkern_core::geometry::curve::CurveEval;
                let p0 = curve.geometry.point_at(curve.t_start);
                let p1 = curve.geometry.point_at(curve.t_end);
                format!(
                    "curve=({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                    p0.x, p0.y, p0.z, p1.x, p1.y, p1.z
                )
            })
            .unwrap_or_else(|| "curve=<none>".to_string());
        eprintln!(
            "[direct-apex-trace]   [{index}] he={cursor:?} edge={:?} origin={:?}@({:.3},{:.3},{:.3}) next={:?} {curve_summary}",
            he.edge, he.origin, origin.x, origin.y, origin.z, he.next
        );
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
}

pub(super) fn native_gregory_cap_outer_g0_error(
    brep: &BRep,
    result: &crate::apex_construction::construct::ApexBuildResult,
    max_gap: f64,
) -> Option<String> {
    use popkern_core::geometry::surface::SurfaceEval;

    let cap_faces = native_gregory_cap_faces(result);
    if cap_faces.is_empty() {
        return Some("no native cap patch faces reported".to_string());
    }

    let mut worst_gap = 0.0f64;
    let mut worst_edge = None;
    let mut measured_edges = 0usize;
    for (edge_handle, edge) in brep.edges.iter() {
        let he_a = edge.half_edge;
        let Some(ha) = brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = brep.half_edges.get(ha.twin) else {
            continue;
        };
        let face_a = brep.loops[ha.loop_].face;
        let face_b = brep.loops[hb.loop_].face;
        if face_a == face_b {
            continue;
        }
        let a_is_cap = cap_faces.contains(&face_a);
        let b_is_cap = cap_faces.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        let (Some(pa), Some(pb)) = (ha.pcurve.as_ref(), hb.pcurve.as_ref()) else {
            continue;
        };
        let Some(sa) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let Some(sb) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
            continue;
        };

        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
        let uv_a = pa.point_at(0.5 * (pa_lo + pa_hi));
        let uv_b = pb.point_at(0.5 * (pb_lo + pb_hi));
        let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
        let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
        let gap = p_a.distance_to(&p_b);
        measured_edges += 1;
        if crate::direct_apex_trace::cap_g0_enabled() {
            trace_native_gregory_cap_seam_g0(
                brep,
                edge_handle,
                he_a,
                ha.twin,
                face_a,
                face_b,
                a_is_cap,
                gap,
            );
        }
        if gap > worst_gap {
            worst_gap = gap;
            worst_edge = Some(edge_handle);
        }
        if gap > max_gap {
            if crate::direct_apex_trace::cap_g0_enabled() {
                continue;
            }
            return Some(format!(
                "cap_outer_g0 {gap:.3e} on {edge_handle:?} exceeds {max_gap:.3e}"
            ));
        }
    }

    measure_exact_boundary_cap_outer_g0(
        brep,
        &cap_faces,
        max_gap,
        &mut worst_gap,
        &mut worst_edge,
        &mut measured_edges,
    );

    if measured_edges == 0 {
        return Some("no cap-vs-noncap seams measured".to_string());
    }
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace] native Gregory cap_outer_g0 postcondition: worst={worst_gap:.3e} over {measured_edges} edge(s)"
        );
    }
    if worst_gap > max_gap {
        let Some(worst_edge) = worst_edge else {
            return Some(format!(
                "cap_outer_g0 {worst_gap:.3e} exceeds {max_gap:.3e}, but no worst edge was recorded"
            ));
        };
        return Some(format!(
            "cap_outer_g0 {worst_gap:.3e} on {:?} exceeds {max_gap:.3e}",
            worst_edge
        ));
    }
    None
}

fn measure_exact_boundary_cap_outer_g0(
    brep: &BRep,
    cap_faces: &std::collections::BTreeSet<popkern_core::topology::FaceHandle>,
    max_gap: f64,
    worst_gap: &mut f64,
    worst_edge: &mut Option<popkern_core::topology::EdgeHandle>,
    measured_edges: &mut usize,
) {
    use popkern_core::geometry::surface::SurfaceEval;

    let mut cap_sides = Vec::new();
    for cap_face in cap_faces {
        let Some(cap_loop) = brep.faces.get(*cap_face).map(|face| face.outer_loop) else {
            continue;
        };
        cap_sides.extend(collect_loop_half_edges_for_acceptance(brep, cap_loop));
    }

    for cap_side in cap_sides {
        let Some(cap_side_data) = brep.half_edges.get(cap_side) else {
            continue;
        };
        let Some(cap_face) = brep
            .loops
            .get(cap_side_data.loop_)
            .map(|loop_data| loop_data.face)
        else {
            continue;
        };
        let Some(cap_surface) = brep
            .faces
            .get(cap_face)
            .and_then(|face| face.surface.as_ref())
        else {
            continue;
        };
        let Some(cap_curve) = brep
            .edges
            .get(cap_side_data.edge)
            .and_then(|edge| edge.curve.as_ref())
        else {
            continue;
        };
        for (cylinder_face, cylinder_face_data) in brep.faces.iter() {
            if cap_faces.contains(&cylinder_face) {
                continue;
            }
            let Some(cylinder_surface) = cylinder_face_data.surface.as_ref() else {
                continue;
            };
            if !matches!(cylinder_surface.geometry, SurfaceGeometry::Cylinder { .. }) {
                continue;
            }
            let cylinder_sides =
                collect_loop_half_edges_for_acceptance(brep, cylinder_face_data.outer_loop);
            if !cylinder_sides.iter().any(|cylinder_side| {
                same_oriented_exact_boundary_side(brep, cap_side, *cylinder_side)
            }) {
                continue;
            }

            let t_mid = 0.5 * (cap_curve.t_start + cap_curve.t_end);
            let seam_point = cap_curve.geometry.point_at(t_mid);
            let Some((cap_u, cap_v)) = project_point_to_uv(&cap_surface.geometry, &seam_point)
            else {
                continue;
            };
            let Some((cylinder_u, cylinder_v)) =
                project_point_to_uv(&cylinder_surface.geometry, &seam_point)
            else {
                continue;
            };
            let cap_point = cap_surface.geometry.point_at(cap_u, cap_v);
            let cylinder_point = cylinder_surface.geometry.point_at(cylinder_u, cylinder_v);
            let gap = cap_point
                .distance_to(&cylinder_point)
                .max(cap_point.distance_to(&seam_point))
                .max(cylinder_point.distance_to(&seam_point));
            *measured_edges += 1;
            if gap > *worst_gap {
                *worst_gap = gap;
                *worst_edge = Some(cap_side_data.edge);
            }
            if gap > max_gap && !crate::direct_apex_trace::cap_g0_enabled() {
                return;
            }
        }
    }
}

pub(super) fn native_gregory_cap_outer_g1_diagnostic(
    brep: &BRep,
    result: &crate::apex_construction::construct::ApexBuildResult,
    max_g0_gap: f64,
    gate: f64,
) -> Option<crate::types::ApexConstructionDiagnostic> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;

    let cap_faces = native_gregory_cap_faces(result);
    if cap_faces.is_empty() {
        return None;
    }

    let mut worst = 0.0f64;
    let mut worst_g0_gap = 0.0f64;
    let mut worst_edge = None;
    let mut worst_faces = None;
    let mut measured_edges = 0usize;
    for (edge_handle, edge) in brep.edges.iter() {
        let he_a = edge.half_edge;
        let Some(ha) = brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = brep.half_edges.get(ha.twin) else {
            continue;
        };
        let face_a = brep.loops[ha.loop_].face;
        let face_b = brep.loops[hb.loop_].face;
        if face_a == face_b {
            continue;
        }
        let a_is_cap = cap_faces.contains(&face_a);
        let b_is_cap = cap_faces.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        let (Some(pa), Some(pb)) = (ha.pcurve.as_ref(), hb.pcurve.as_ref()) else {
            continue;
        };
        let Some(sa) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let Some(sb) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
            continue;
        };

        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
        let pa_at = |t: f64| {
            let uv = pa.point_at(t);
            sa.geometry.point_at(uv.x, uv.y)
        };
        let pb_at = |t: f64| {
            let uv = pb.point_at(t);
            sb.geometry.point_at(uv.x, uv.y)
        };
        let Some(coe) = edge.curve.as_ref() else {
            continue;
        };
        let (t_lo, t_hi) = (coe.t_start, coe.t_end);
        let p3_lo = coe.geometry.point_at(t_lo);
        let p3_hi = coe.geometry.point_at(t_hi);
        let pa_forward = p3_lo.distance_to(&pa_at(pa_lo)) + p3_hi.distance_to(&pa_at(pa_hi))
            <= p3_lo.distance_to(&pa_at(pa_hi)) + p3_hi.distance_to(&pa_at(pa_lo));
        let pb_forward = p3_lo.distance_to(&pb_at(pb_lo)) + p3_hi.distance_to(&pb_at(pb_hi))
            <= p3_lo.distance_to(&pb_at(pb_hi)) + p3_hi.distance_to(&pb_at(pb_lo));

        measured_edges += 1;
        for i in 0..100 {
            let s = 0.05 + 0.90 * (i as f64) / 99.0;
            let s_a = if pa_forward { s } else { 1.0 - s };
            let s_b = if pb_forward { s } else { 1.0 - s };
            let uv_a = pa.point_at(pa_lo + s_a * (pa_hi - pa_lo));
            let uv_b = pb.point_at(pb_lo + s_b * (pb_hi - pb_lo));
            let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
            let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
            worst_g0_gap = worst_g0_gap.max(p_a.distance_to(&p_b));
            let n_a = sa
                .geometry
                .partial_u(uv_a.x, uv_a.y)
                .cross(&sa.geometry.partial_v(uv_a.x, uv_a.y));
            let n_b = sb
                .geometry
                .partial_u(uv_b.x, uv_b.y)
                .cross(&sb.geometry.partial_v(uv_b.x, uv_b.y));
            let nn_a = (n_a.x * n_a.x + n_a.y * n_a.y + n_a.z * n_a.z).sqrt();
            let nn_b = (n_b.x * n_b.x + n_b.y * n_b.y + n_b.z * n_b.z).sqrt();
            if nn_a <= 1e-300 || nn_b <= 1e-300 {
                continue;
            }
            let deficit =
                1.0 - ((n_a.x * n_b.x + n_a.y * n_b.y + n_a.z * n_b.z) / (nn_a * nn_b)).abs();
            if deficit > worst {
                worst = deficit;
                worst_edge = Some(edge_handle);
                let (cap_face, support_face) = if a_is_cap {
                    (face_a, face_b)
                } else {
                    (face_b, face_a)
                };
                worst_faces = Some((cap_face, support_face));
            }
        }
    }

    if worst <= gate {
        return None;
    }
    let (cap_face, support_face) = worst_faces?;
    let target = crate::types::CapContinuityTarget::native_cap_g1(max_g0_gap, gate);
    Some(
        crate::types::ApexConstructionDiagnostic::CapSupportContinuity(
            crate::types::CapSupportContinuityDiagnostic {
                continuity: crate::types::CapContinuityReport::declined_g1(
                    target,
                    measured_edges,
                    worst_g0_gap,
                    worst,
                ),
                cap_face,
                support_face,
                edge: worst_edge?,
                worst_g1_coplanarity: worst,
                gate,
                edges_measured: measured_edges,
            },
        ),
    )
}

#[allow(clippy::too_many_arguments)]
fn trace_native_gregory_cap_seam_g0(
    brep: &BRep,
    edge_handle: popkern_core::topology::EdgeHandle,
    he_a: popkern_core::topology::HalfEdgeHandle,
    he_b: popkern_core::topology::HalfEdgeHandle,
    face_a: popkern_core::topology::FaceHandle,
    face_b: popkern_core::topology::FaceHandle,
    a_is_cap: bool,
    midpoint_gap: f64,
) {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

    let surface_kind = |face| {
        brep.faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| match &s.geometry {
                SurfaceGeometry::Plane { .. } => "Plane",
                SurfaceGeometry::Cylinder { .. } => "Cylinder",
                SurfaceGeometry::Cone { .. } => "Cone",
                SurfaceGeometry::Sphere { .. } => "Sphere",
                SurfaceGeometry::Torus { .. } => "Torus",
                SurfaceGeometry::Nurbs(_) => "Nurbs",
                SurfaceGeometry::Gregory(_) => "Gregory",
            })
            .unwrap_or("<none>")
    };
    let he_a_data = &brep.half_edges[he_a];
    let he_b_data = &brep.half_edges[he_b];
    eprintln!(
        "[direct-apex-trace cap-g0] seam edge={edge_handle:?} he=({he_a:?},{he_b:?}) \
         faces=({face_a:?}:{}, {face_b:?}:{}) cap_side={} midpoint_gap={midpoint_gap:.6e}",
        surface_kind(face_a),
        surface_kind(face_b),
        if a_is_cap { "A" } else { "B" }
    );

    let Some(edge_curve) = brep
        .edges
        .get(edge_handle)
        .and_then(|edge| edge.curve.as_ref())
    else {
        eprintln!("[direct-apex-trace cap-g0]   edge has no 3D curve");
        return;
    };
    let Some(surface_a) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
        return;
    };
    let Some(surface_b) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
        return;
    };
    let (Some(pcurve_a), Some(pcurve_b)) = (he_a_data.pcurve.as_ref(), he_b_data.pcurve.as_ref())
    else {
        eprintln!("[direct-apex-trace cap-g0]   missing pcurve on one side");
        return;
    };

    let edge_point = |s: f64| {
        let t = edge_curve.t_start + s * (edge_curve.t_end - edge_curve.t_start);
        edge_curve.geometry.point_at(t)
    };
    let side_point = |pcurve: &popkern_core::geometry::pcurve::Pcurve,
                      surface: &popkern_core::geometry::surface::SurfaceOnFace,
                      s: f64| {
        let (lo, hi) = pcurve.parameter_range();
        let uv = pcurve.point_at(lo + s * (hi - lo));
        surface.geometry.point_at(uv.x, uv.y)
    };
    let forward = |pcurve: &popkern_core::geometry::pcurve::Pcurve,
                   surface: &popkern_core::geometry::surface::SurfaceOnFace| {
        let p0 = side_point(pcurve, surface, 0.0);
        let p1 = side_point(pcurve, surface, 1.0);
        p0.distance_to(&edge_point(0.0)) + p1.distance_to(&edge_point(1.0))
            <= p0.distance_to(&edge_point(1.0)) + p1.distance_to(&edge_point(0.0))
    };
    let a_forward = forward(pcurve_a, surface_a);
    let b_forward = forward(pcurve_b, surface_b);

    for sample in [0.0, 0.5, 1.0] {
        let s_a = if a_forward { sample } else { 1.0 - sample };
        let s_b = if b_forward { sample } else { 1.0 - sample };
        let p_edge = edge_point(sample);
        let p_a = side_point(pcurve_a, surface_a, s_a);
        let p_b = side_point(pcurve_b, surface_b, s_b);
        eprintln!(
            "[direct-apex-trace cap-g0]   s={sample:.1} cap_noncap={:.6e} \
             sideA_edge={:.6e} sideB_edge={:.6e}",
            p_a.distance_to(&p_b),
            p_a.distance_to(&p_edge),
            p_b.distance_to(&p_edge)
        );
    }
}

pub(super) fn native_gregory_has_expected_cylinder_seams(
    brep: &BRep,
    result: &crate::apex_construction::construct::ApexBuildResult,
    n_fillets: usize,
) -> bool {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let cap_faces = native_gregory_cap_faces(result);
    if cap_faces.is_empty() {
        return false;
    }

    let mut cylinder_seams = 0usize;
    let mut cylinder_faces = std::collections::BTreeSet::new();
    for (_, edge) in brep.edges.iter() {
        let he_a = edge.half_edge;
        let Some(ha) = brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = brep.half_edges.get(ha.twin) else {
            continue;
        };
        let face_a = brep.loops[ha.loop_].face;
        let face_b = brep.loops[hb.loop_].face;
        if face_a == face_b {
            continue;
        }
        let a_is_cap = cap_faces.contains(&face_a);
        let b_is_cap = cap_faces.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        let neighbor = if a_is_cap { face_b } else { face_a };
        if matches!(
            brep.faces
                .get(neighbor)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry),
            Some(SurfaceGeometry::Cylinder { .. })
        ) {
            cylinder_seams += 1;
            cylinder_faces.insert(neighbor);
        }
    }

    if cylinder_faces.len() < n_fillets {
        count_exact_boundary_cylinder_seams(
            brep,
            &cap_faces,
            &mut cylinder_faces,
            &mut cylinder_seams,
        );
    }

    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace] native Gregory seam postcondition: {cylinder_seams} seams across {} distinct cylinder face(s), expected {n_fillets}",
            cylinder_faces.len()
        );
    }
    cylinder_faces.len() >= n_fillets
}

fn native_gregory_cap_faces(
    result: &crate::apex_construction::construct::ApexBuildResult,
) -> std::collections::BTreeSet<popkern_core::topology::FaceHandle> {
    result
        .corner_patches
        .iter()
        .map(|patch| patch.face)
        .chain(result.exact_boundary_report.committed_faces.iter().copied())
        .collect()
}

fn count_exact_boundary_cylinder_seams(
    brep: &BRep,
    cap_faces: &std::collections::BTreeSet<popkern_core::topology::FaceHandle>,
    cylinder_faces: &mut std::collections::BTreeSet<popkern_core::topology::FaceHandle>,
    cylinder_seams: &mut usize,
) {
    let mut cap_sides = Vec::new();
    for cap_face in cap_faces {
        let Some(cap_loop) = brep.faces.get(*cap_face).map(|face| face.outer_loop) else {
            continue;
        };
        cap_sides.extend(collect_loop_half_edges_for_acceptance(brep, cap_loop));
    }

    for (face, face_data) in brep.faces.iter() {
        if cap_faces.contains(&face) || cylinder_faces.contains(&face) {
            continue;
        }
        let Some(surface) = face_data.surface.as_ref() else {
            continue;
        };
        if !matches!(surface.geometry, SurfaceGeometry::Cylinder { .. }) {
            continue;
        }
        let cylinder_sides = collect_loop_half_edges_for_acceptance(brep, face_data.outer_loop);
        if cap_sides.iter().any(|cap_side| {
            cylinder_sides.iter().any(|cylinder_side| {
                same_oriented_exact_boundary_side(brep, *cap_side, *cylinder_side)
            })
        }) {
            *cylinder_seams += 1;
            cylinder_faces.insert(face);
        }
    }
}

fn collect_loop_half_edges_for_acceptance(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Vec<popkern_core::topology::HalfEdgeHandle> {
    let Some(loop_data) = brep.loops.get(loop_handle) else {
        return Vec::new();
    };
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return Vec::new();
    }
    let mut cycle = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        cycle.push(cursor);
        let Some(next) = brep.half_edges.get(cursor).map(|half_edge| half_edge.next) else {
            return Vec::new();
        };
        cursor = next;
        if cursor == start {
            return cycle;
        }
    }
    Vec::new()
}

fn same_oriented_exact_boundary_side(
    brep: &BRep,
    cap_side: popkern_core::topology::HalfEdgeHandle,
    cylinder_side: popkern_core::topology::HalfEdgeHandle,
) -> bool {
    let Some(cap_start) = brep.half_edges.get(cap_side).map(|data| data.origin) else {
        return false;
    };
    let Some(cylinder_start) = brep.half_edges.get(cylinder_side).map(|data| data.origin) else {
        return false;
    };
    if cap_start != cylinder_start {
        return false;
    }
    let Some(cap_end) = half_edge_destination_vertex_for_acceptance(brep, cap_side) else {
        return false;
    };
    let Some(cylinder_end) = half_edge_destination_vertex_for_acceptance(brep, cylinder_side)
    else {
        return false;
    };
    cap_end == cylinder_end
}

fn half_edge_destination_vertex_for_acceptance(
    brep: &BRep,
    half_edge: popkern_core::topology::HalfEdgeHandle,
) -> Option<popkern_core::topology::VertexHandle> {
    let data = brep.half_edges.get(half_edge)?;
    if let Some(twin) = brep.half_edges.get(data.twin) {
        return Some(twin.origin);
    }
    if !brep.boundary_edges.contains(&data.edge) {
        return brep
            .half_edges
            .get(data.next)
            .map(|next_half_edge| next_half_edge.origin);
    }
    let edge = brep.edges.get(data.edge)?;
    let Some(curve) = edge.curve.as_ref() else {
        return brep
            .half_edges
            .get(data.next)
            .map(|next_half_edge| next_half_edge.origin);
    };
    let start_position = brep.vertices.get(data.origin)?.position;
    let p0 = curve.geometry.point_at(curve.t_start);
    let p1 = curve.geometry.point_at(curve.t_end);
    let destination = if start_position.distance_to(&p0) <= start_position.distance_to(&p1) {
        p1
    } else {
        p0
    };
    let tolerance = brep.tolerance.linear().max(1.0e-9);
    brep.vertices
        .iter()
        .filter_map(|(vertex, data)| {
            let distance = data.position.distance_to(&destination);
            (distance <= tolerance).then_some((vertex, distance))
        })
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(vertex, _)| vertex)
}
