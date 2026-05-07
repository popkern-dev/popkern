//! Stage C and C.5 cap-loop construction for constructive Gregory apexes.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::fillet_surface::FilletSurfaceResult;

use super::constructive_geometry::{
    assign_fillet_strip_surface, attach_fillet_end_curve_to_edge, planned_fillet_index_for_sector,
};
use super::constructive_topology::{
    find_outgoing_half_edge, find_two_origins_on_all_faces, find_two_origins_on_face,
    loop_segment_contains_origin,
};
use super::error::ApexBuildError;
use super::plan::ApexPlan;
use super::repeated_apex_loop_plan::plan_repeated_apex_loop_normalization;
use super::repeated_apex_partition::RepeatedApexPartitionExecutionOptions;
use super::repeated_apex_serial::try_normalize_repeated_apex_loop;

pub(super) fn stage_c_kef_cube_edges_interleaved(
    brep: &mut BRep,
    plan: &ApexPlan,
    cube_edges_pre_a: &[EdgeHandle],
    apex_vertex_handles: &[VertexHandle],
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<(FaceHandle, Vec<Option<EdgeHandle>>), ApexBuildError> {
    let trace = std::env::var_os("POPKERN_A2_TRACE").is_some();
    let n = plan.cap_kind.n_fillets();
    let mut fillet_ends = vec![None; n];

    for &edge in cube_edges_pre_a {
        if !brep.edges.contains_key(edge) {
            continue;
        }
        let he = brep.edges[edge].half_edge;
        let twin = brep.half_edges[he].twin;
        let he_loop = brep.half_edges[he].loop_;
        let twin_loop = brep.half_edges[twin].loop_;
        if he_loop == twin_loop {
            if trace {
                eprintln!(
                    "[a2-trace §14.8/C] cube edge {edge:?} is same-face bridge, deferring to Stage D"
                );
            }
            continue;
        }
        if trace {
            let he_face = brep.loops[he_loop].face;
            let twin_face = brep.loops[twin_loop].face;
            eprintln!(
                "[a2-trace §14.8/C] kef cube edge {edge:?} (he in {he_face:?}, twin in {twin_face:?})"
            );
        }
        brep.kef(he).map_err(|e| ApexBuildError::EulerOpFailed {
            op: "kef",
            step: 200,
            cause: format!("stage_c_kef_cube_edges_interleaved: kef edge {edge:?}: {e}"),
        })?;

        for (k, fillet_end) in fillet_ends.iter_mut().enumerate().take(n) {
            if fillet_end.is_some() {
                continue;
            }
            if let Some(created) = try_stage_c5_mef_single_fillet_end(
                brep,
                plan,
                k,
                apex_vertex_handles,
                fillet_surfaces,
            )? {
                if trace {
                    eprintln!(
                        "[a2-trace §14.8/C.5] interleaved FilletEnd[{k}] after cube edge {edge:?}: {:?}",
                        created
                    );
                }
                *fillet_end = Some(created);
            }
        }
    }

    let cap_face = find_cap_face_from_apex_outgoing(brep, plan.apex_vertex).ok_or_else(|| {
        ApexBuildError::PlanMismatch {
            reason: "stage_c_kef_cube_edges_interleaved: V has no outgoing half-edges".to_string(),
        }
    })?;
    trace_stage_c_result(brep, cap_face, plan.apex_vertex);
    brep.vertices[plan.apex_vertex].outgoing = find_outgoing_half_edge(brep, plan.apex_vertex)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: "stage_c_kef_cube_edges_interleaved: V outgoing repair failed".to_string(),
        })?;

    Ok((cap_face, fillet_ends))
}

fn find_cap_face_from_apex_outgoing(brep: &BRep, v: VertexHandle) -> Option<FaceHandle> {
    let cap_he = brep
        .half_edges
        .iter()
        .find(|(_, h)| h.origin == v)
        .map(|(h, _)| h)?;
    let cap_loop = brep.half_edges[cap_he].loop_;
    Some(brep.loops[cap_loop].face)
}

fn trace_stage_c_result(brep: &BRep, cap_face: FaceHandle, v: VertexHandle) {
    let trace = std::env::var_os("POPKERN_A2_TRACE").is_some();
    if !trace {
        return;
    }
    let outgoings: Vec<popkern_core::topology::HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, h)| h.origin == v)
        .map(|(h, _)| h)
        .collect();
    eprintln!(
        "[a2-trace §14.8/C] selective C complete: cap_face={cap_face:?}, V has {} outgoings",
        outgoings.len()
    );
    for he in &outgoings {
        let l = brep.half_edges[*he].loop_;
        let f = brep.loops[l].face;
        let dest = brep.half_edges[brep.half_edges[*he].twin].origin;
        eprintln!("[a2-trace §14.8/C]   he={he:?} -> {dest:?} face={f:?}");
    }
    let cap_loop_handle = brep.faces[cap_face].outer_loop;
    let start = brep.loops[cap_loop_handle].half_edge;
    let mut cursor = start;
    let mut chain: Vec<String> = Vec::new();
    for _ in 0..(brep.half_edges.len() + 1) {
        let he = &brep.half_edges[cursor];
        chain.push(format!("{:?}", he.origin));
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    eprintln!(
        "[a2-trace §14.8/C]   cap_face={cap_face:?} loop ({} half-edges): {}",
        chain.len(),
        chain.join(" -> ")
    );
}

pub(super) fn stage_c5_mef_remaining_fillet_ends(
    brep: &mut BRep,
    plan: &ApexPlan,
    cap_face: &mut FaceHandle,
    apex_vertex_handles: &[VertexHandle],
    face_edge_handles: &[EdgeHandle],
    fillet_surfaces: &[FilletSurfaceResult],
    fillet_ends: &mut [Option<EdgeHandle>],
) -> Result<bool, ApexBuildError> {
    let n = plan.cap_kind.n_fillets();
    let trace = std::env::var_os("POPKERN_A2_TRACE").is_some();

    if trace {
        trace_cap_loop_apex_arcs(brep, *cap_face, plan.apex_vertex, apex_vertex_handles, n);
    }

    if fillet_ends.iter().all(Option::is_none)
        && let Some((normalized_cap_face, created)) = try_normalize_repeated_apex_loop(
            brep,
            plan,
            *cap_face,
            apex_vertex_handles,
            face_edge_handles,
            fillet_surfaces,
            repeated_apex_partition_execution_options_from_environment(),
        )?
    {
        *cap_face = normalized_cap_face;
        for (slot, edge) in fillet_ends.iter_mut().zip(created) {
            *slot = Some(edge);
        }
        return Ok(true);
    }

    for (k, fillet_end) in fillet_ends.iter_mut().enumerate().take(n) {
        if fillet_end.is_some() {
            continue;
        }
        let created = mef_fillet_end_on_face(
            brep,
            plan,
            *cap_face,
            k,
            apex_vertex_handles,
            fillet_surfaces,
        )?;
        *fillet_end = Some(created);
    }

    Ok(false)
}

fn repeated_apex_partition_execution_options_from_environment()
-> RepeatedApexPartitionExecutionOptions {
    RepeatedApexPartitionExecutionOptions::new(
        std::env::var_os("POPKERN_A2_APPLY_CYLINDER_TRIM_REWRITE").is_some(),
        std::env::var_os("POPKERN_A2_REQUIRE_CYLINDER_TRIM_REWRITE_READY").is_some(),
        std::env::var_os("POPKERN_A2_REQUIRE_CYLINDER_TRIM_REBUILD_PLAN").is_some(),
        std::env::var_os("POPKERN_A2_ADOPT_CYLINDER_FILLET_END_ENDPOINTS").is_some(),
        std::env::var_os("POPKERN_A2_REQUIRE_SHARED_CYLINDER_FILLET_END").is_some(),
    )
}

fn try_stage_c5_mef_single_fillet_end(
    brep: &mut BRep,
    plan: &ApexPlan,
    k: usize,
    apex_vertex_handles: &[VertexHandle],
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<Option<EdgeHandle>, ApexBuildError> {
    let trace = std::env::var_os("POPKERN_A2_TRACE").is_some();
    let t_a_k = apex_vertex_handles[2 * k];
    let t_b_k = apex_vertex_handles[2 * k + 1];
    for (face, he_origin, he_dest) in find_two_origins_on_all_faces(brep, t_a_k, t_b_k) {
        let forward_contains_v =
            loop_segment_contains_origin(brep, he_origin, he_dest, plan.apex_vertex);
        let reverse_contains_v =
            loop_segment_contains_origin(brep, he_dest, he_origin, plan.apex_vertex);
        if trace {
            eprintln!(
                "[a2-trace §14.8/C.5] interleaved candidate FilletEnd[{k}] on {face:?}: forward_contains_v={forward_contains_v}, reverse_contains_v={reverse_contains_v}"
            );
        }
        if forward_contains_v && reverse_contains_v {
            continue;
        }
        return mef_fillet_end_on_face(brep, plan, face, k, apex_vertex_handles, fillet_surfaces)
            .map(Some);
    }
    Ok(None)
}

fn mef_fillet_end_on_face(
    brep: &mut BRep,
    plan: &ApexPlan,
    face: FaceHandle,
    k: usize,
    apex_vertex_handles: &[VertexHandle],
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<EdgeHandle, ApexBuildError> {
    let trace = std::env::var_os("POPKERN_A2_TRACE").is_some();
    let t_a_k = apex_vertex_handles[2 * k];
    let t_b_k = apex_vertex_handles[2 * k + 1];

    let (he_origin, he_dest) =
        find_two_origins_on_face(brep, face, t_a_k, t_b_k).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "stage_c5_mef_fillet_ends: T_a[{k}]={t_a_k:?} and/or T_b[{k}]={t_b_k:?} \
                     not on the same loop of face={face:?}"
                ),
            }
        })?;

    let forward_contains_v =
        loop_segment_contains_origin(brep, he_origin, he_dest, plan.apex_vertex);
    let reverse_contains_v =
        loop_segment_contains_origin(brep, he_dest, he_origin, plan.apex_vertex);
    let (strip_start, strip_end) = match (forward_contains_v, reverse_contains_v) {
        (true, false) => (he_dest, he_origin),
        (false, _) => (he_origin, he_dest),
        (true, true) if std::env::var_os("POPKERN_A2_ALLOW_REPEATED_APEX_LOOP").is_some() => {
            (he_origin, he_dest)
        }
        (true, true) => {
            let diagnosis = plan_repeated_apex_loop_normalization(
                brep,
                face,
                plan.apex_vertex,
                apex_vertex_handles,
                plan.cap_kind.n_fillets(),
            )
            .map(|normalization_plan| {
                normalization_plan.diagnostic_summary(plan.cap_kind.n_fillets())
            })
            .unwrap_or_else(|| "no repeated-apex loop normalization plan available".to_string());
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "stage_c5_mef_fillet_ends: FilletEnd[{k}] has no V-free cap-face arc \
                     between T_a[{k}]={t_a_k:?} and T_b[{k}]={t_b_k:?}; Stage C produced \
                     a repeated-apex cap loop that cannot isolate a commercial-grade \
                     cylinder strip before apex removal; {diagnosis}"
                ),
            });
        }
    };

    let mef_result =
        brep.mef(strip_start, strip_end)
            .map_err(|e| ApexBuildError::EulerOpFailed {
                op: "mef",
                step: 200 + k,
                cause: format!("FilletEnd[{k}]: {e}"),
            })?;

    let fillet_index = planned_fillet_index_for_sector(plan, k)?;
    let curve = &plan.target_topology.apex_edges[2 * k].curve;
    attach_fillet_end_curve_to_edge(brep, mef_result.edge, curve, fillet_index, fillet_surfaces)?;
    assign_fillet_strip_surface(brep, mef_result.face, fillet_index, fillet_surfaces)?;

    if trace {
        eprintln!(
            "[a2-trace §14.8/C.5] mef FilletEnd[{k}] T_a[{k}]={t_a_k:?} -> T_b[{k}]={t_b_k:?} \
             on face={face:?} -> edge={:?}, cyl-strip={:?}, forward_contains_v={}, reverse_contains_v={}",
            mef_result.edge, mef_result.face, forward_contains_v, reverse_contains_v
        );
    }
    Ok(mef_result.edge)
}

pub(super) fn collect_required_fillet_ends(
    fillet_ends: Vec<Option<EdgeHandle>>,
) -> Result<Vec<EdgeHandle>, ApexBuildError> {
    fillet_ends
        .into_iter()
        .enumerate()
        .map(|(k, edge)| {
            edge.ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!("FilletEnd[{k}] was not constructed"),
            })
        })
        .collect()
}

fn trace_cap_loop_apex_arcs(
    brep: &BRep,
    cap_face: FaceHandle,
    apex_vertex: VertexHandle,
    apex_vertex_handles: &[VertexHandle],
    n_fillets: usize,
) {
    let Some(face_data) = brep.faces.get(cap_face) else {
        return;
    };
    let mut loops = Vec::with_capacity(1 + face_data.inner_loops.len());
    loops.push(("outer".to_string(), face_data.outer_loop));
    for (index, loop_handle) in face_data.inner_loops.iter().copied().enumerate() {
        loops.push((format!("inner[{index}]"), loop_handle));
    }

    for (loop_label, loop_handle) in loops {
        let Some(loop_data) = brep.loops.get(loop_handle) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut vertices = Vec::new();
        let mut guard = 0usize;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(half_edge) = brep.half_edges.get(cursor) {
            vertices.push(half_edge.origin);
            cursor = half_edge.next;
            guard += 1;
            if cursor == start || guard >= max_iters {
                break;
            }
        }

        let apex_positions: Vec<usize> = vertices
            .iter()
            .enumerate()
            .filter_map(|(index, vertex)| (*vertex == apex_vertex).then_some(index))
            .collect();
        eprintln!(
            "[a2-trace §14.8/C.5] cap loop {loop_label}={loop_handle:?}: {} vertices, {} apex occurrences",
            vertices.len(),
            apex_positions.len()
        );

        for k in 0..n_fillets {
            let t_a = apex_vertex_handles[2 * k];
            let t_b = apex_vertex_handles[2 * k + 1];
            let t_a_positions = vertex_positions(&vertices, t_a);
            let t_b_positions = vertex_positions(&vertices, t_b);
            eprintln!(
                "[a2-trace §14.8/C.5]   T_a[{k}]={t_a:?}@{:?}, T_b[{k}]={t_b:?}@{:?}",
                t_a_positions, t_b_positions
            );
        }

        if apex_positions.len() < 2 || vertices.is_empty() {
            continue;
        }

        for (arc_index, &start_pos) in apex_positions.iter().enumerate() {
            let end_pos = apex_positions[(arc_index + 1) % apex_positions.len()];
            let mut pos = (start_pos + 1) % vertices.len();
            let mut labels = Vec::new();
            while pos != end_pos {
                labels.push(trace_vertex_label(
                    vertices[pos],
                    apex_vertex,
                    apex_vertex_handles,
                    n_fillets,
                ));
                pos = (pos + 1) % vertices.len();
            }
            let contents = if labels.is_empty() {
                "<empty>".to_string()
            } else {
                labels.join(" -> ")
            };
            eprintln!(
                "[a2-trace §14.8/C.5]   apex-arc[{arc_index}] V@{start_pos} -> V@{end_pos}: {contents}"
            );
        }
    }
}

fn vertex_positions(vertices: &[VertexHandle], needle: VertexHandle) -> Vec<usize> {
    vertices
        .iter()
        .enumerate()
        .filter_map(|(index, vertex)| (*vertex == needle).then_some(index))
        .collect()
}

fn trace_vertex_label(
    vertex: VertexHandle,
    apex_vertex: VertexHandle,
    apex_vertex_handles: &[VertexHandle],
    n_fillets: usize,
) -> String {
    if vertex == apex_vertex {
        return format!("V({vertex:?})");
    }
    for k in 0..n_fillets {
        if apex_vertex_handles[2 * k] == vertex {
            return format!("T_a[{k}]({vertex:?})");
        }
        if apex_vertex_handles[2 * k + 1] == vertex {
            return format!("T_b[{k}]({vertex:?})");
        }
    }
    format!("{vertex:?}")
}
