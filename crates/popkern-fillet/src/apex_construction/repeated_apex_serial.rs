//! Serial repeated-apex normalization fallback for constructive Gregory apexes.
//!
//! The planned partition is the preferred repeated-apex path. This module keeps
//! the older serial `mef` scheduler isolated so its preflight, rollback, and
//! diagnostic behavior can be reviewed independently.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::fillet_surface::FilletSurfaceResult;

use super::constructive_geometry::{
    assign_fillet_strip_surface, attach_fillet_end_curve_to_edge, planned_fillet_index_for_sector,
};
use super::constructive_topology::{
    assert_normalized_cap_boundary, collect_outer_loop_vertices, cyclic_vertices_match,
    find_two_origins_on_all_faces, find_two_origins_on_face,
    validate_repeated_apex_normalization_structure,
};
use super::error::ApexBuildError;
use super::plan::ApexPlan;
use super::repeated_apex_loop_plan::plan_repeated_apex_loop_normalization;
use super::repeated_apex_partition::try_build_planned_repeated_apex_partition;

pub(super) fn try_normalize_repeated_apex_loop(
    brep: &mut BRep,
    plan: &ApexPlan,
    cap_face: FaceHandle,
    apex_vertex_handles: &[VertexHandle],
    face_edge_handles: &[EdgeHandle],
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<Option<(FaceHandle, Vec<EdgeHandle>)>, ApexBuildError> {
    let n = plan.cap_kind.n_fillets();
    let Some(normalization_plan) = plan_repeated_apex_loop_normalization(
        brep,
        cap_face,
        plan.apex_vertex,
        apex_vertex_handles,
        n,
    ) else {
        return Ok(None);
    };

    if !normalization_plan.is_commercial_grade_candidate(n) {
        return Ok(None);
    }

    let snapshot = brep.clone();
    match normalize_repeated_apex_loop_in_place(
        brep,
        plan,
        cap_face,
        apex_vertex_handles,
        face_edge_handles,
        fillet_surfaces,
        n,
    ) {
        Ok(created) => Ok(Some(created)),
        Err(err) => {
            *brep = snapshot;
            Err(err)
        }
    }
}

fn normalize_repeated_apex_loop_in_place(
    brep: &mut BRep,
    plan: &ApexPlan,
    cap_face: FaceHandle,
    apex_vertex_handles: &[VertexHandle],
    face_edge_handles: &[EdgeHandle],
    fillet_surfaces: &[FilletSurfaceResult],
    n_fillets: usize,
) -> Result<(FaceHandle, Vec<EdgeHandle>), ApexBuildError> {
    let trace = std::env::var_os("POPKERN_A2_TRACE").is_some();
    if trace {
        trace_repeated_apex_partition_inputs(
            brep,
            cap_face,
            apex_vertex_handles,
            face_edge_handles,
            n_fillets,
        );
    }
    if let Some(partition) = try_build_planned_repeated_apex_partition(
        brep,
        plan,
        cap_face,
        apex_vertex_handles,
        face_edge_handles,
        fillet_surfaces,
        n_fillets,
    )? {
        return Ok(partition);
    }

    let mut created = vec![None; n_fillets];
    let mut current_cap_face = cap_face;
    let mut created_sectors = vec![false; n_fillets];

    for step in 0..n_fillets {
        let (k, construction_face, he_origin, he_dest) = plan_next_repeated_apex_sector(
            brep,
            current_cap_face,
            apex_vertex_handles,
            &created_sectors,
            n_fillets,
        )?;
        let t_a_k = apex_vertex_handles[2 * k];
        let t_b_k = apex_vertex_handles[2 * k + 1];
        if !half_edges_have_origins(brep, he_origin, he_dest, t_a_k, t_b_k) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "repeated-apex normalization: scheduled FilletEnd[{k}] half-edges no longer \
                     originate at T_a[{k}]={t_a_k:?} and T_b[{k}]={t_b_k:?}"
                ),
            });
        }

        let mef_result =
            brep.mef(he_origin, he_dest)
                .map_err(|e| ApexBuildError::EulerOpFailed {
                    op: "mef",
                    step: 300 + step,
                    cause: format!("repeated-apex FilletEnd[{k}]: {e}"),
                })?;

        let fillet_index = planned_fillet_index_for_sector(plan, k)?;
        let curve = &plan.target_topology.apex_edges[2 * k].curve;
        attach_fillet_end_curve_to_edge(
            brep,
            mef_result.edge,
            curve,
            fillet_index,
            fillet_surfaces,
        )?;
        assign_fillet_strip_surface(brep, mef_result.face, fillet_index, fillet_surfaces)?;
        created[k] = Some(mef_result.edge);
        created_sectors[k] = true;

        if let Some(normalized_cap_face) =
            find_normalized_cap_face(brep, apex_vertex_handles, n_fillets)
        {
            current_cap_face = normalized_cap_face;
        } else if let Some(residual_cap_face) = find_residual_face_with_constructible_sectors(
            brep,
            apex_vertex_handles,
            &created_sectors,
            n_fillets,
            mef_result.face,
        ) {
            current_cap_face = residual_cap_face;
        }

        if trace {
            eprintln!(
                "[a2-trace §14.9/normalize] step[{step}] mef FilletEnd[{k}] on \
                 repeated-apex loop: face={construction_face:?}, edge={:?}, strip_face={:?}, \
                 residual_cap_face={:?}",
                mef_result.edge, mef_result.face, current_cap_face
            );
        }
    }

    let current_cap_face = find_normalized_cap_face(brep, apex_vertex_handles, n_fillets)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason:
                "repeated-apex normalization did not produce a face with the expected cap boundary"
                    .to_string(),
        })?;
    assert_normalized_cap_boundary(brep, current_cap_face, apex_vertex_handles, n_fillets)?;
    validate_repeated_apex_normalization_structure(brep)?;
    let created = created
        .into_iter()
        .enumerate()
        .map(|(k, edge)| {
            edge.ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!("repeated-apex normalization did not create FilletEnd[{k}]"),
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok((current_cap_face, created))
}

fn plan_next_repeated_apex_sector(
    brep: &BRep,
    preferred_face: FaceHandle,
    apex_vertex_handles: &[VertexHandle],
    created_sectors: &[bool],
    n_fillets: usize,
) -> Result<(usize, FaceHandle, HalfEdgeHandle, HalfEdgeHandle), ApexBuildError> {
    let mut blocked = Vec::new();
    for k in 0..n_fillets {
        if created_sectors[k] {
            continue;
        }
        let t_a_k = apex_vertex_handles[2 * k];
        let t_b_k = apex_vertex_handles[2 * k + 1];
        let mut candidates = Vec::new();
        if let Some((he_origin, he_dest)) =
            find_two_origins_on_face(brep, preferred_face, t_a_k, t_b_k)
        {
            candidates.push((preferred_face, he_origin, he_dest));
            candidates.push((preferred_face, he_dest, he_origin));
        }
        for (face, he_origin, he_dest) in find_two_origins_on_all_faces(brep, t_a_k, t_b_k)
            .into_iter()
            .filter(|(face, _, _)| *face != preferred_face)
        {
            candidates.push((face, he_origin, he_dest));
            candidates.push((face, he_dest, he_origin));
        }
        if candidates.is_empty() {
            blocked.push(format!("FilletEnd[{k}] endpoints are not both on any face"));
            continue;
        }

        for (face, he_origin, he_dest) in candidates {
            let mut candidate = brep.clone();
            if candidate.mef(he_origin, he_dest).is_err() {
                blocked.push(format!("FilletEnd[{k}] candidate mef failed on {face:?}"));
                continue;
            }
            let mut candidate_created = created_sectors.to_vec();
            candidate_created[k] = true;
            if remaining_sectors_are_constructible(
                &candidate,
                apex_vertex_handles,
                &candidate_created,
                n_fillets,
            ) {
                return Ok((k, face, he_origin, he_dest));
            }
            blocked.push(format!(
                "FilletEnd[{k}] on {face:?} strands at least one remaining sector"
            ));
        }
    }

    Err(ApexBuildError::PlanMismatch {
        reason: format!(
            "repeated-apex normalization: no next FilletEnd preserves constructibility \
             from preferred face {preferred_face:?}; {}",
            blocked.join("; ")
        ),
    })
}

fn remaining_sectors_are_constructible(
    brep: &BRep,
    apex_vertex_handles: &[VertexHandle],
    created_sectors: &[bool],
    n_fillets: usize,
) -> bool {
    (0..n_fillets).all(|k| {
        if created_sectors[k] {
            return true;
        }
        let t_a_k = apex_vertex_handles[2 * k];
        let t_b_k = apex_vertex_handles[2 * k + 1];
        !find_two_origins_on_all_faces(brep, t_a_k, t_b_k).is_empty()
    })
}

fn half_edges_have_origins(
    brep: &BRep,
    he_a: HalfEdgeHandle,
    he_b: HalfEdgeHandle,
    vertex_a: VertexHandle,
    vertex_b: VertexHandle,
) -> bool {
    let Some(origin_a) = brep.half_edges.get(he_a).map(|half_edge| half_edge.origin) else {
        return false;
    };
    let Some(origin_b) = brep.half_edges.get(he_b).map(|half_edge| half_edge.origin) else {
        return false;
    };
    (origin_a == vertex_a && origin_b == vertex_b) || (origin_a == vertex_b && origin_b == vertex_a)
}

fn find_residual_face_with_constructible_sectors(
    brep: &BRep,
    apex_vertex_handles: &[VertexHandle],
    created_sectors: &[bool],
    n_fillets: usize,
    excluded_face: FaceHandle,
) -> Option<FaceHandle> {
    for (face, _) in &brep.faces {
        if face == excluded_face {
            continue;
        }
        let supports_remaining_sector = (0..n_fillets).any(|k| {
            !created_sectors[k]
                && find_two_origins_on_face(
                    brep,
                    face,
                    apex_vertex_handles[2 * k],
                    apex_vertex_handles[2 * k + 1],
                )
                .is_some()
        });
        if supports_remaining_sector {
            return Some(face);
        }
    }
    None
}

fn find_normalized_cap_face(
    brep: &BRep,
    apex_vertex_handles: &[VertexHandle],
    n_fillets: usize,
) -> Option<FaceHandle> {
    let mut expected = Vec::with_capacity(2 * n_fillets);
    for k in 0..n_fillets {
        expected.push(apex_vertex_handles[2 * k]);
        expected.push(apex_vertex_handles[2 * k + 1]);
    }
    for (face, _) in &brep.faces {
        let Some(vertices) = collect_outer_loop_vertices(brep, face) else {
            continue;
        };
        if vertices.len() == expected.len() && cyclic_vertices_match(&vertices, &expected) {
            return Some(face);
        }
    }
    None
}

fn trace_repeated_apex_partition_inputs(
    brep: &BRep,
    cap_face: FaceHandle,
    apex_vertex_handles: &[VertexHandle],
    face_edge_handles: &[EdgeHandle],
    n_fillets: usize,
) {
    eprintln!("[a2-trace §14.9/partition] cap_face={cap_face:?}");
    for k in 0..n_fillets {
        let edge = face_edge_handles.get(k).copied();
        let Some(edge) = edge else {
            eprintln!("[a2-trace §14.9/partition]   FaceEdge[{k}] missing");
            continue;
        };
        let Some(edge_data) = brep.edges.get(edge) else {
            eprintln!("[a2-trace §14.9/partition]   FaceEdge[{k}] {edge:?} dead");
            continue;
        };
        let he = edge_data.half_edge;
        let Some(he_data) = brep.half_edges.get(he) else {
            continue;
        };
        let twin = he_data.twin;
        let he_face = brep
            .loops
            .get(he_data.loop_)
            .map(|loop_data| loop_data.face);
        let twin_face = brep.half_edges.get(twin).and_then(|twin_data| {
            brep.loops
                .get(twin_data.loop_)
                .map(|loop_data| loop_data.face)
        });
        let he_origin = he_data.origin;
        let twin_origin = brep.half_edges.get(twin).map(|twin_data| twin_data.origin);
        eprintln!(
            "[a2-trace §14.9/partition]   FaceEdge[{k}] edge={edge:?} he={he:?} \
             {he_origin:?}->{twin_origin:?} face={he_face:?}; twin={twin:?} face={twin_face:?}; \
             expected {:?}->{:?}",
            apex_vertex_handles[2 * k + 1],
            apex_vertex_handles[2 * ((k + 1) % n_fillets)]
        );
    }
}
