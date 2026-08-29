//! Apex-vertex removal for constructive apex building.

use popkern_core::topology::{BRep, FaceHandle};

use super::constructive_topology::{
    collect_outer_loop_vertices, production_faces_containing_vertex,
};
use super::error::ApexBuildError;
use super::plan::ApexPlan;
use super::repeated_apex_partition_topology::repair_stale_shell_anchors;

/// Remove V via residual same-face bridges. For Sphere, V has a single
/// same-face bridge after V-incident edge removal; `kev_v_bridge` removes it
/// plus V in one op. For Gregory, V has multiple same-face MEV spurs to each
/// T-vertex plus the residual preexisting cross-face edge bridge from cap-loop
/// construction; we dissolve the spurs first, then `kev_v_bridge` consumes
/// the last bridge.
pub(super) fn remove_apex_vertex(
    brep: &mut BRep,
    plan: &ApexPlan,
    cap_face: FaceHandle,
    repeated_apex_normalized: bool,
    require_v_free_repeated_apex: bool,
) -> Result<(), ApexBuildError> {
    let v = plan.apex_vertex;
    if !brep.vertices.contains_key(v) {
        return Ok(());
    }

    if repeated_apex_normalized {
        return remove_apex_vertex_after_repeated_apex_normalization(
            brep,
            plan,
            cap_face,
            require_v_free_repeated_apex,
        );
    }

    let trace = crate::direct_apex_trace::enabled();
    if trace {
        let outgoings: Vec<_> = brep
            .half_edges
            .iter()
            .filter(|(_, h)| h.origin == v)
            .map(|(h, hd)| (h, hd.loop_, hd.twin))
            .collect();
        eprintln!(
            "[direct-apex-trace §14.8/D] enter remove_apex_vertex: V={v:?} has {} outgoings",
            outgoings.len()
        );
        for (h, l, t) in &outgoings {
            let twin_loop = brep.half_edges[*t].loop_;
            let face = brep.loops[*l].face;
            let twin_face = brep.loops[twin_loop].face;
            let dest = brep.half_edges[*t].origin;
            eprintln!(
                "[direct-apex-trace §14.8/D]   he={h:?} -> {dest:?} loop={l:?} (face={face:?}) twin_loop={twin_loop:?} (face={twin_face:?}) {}",
                if l == &twin_loop {
                    "SAME-LOOP"
                } else {
                    "cross-face"
                }
            );
        }
    }

    let max_iters = brep.half_edges.len() + 16;
    for _ in 0..max_iters {
        if !brep.vertices.contains_key(v) {
            return Ok(());
        }
        let outgoings: Vec<popkern_core::topology::HalfEdgeHandle> = brep
            .half_edges
            .iter()
            .filter(|(_, h)| h.origin == v)
            .map(|(h, _)| h)
            .collect();
        if outgoings.len() <= 1 {
            break;
        }

        let mut adjacent_same_loop: Option<popkern_core::topology::HalfEdgeHandle> = None;
        let mut slit_same_loop: Option<popkern_core::topology::HalfEdgeHandle> = None;
        let mut cross_face_in_cap: Option<popkern_core::topology::HalfEdgeHandle> = None;
        let mut cross_face_other: Option<popkern_core::topology::HalfEdgeHandle> = None;
        for &he in &outgoings {
            let twin = brep.half_edges[he].twin;
            let he_loop = brep.half_edges[he].loop_;
            let twin_loop = brep.half_edges[twin].loop_;
            if he_loop == twin_loop {
                let he_next = brep.half_edges[he].next;
                let twin_next = brep.half_edges[twin].next;
                let adjacent = he_next == twin || twin_next == he;
                if adjacent {
                    if adjacent_same_loop.is_none() {
                        adjacent_same_loop = Some(he);
                    }
                } else if slit_same_loop.is_none() {
                    slit_same_loop = Some(he);
                }
                continue;
            }
            let he_face = brep.loops[he_loop].face;
            if he_face == cap_face && cross_face_in_cap.is_none() {
                cross_face_in_cap = Some(he);
            } else if cross_face_other.is_none() {
                cross_face_other = Some(he);
            }
        }

        if let Some(he) = cross_face_in_cap {
            if trace {
                let twin = brep.half_edges[he].twin;
                let twin_face = brep.loops[brep.half_edges[twin].loop_].face;
                eprintln!(
                    "[direct-apex-trace §14.8/D] kef cross-face he={he:?} (in cap_face={cap_face:?}, killing {twin_face:?}; V outgoings={})",
                    outgoings.len()
                );
            }
            brep.kef(he).map_err(|e| ApexBuildError::EulerOpFailed {
                op: "kef",
                step: 301,
                cause: format!("V cross-face kef (cap-side): {e}"),
            })?;
        } else if let Some(he) = cross_face_other {
            let twin = brep.half_edges[he].twin;
            if trace {
                let twin_face = brep.loops[brep.half_edges[twin].loop_].face;
                eprintln!(
                    "[direct-apex-trace §14.8/D] kef cross-face TWIN he={twin:?} (in cap_face={twin_face:?}; V outgoings={})",
                    outgoings.len()
                );
            }
            brep.kef(twin).map_err(|e| ApexBuildError::EulerOpFailed {
                op: "kef",
                step: 302,
                cause: format!("V cross-face kef (twin-side): {e}"),
            })?;
        } else if let Some(bridge) = adjacent_same_loop {
            if trace {
                eprintln!(
                    "[direct-apex-trace §14.8/D] dissolve_same_loop_bridge (adjacent) he={bridge:?} (V outgoings={})",
                    outgoings.len()
                );
            }
            brep.dissolve_same_loop_bridge(bridge)
                .map_err(|e| ApexBuildError::EulerOpFailed {
                    op: "dissolve_same_loop_bridge",
                    step: 300,
                    cause: format!("V-adjacent dissolve: {e}"),
                })?;
        } else if let Some(slit) = slit_same_loop {
            if outgoings.len() > 1 {
                if trace {
                    eprintln!(
                        "[direct-apex-trace §14.8/D] kemr slit he={slit:?} (V outgoings={})",
                        outgoings.len()
                    );
                }
                brep.kemr(slit).map_err(|e| ApexBuildError::EulerOpFailed {
                    op: "kemr",
                    step: 303,
                    cause: format!("V slit kemr: {e}"),
                })?;
            } else {
                break;
            }
        } else {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "remove_apex_vertex: V has {} outgoings, no usable same-loop bridge \
                     or cross-face incidence found - pipeline corruption",
                    outgoings.len()
                ),
            });
        }
    }

    crate::construct::blend_insert::kev_v_bridge(brep, v).map_err(|e| {
        ApexBuildError::PlanMismatch {
            reason: format!("remove_apex_vertex: kev_v_bridge: {e}"),
        }
    })?;
    crate::construct::blend_insert::remove_orphaned_vertex(brep, v);
    repair_stale_shell_anchors(brep)?;
    Ok(())
}

fn remove_apex_vertex_after_repeated_apex_normalization(
    brep: &mut BRep,
    plan: &ApexPlan,
    cap_face: FaceHandle,
    require_v_free_repeated_apex: bool,
) -> Result<(), ApexBuildError> {
    let v = plan.apex_vertex;
    let trace = crate::direct_apex_trace::enabled();

    if collect_outer_loop_vertices(brep, cap_face).is_some_and(|vertices| vertices.contains(&v)) {
        return Err(ApexBuildError::PlanMismatch {
            reason: "repeated-apex normalization left V on the cap loop before apex removal"
                .to_string(),
        });
    }

    if trace {
        let outgoings = brep
            .half_edges
            .iter()
            .filter(|(_, h)| h.origin == v)
            .count();
        eprintln!(
            "[direct-apex-trace §14.9/D] repeated-apex normalized cap excludes V; preserving strip ownership with V still present in non-cap topology ({outgoings} outgoing half-edge(s))"
        );
    }

    if require_v_free_repeated_apex {
        let production_faces = production_faces_containing_vertex(brep, v);
        if !production_faces.is_empty() {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "repeated-apex construction retained apex vertex {v:?} on production \
                     surfaced face(s) {production_faces:?}; accepted construction must author \
                     V-free cylinder/support trims rather than carry retired apex cells"
                ),
            });
        }
    }

    repair_stale_shell_anchors(brep)?;
    Ok(())
}
