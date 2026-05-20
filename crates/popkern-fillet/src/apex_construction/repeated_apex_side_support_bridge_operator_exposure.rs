use std::collections::{BTreeMap, BTreeSet};

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::bridge_operator_exposure_walk::{
    aggregate_bridge_half_edges, bridge_half_edge_is_boundary_ready,
    is_non_blocking_exposure_error, live_contiguous_fragments, live_face_outer_half_edges,
};
use super::*;

pub(in crate::apex_construction) fn expose_aggregate_bridge_boundaries(
    primitive: &FaceFamilyBridgePrimitivePlan,
    candidate: &mut BRep,
) -> Result<(), String> {
    let mut errors = Vec::new();
    let mut retired_faces = BTreeSet::<FaceHandle>::new();
    for region in primitive.regions() {
        for plan in region.boundary_exposure().retire_faces() {
            retired_faces.insert(plan.face());
        }
    }
    if !retired_faces.is_empty()
        && let Err(err) =
            candidate.retire_faces_to_boundary(&retired_faces.iter().copied().collect::<Vec<_>>())
    {
        errors.push(format!("aggregate whole face retire rejected: {err}"));
    }

    expose_self_owned_pairs_once(primitive, candidate, &mut errors);
    expose_planned_partial_fragments(primitive, candidate, &mut errors);
    complete_remaining_owned_bridge_sides(primitive, candidate, &mut errors);
    complete_live_owned_bridge_sides(primitive, candidate, &mut errors);

    if errors.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "aggregate bridge boundary exposure rejected: [{}]",
            errors.join("; ")
        ))
    }
}

pub(in crate::apex_construction) fn expose_bridge_boundaries(
    region: &FaceFamilyBridgeRegionPrimitivePlan,
    candidate: &mut BRep,
) -> Option<String> {
    match region.boundary_exposure().status() {
        BridgeBoundaryExposureStatus::AlreadyBoundaryReady => None,
        BridgeBoundaryExposureStatus::PartialFaceOwnership => {
            let mut errors = Vec::new();
            let retire_faces = region
                .boundary_exposure()
                .retire_faces()
                .iter()
                .map(|plan| plan.face())
                .collect::<Vec<_>>();
            if !retire_faces.is_empty()
                && let Err(err) = candidate.retire_faces_to_boundary(&retire_faces)
            {
                errors.push(format!("whole face retire rejected: {err}"));
            }
            for partial_face in region.boundary_exposure().partial_faces() {
                for pair in partial_face.self_owned_protected_pairs() {
                    let edge = candidate.half_edges[pair.half_edge()].edge;
                    if !candidate.boundary_edges.contains(&edge)
                        && let Err(err) =
                            candidate.expose_same_loop_edge_as_boundary(pair.half_edge())
                    {
                        errors.push(format!(
                            "partial face {:?} self-owned pair exposure rejected: {err}",
                            partial_face.face()
                        ));
                    }
                }
                if partial_face.is_fully_self_owned_pair_exposure() {
                    continue;
                }
                let Some(fragment) = partial_face.fragment_retire() else {
                    errors.push(format!(
                        "partial face {:?} has no contiguous fragment-retire interval",
                        partial_face.face()
                    ));
                    continue;
                };
                if let Err(err) = candidate.retire_face_boundary_fragment(
                    partial_face.face(),
                    fragment.first_half_edge(),
                    fragment.after_last_half_edge(),
                ) {
                    errors.push(format!(
                        "partial face {:?} fragment retire rejected: {err}",
                        partial_face.face()
                    ));
                }
            }
            (!errors.is_empty()).then(|| {
                format!(
                    "{} errors=[{}]",
                    region.boundary_exposure().diagnostic_label(),
                    errors.join("; ")
                )
            })
        }
        BridgeBoundaryExposureStatus::RetireWholeFaces => {
            let retire_faces = region
                .boundary_exposure()
                .retire_faces()
                .iter()
                .map(|plan| plan.face())
                .collect::<Vec<_>>();
            candidate
                .retire_faces_to_boundary(&retire_faces)
                .err()
                .map(|err| format!("{err}"))
        }
    }
}

fn expose_self_owned_pairs_once(
    primitive: &FaceFamilyBridgePrimitivePlan,
    candidate: &mut BRep,
    errors: &mut Vec<String>,
) {
    let mut exposed_self_owned_edges = BTreeSet::<EdgeHandle>::new();
    for region in primitive.regions() {
        for partial_face in region.boundary_exposure().partial_faces() {
            for pair in partial_face.self_owned_protected_pairs() {
                let Some(data) = candidate.half_edges.get(pair.half_edge()) else {
                    errors.push(format!(
                        "partial face {:?} self-owned pair exposure has stale half-edge {:?}",
                        partial_face.face(),
                        pair.half_edge()
                    ));
                    continue;
                };
                let edge = data.edge;
                if candidate.boundary_edges.contains(&edge)
                    || !exposed_self_owned_edges.insert(edge)
                {
                    continue;
                }
                if let Err(err) = candidate.expose_same_loop_edge_as_boundary(pair.half_edge()) {
                    errors.push(format!(
                        "partial face {:?} self-owned pair exposure rejected: {err}",
                        partial_face.face()
                    ));
                }
            }
        }
    }
}

fn expose_planned_partial_fragments(
    primitive: &FaceFamilyBridgePrimitivePlan,
    candidate: &mut BRep,
    errors: &mut Vec<String>,
) {
    for region in primitive.regions() {
        for partial_face in region.boundary_exposure().partial_faces() {
            if partial_face.is_fully_self_owned_pair_exposure() {
                continue;
            }
            let Some(fragment) = partial_face.fragment_retire() else {
                errors.push(format!(
                    "partial face {:?} has no contiguous fragment-retire interval",
                    partial_face.face()
                ));
                continue;
            };
            let Some(first) = candidate.half_edges.get(fragment.first_half_edge()) else {
                continue;
            };
            if candidate.boundary_edges.contains(&first.edge) {
                continue;
            }
            if let Err(err) = candidate.retire_face_boundary_fragment(
                partial_face.face(),
                fragment.first_half_edge(),
                fragment.after_last_half_edge(),
            ) {
                let message = err.to_string();
                if is_non_blocking_exposure_error(&message) {
                    continue;
                }
                errors.push(format!(
                    "partial face {:?} fragment retire rejected: {message}",
                    partial_face.face()
                ));
            }
        }
    }
}

fn complete_remaining_owned_bridge_sides(
    primitive: &FaceFamilyBridgePrimitivePlan,
    candidate: &mut BRep,
    errors: &mut Vec<String>,
) {
    for region in primitive.regions() {
        for blocker in region.boundary_exposure().blockers() {
            let Some(source) = candidate.half_edges.get(blocker.source_half_edge()) else {
                continue;
            };
            if candidate.boundary_edges.contains(&source.edge) {
                continue;
            }
            let Some(blocking) = candidate.half_edges.get(blocker.blocking_half_edge()) else {
                continue;
            };
            let blocking_face = candidate
                .loops
                .get(blocking.loop_)
                .map(|loop_data| loop_data.face);
            let after_last = blocking.next;
            let Some(blocking_face) = blocking_face else {
                continue;
            };
            if let Err(err) = candidate.retire_face_boundary_fragment(
                blocking_face,
                blocker.blocking_half_edge(),
                after_last,
            ) {
                let message = err.to_string();
                if is_non_blocking_exposure_error(&message)
                    || message.contains("interval must not cover the whole loop")
                {
                    continue;
                }
                errors.push(format!(
                    "remaining bridge side {:?} blocked by {:?} rejected: {message}",
                    blocker.source_half_edge(),
                    blocker.blocking_half_edge()
                ));
            }
        }
    }
}

fn complete_live_owned_bridge_sides(
    primitive: &FaceFamilyBridgePrimitivePlan,
    candidate: &mut BRep,
    errors: &mut Vec<String>,
) {
    let protected_half_edges = aggregate_bridge_half_edges(primitive);
    let max_passes = protected_half_edges.len().saturating_add(1);
    for _ in 0..max_passes {
        let blockers = live_bridge_blockers(primitive, candidate, errors);
        if blockers.is_empty() {
            return;
        }
        let mut made_progress = false;
        for (face, blocked) in blockers {
            match expose_live_blocker_interval(candidate, face, &blocked, &protected_half_edges) {
                LiveExposureProgress::Changed => made_progress = true,
                LiveExposureProgress::Unchanged => {}
                LiveExposureProgress::Rejected(message) => errors.push(message),
            }
        }
        if !made_progress {
            break;
        }
    }

    for half_edge in aggregate_bridge_half_edges(primitive) {
        if !bridge_half_edge_is_boundary_ready(candidate, half_edge) {
            errors.push(format!(
                "aggregate bridge half-edge {:?} remains owned after live exposure completion",
                half_edge
            ));
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum LiveExposureProgress {
    Changed,
    Unchanged,
    Rejected(String),
}

fn live_bridge_blockers(
    primitive: &FaceFamilyBridgePrimitivePlan,
    candidate: &BRep,
    errors: &mut Vec<String>,
) -> BTreeMap<FaceHandle, BTreeSet<HalfEdgeHandle>> {
    let mut blockers_by_face = BTreeMap::<FaceHandle, BTreeSet<HalfEdgeHandle>>::new();
    for half_edge in aggregate_bridge_half_edges(primitive) {
        if bridge_half_edge_is_boundary_ready(candidate, half_edge) {
            continue;
        }
        let Some(source) = candidate.half_edges.get(half_edge) else {
            errors.push(format!(
                "aggregate bridge half-edge {:?} is stale during boundary exposure",
                half_edge
            ));
            continue;
        };
        let blocking_half_edge = source.twin;
        let Some(blocking) = candidate.half_edges.get(blocking_half_edge) else {
            errors.push(format!(
                "aggregate bridge half-edge {:?} is not boundary but has no live blocking twin",
                half_edge
            ));
            continue;
        };
        let Some(blocking_face) = candidate
            .loops
            .get(blocking.loop_)
            .map(|loop_data| loop_data.face)
        else {
            errors.push(format!(
                "aggregate bridge half-edge {:?} is blocked by stale-loop half-edge {:?}",
                half_edge, blocking_half_edge
            ));
            continue;
        };
        blockers_by_face
            .entry(blocking_face)
            .or_default()
            .insert(blocking_half_edge);
    }
    blockers_by_face
}

fn expose_live_blocker_interval(
    candidate: &mut BRep,
    face: FaceHandle,
    blocked: &BTreeSet<HalfEdgeHandle>,
    protected_half_edges: &BTreeSet<HalfEdgeHandle>,
) -> LiveExposureProgress {
    let mut blocked = blocked
        .iter()
        .filter(|half_edge| {
            candidate
                .half_edges
                .get(**half_edge)
                .is_some_and(|data| !candidate.boundary_edges.contains(&data.edge))
        })
        .copied()
        .collect::<BTreeSet<_>>();
    if blocked.is_empty() {
        return LiveExposureProgress::Unchanged;
    }

    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace §14.9/bridge-exposure] live blockers on face {:?}: {:?}",
            face, blocked
        );
    }

    for half_edge in blocked.clone() {
        if !protected_half_edges.contains(&half_edge) {
            continue;
        }
        let Some(data) = candidate.half_edges.get(half_edge) else {
            blocked.remove(&half_edge);
            continue;
        };
        let twin = data.twin;
        if !protected_half_edges.contains(&twin) || data.loop_ != candidate.half_edges[twin].loop_ {
            return LiveExposureProgress::Rejected(format!(
                "live bridge boundary exposure would retire protected half-edge {:?} on face {:?}",
                half_edge, face
            ));
        }
        if let Err(err) = candidate.expose_same_loop_edge_as_boundary(half_edge) {
            return LiveExposureProgress::Rejected(format!(
                "live bridge boundary self-owned exposure rejected for {:?}: {err}",
                half_edge
            ));
        }
        blocked.remove(&half_edge);
        blocked.remove(&twin);
    }
    if blocked.is_empty() {
        return LiveExposureProgress::Changed;
    }

    let Some(face_half_edges) = live_face_outer_half_edges(candidate, face) else {
        return LiveExposureProgress::Rejected(format!(
            "live bridge boundary exposure could not walk face {:?}",
            face
        ));
    };
    if blocked.len() == face_half_edges.len() {
        return match candidate.retire_faces_to_boundary(&[face]) {
            Ok(_) => LiveExposureProgress::Changed,
            Err(err) => LiveExposureProgress::Rejected(format!(
                "live bridge boundary whole-face retire rejected on face {:?}: {err}",
                face
            )),
        };
    }
    let fragments = live_contiguous_fragments(&face_half_edges, &blocked);
    if crate::direct_apex_trace::enabled() {
        eprintln!(
            "[direct-apex-trace §14.9/bridge-exposure] face {:?} live loop {:?} fragments {:?}",
            face, face_half_edges, fragments
        );
    }
    let Some((first, after_last)) = fragments.first().copied() else {
        return LiveExposureProgress::Unchanged;
    };
    match candidate.retire_face_boundary_fragment(face, first, after_last) {
        Ok(_) => {
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace §14.9/bridge-exposure] retired live fragment {:?}->{:?} on face {:?}",
                    first, after_last, face
                );
            }
            LiveExposureProgress::Changed
        }
        Err(err) => {
            let message = err.to_string();
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace §14.9/bridge-exposure] live fragment {:?}->{:?} on face {:?} rejected: {message}",
                    first, after_last, face
                );
            }
            if is_non_blocking_exposure_error(&message)
                || message.contains("interval must not cover the whole loop")
            {
                LiveExposureProgress::Unchanged
            } else {
                LiveExposureProgress::Rejected(format!(
                    "live bridge boundary fragment retire rejected on face {:?}: {message}",
                    face
                ))
            }
        }
    }
}
