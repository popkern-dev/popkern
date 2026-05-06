use std::collections::BTreeSet;

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use super::*;

pub(in crate::apex_construction) fn expose_aggregate_bridge_boundaries(
    primitive: &FaceFamilyBridgePrimitivePlan,
    candidate: &mut BRep,
) -> Result<(), String> {
    let mut errors = Vec::new();
    let mut retired_faces = BTreeSet::<FaceHandle>::new();
    for region in &primitive.regions {
        for plan in &region.boundary_exposure.retire_faces {
            retired_faces.insert(plan.face);
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
    match region.boundary_exposure.status {
        BridgeBoundaryExposureStatus::AlreadyBoundaryReady => None,
        BridgeBoundaryExposureStatus::PartialFaceOwnership => {
            let mut errors = Vec::new();
            let retire_faces = region
                .boundary_exposure
                .retire_faces
                .iter()
                .map(|plan| plan.face)
                .collect::<Vec<_>>();
            if !retire_faces.is_empty()
                && let Err(err) = candidate.retire_faces_to_boundary(&retire_faces)
            {
                errors.push(format!("whole face retire rejected: {err}"));
            }
            for partial_face in &region.boundary_exposure.partial_faces {
                for pair in &partial_face.self_owned_protected_pairs {
                    let edge = candidate.half_edges[pair.half_edge].edge;
                    if !candidate.boundary_edges.contains(&edge)
                        && let Err(err) =
                            candidate.expose_same_loop_edge_as_boundary(pair.half_edge)
                    {
                        errors.push(format!(
                            "partial face {:?} self-owned pair exposure rejected: {err}",
                            partial_face.face
                        ));
                    }
                }
                if partial_face.is_fully_self_owned_pair_exposure() {
                    continue;
                }
                let Some(fragment) = &partial_face.fragment_retire else {
                    errors.push(format!(
                        "partial face {:?} has no contiguous fragment-retire interval",
                        partial_face.face
                    ));
                    continue;
                };
                if let Err(err) = candidate.retire_face_boundary_fragment(
                    partial_face.face,
                    fragment.first_half_edge,
                    fragment.after_last_half_edge,
                ) {
                    errors.push(format!(
                        "partial face {:?} fragment retire rejected: {err}",
                        partial_face.face
                    ));
                }
            }
            (!errors.is_empty()).then(|| {
                format!(
                    "{} errors=[{}]",
                    region.boundary_exposure.diagnostic_label(),
                    errors.join("; ")
                )
            })
        }
        BridgeBoundaryExposureStatus::RetireWholeFaces => {
            let retire_faces = region
                .boundary_exposure
                .retire_faces
                .iter()
                .map(|plan| plan.face)
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
    for region in &primitive.regions {
        for partial_face in &region.boundary_exposure.partial_faces {
            for pair in &partial_face.self_owned_protected_pairs {
                let Some(data) = candidate.half_edges.get(pair.half_edge) else {
                    errors.push(format!(
                        "partial face {:?} self-owned pair exposure has stale half-edge {:?}",
                        partial_face.face, pair.half_edge
                    ));
                    continue;
                };
                let edge = data.edge;
                if candidate.boundary_edges.contains(&edge)
                    || !exposed_self_owned_edges.insert(edge)
                {
                    continue;
                }
                if let Err(err) = candidate.expose_same_loop_edge_as_boundary(pair.half_edge) {
                    errors.push(format!(
                        "partial face {:?} self-owned pair exposure rejected: {err}",
                        partial_face.face
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
    for region in &primitive.regions {
        for partial_face in &region.boundary_exposure.partial_faces {
            if partial_face.is_fully_self_owned_pair_exposure() {
                continue;
            }
            let Some(fragment) = &partial_face.fragment_retire else {
                errors.push(format!(
                    "partial face {:?} has no contiguous fragment-retire interval",
                    partial_face.face
                ));
                continue;
            };
            let Some(first) = candidate.half_edges.get(fragment.first_half_edge) else {
                continue;
            };
            if candidate.boundary_edges.contains(&first.edge) {
                continue;
            }
            if let Err(err) = candidate.retire_face_boundary_fragment(
                partial_face.face,
                fragment.first_half_edge,
                fragment.after_last_half_edge,
            ) {
                let message = err.to_string();
                if message.contains("already contains a boundary edge")
                    || message.contains("after-last half-edge is invalid")
                {
                    continue;
                }
                errors.push(format!(
                    "partial face {:?} fragment retire rejected: {message}",
                    partial_face.face
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
    for region in &primitive.regions {
        for blocker in &region.boundary_exposure.blockers {
            let Some(source) = candidate.half_edges.get(blocker.source_half_edge) else {
                continue;
            };
            if candidate.boundary_edges.contains(&source.edge) {
                continue;
            }
            let Some(blocking) = candidate.half_edges.get(blocker.blocking_half_edge) else {
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
                blocker.blocking_half_edge,
                after_last,
            ) {
                let message = err.to_string();
                if message.contains("already contains a boundary edge")
                    || message.contains("interval must not cover the whole loop")
                    || message.contains("after-last half-edge is invalid")
                {
                    continue;
                }
                errors.push(format!(
                    "remaining bridge side {:?} blocked by {:?} rejected: {message}",
                    blocker.source_half_edge, blocker.blocking_half_edge
                ));
            }
        }
    }
}
