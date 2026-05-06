//! Boundary-exposure planning for face-family bridge operands.

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::repeated_apex_side_support_bridge::FaceFamilyBridgeBoundaryCycle;

#[path = "repeated_apex_side_support_boundary_diagnostics.rs"]
mod boundary_diagnostics;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum FaceFamilyBridgeBoundaryReadiness {
    BoundaryReady,
    InteriorStillOwned { blocking_half_edge: HalfEdgeHandle },
}

#[derive(Clone, Debug)]
pub(super) struct FaceFamilyBridgeBoundaryExposurePlan {
    pub(super) readiness: FaceFamilyBridgeBoundaryReadiness,
    pub(super) blockers: Vec<BridgeBoundaryOwnershipBlocker>,
    pub(super) retire_faces: Vec<BridgeBoundaryRetireFacePlan>,
    pub(super) partial_faces: Vec<BridgeBoundaryPartialFacePlan>,
    pub(super) status: BridgeBoundaryExposureStatus,
}

#[derive(Clone, Debug)]
pub(super) struct BridgeBoundaryOwnershipBlocker {
    pub(super) source_half_edge: HalfEdgeHandle,
    pub(super) edge: EdgeHandle,
    pub(super) source_face: FaceHandle,
    pub(super) blocking_half_edge: HalfEdgeHandle,
    pub(super) blocking_face: FaceHandle,
}

#[derive(Clone, Debug)]
pub(super) struct BridgeBoundaryRetireFacePlan {
    pub(super) face: FaceHandle,
    pub(super) half_edges: Vec<HalfEdgeHandle>,
}

#[derive(Clone, Debug)]
pub(super) struct BridgeBoundaryPartialFacePlan {
    pub(super) face: FaceHandle,
    blocked_half_edges: Vec<HalfEdgeHandle>,
    face_half_edges: Vec<HalfEdgeHandle>,
    pub(super) fragment_retire: Option<BridgeBoundaryFragmentRetirePlan>,
    pub(super) self_owned_protected_pairs: Vec<BridgeBoundarySelfOwnedPairPlan>,
    fragment_retire_refusal: Option<BridgeBoundaryFragmentRetireRefusal>,
}

#[derive(Clone, Debug)]
pub(super) struct BridgeBoundaryFragmentRetirePlan {
    pub(super) first_half_edge: HalfEdgeHandle,
    pub(super) after_last_half_edge: HalfEdgeHandle,
}

#[derive(Clone, Debug)]
pub(super) struct BridgeBoundarySelfOwnedPairPlan {
    pub(super) half_edge: HalfEdgeHandle,
    pub(super) twin: HalfEdgeHandle,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BridgeBoundaryFragmentRetireRefusal {
    ProtectedBridgeOperandOverlap,
    SelfOwnedProtectedBridgeOperandPair,
    NonContiguousFragment,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum BridgeBoundaryExposureStatus {
    AlreadyBoundaryReady,
    RetireWholeFaces,
    PartialFaceOwnership,
}

impl FaceFamilyBridgeBoundaryExposurePlan {
    pub(super) fn build(brep: &BRep, cycles: &[FaceFamilyBridgeBoundaryCycle]) -> Self {
        let mut blockers = Vec::new();
        for cycle in cycles {
            for half_edge in cycle
                .support_interval()
                .iter()
                .chain(cycle.outside_interval().iter())
            {
                if bridge_interval_half_edge_is_boundary_ready(brep, *half_edge) {
                    continue;
                }
                if let Some(blocker) = BridgeBoundaryOwnershipBlocker::build(brep, *half_edge) {
                    blockers.push(blocker);
                } else {
                    return Self::partial(blockers, *half_edge);
                }
            }
        }

        if blockers.is_empty() {
            return Self {
                readiness: FaceFamilyBridgeBoundaryReadiness::BoundaryReady,
                blockers,
                retire_faces: Vec::new(),
                partial_faces: Vec::new(),
                status: BridgeBoundaryExposureStatus::AlreadyBoundaryReady,
            };
        }

        let protected_half_edges = bridge_cycle_half_edges(cycles);
        let retire_faces = whole_retire_faces(brep, &blockers);
        let partial_faces =
            partial_exposure_faces(brep, &blockers, &retire_faces, &protected_half_edges);
        let status = if retire_faces_ready(&retire_faces, &blockers) {
            BridgeBoundaryExposureStatus::RetireWholeFaces
        } else {
            BridgeBoundaryExposureStatus::PartialFaceOwnership
        };
        let readiness = FaceFamilyBridgeBoundaryReadiness::InteriorStillOwned {
            blocking_half_edge: blockers[0].source_half_edge,
        };
        Self {
            readiness,
            blockers,
            retire_faces,
            partial_faces,
            status,
        }
    }

    pub(super) fn is_partial_exposure_operator_ready(&self) -> bool {
        self.status == BridgeBoundaryExposureStatus::PartialFaceOwnership
            && !self.partial_faces.is_empty()
            && self.partial_faces.iter().all(|plan| {
                plan.fragment_retire.is_some() || plan.is_fully_self_owned_pair_exposure()
            })
    }

    pub(super) fn fragment_retire_count(&self) -> usize {
        self.partial_faces
            .iter()
            .filter(|plan| plan.fragment_retire.is_some())
            .count()
    }

    pub(super) fn self_owned_exposure_count(&self) -> usize {
        self.partial_faces
            .iter()
            .map(|plan| plan.self_owned_protected_pairs.len())
            .sum()
    }

    pub(super) fn has_self_owned_protected_operand_pair(&self) -> bool {
        self.partial_faces.iter().any(|plan| {
            plan.fragment_retire_refusal
                == Some(BridgeBoundaryFragmentRetireRefusal::SelfOwnedProtectedBridgeOperandPair)
        })
    }

    fn partial(
        blockers: Vec<BridgeBoundaryOwnershipBlocker>,
        blocking_half_edge: HalfEdgeHandle,
    ) -> Self {
        Self {
            readiness: FaceFamilyBridgeBoundaryReadiness::InteriorStillOwned { blocking_half_edge },
            blockers,
            retire_faces: Vec::new(),
            partial_faces: Vec::new(),
            status: BridgeBoundaryExposureStatus::PartialFaceOwnership,
        }
    }
}

impl BridgeBoundaryOwnershipBlocker {
    fn build(brep: &BRep, source_half_edge: HalfEdgeHandle) -> Option<Self> {
        let source = brep.half_edges.get(source_half_edge)?;
        let blocking_half_edge = source.twin;
        let blocking = brep.half_edges.get(blocking_half_edge)?;
        Some(Self {
            source_half_edge,
            edge: source.edge,
            source_face: brep.loops.get(source.loop_)?.face,
            blocking_half_edge,
            blocking_face: brep.loops.get(blocking.loop_)?.face,
        })
    }
}

impl BridgeBoundaryPartialFacePlan {
    pub(super) fn is_fully_self_owned_pair_exposure(&self) -> bool {
        self.fragment_retire.is_none()
            && !self.self_owned_protected_pairs.is_empty()
            && self.blocked_half_edges.len() == self.self_owned_protected_pairs.len() * 2
    }
}

fn bridge_interval_half_edge_is_boundary_ready(brep: &BRep, half_edge: HalfEdgeHandle) -> bool {
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return false;
    };
    brep.boundary_edges.contains(&half_edge_data.edge)
        && !brep.half_edges.contains_key(half_edge_data.twin)
}

fn whole_retire_faces(
    brep: &BRep,
    blockers: &[BridgeBoundaryOwnershipBlocker],
) -> Vec<BridgeBoundaryRetireFacePlan> {
    let mut blocked_by_face = BTreeMap::<FaceHandle, BTreeSet<HalfEdgeHandle>>::new();
    for blocker in blockers {
        blocked_by_face
            .entry(blocker.blocking_face)
            .or_default()
            .insert(blocker.blocking_half_edge);
    }
    blocked_by_face
        .into_iter()
        .filter_map(|(face, blocked)| {
            let face_half_edges = face_outer_half_edges(brep, face)?;
            let face_set = face_half_edges.iter().copied().collect::<BTreeSet<_>>();
            (face_set == blocked).then_some(BridgeBoundaryRetireFacePlan {
                face,
                half_edges: face_half_edges,
            })
        })
        .collect()
}

fn retire_faces_ready(
    retire_faces: &[BridgeBoundaryRetireFacePlan],
    blockers: &[BridgeBoundaryOwnershipBlocker],
) -> bool {
    let planned_faces = retire_faces
        .iter()
        .map(|plan| plan.face)
        .collect::<BTreeSet<_>>();
    blockers
        .iter()
        .all(|blocker| planned_faces.contains(&blocker.blocking_face))
}

fn partial_exposure_faces(
    brep: &BRep,
    blockers: &[BridgeBoundaryOwnershipBlocker],
    retire_faces: &[BridgeBoundaryRetireFacePlan],
    protected_half_edges: &BTreeSet<HalfEdgeHandle>,
) -> Vec<BridgeBoundaryPartialFacePlan> {
    let whole_retire_faces = retire_faces
        .iter()
        .map(|plan| plan.face)
        .collect::<BTreeSet<_>>();
    let mut blocked_by_face = BTreeMap::<FaceHandle, BTreeSet<HalfEdgeHandle>>::new();
    for blocker in blockers {
        if whole_retire_faces.contains(&blocker.blocking_face) {
            continue;
        }
        blocked_by_face
            .entry(blocker.blocking_face)
            .or_default()
            .insert(blocker.blocking_half_edge);
    }
    blocked_by_face
        .into_iter()
        .filter_map(|(face, blocked)| {
            let face_half_edges = face_outer_half_edges(brep, face)?;
            let self_owned_protected_pairs =
                self_owned_protected_pairs(brep, &blocked, protected_half_edges);
            let fragment_retire = if blocked.is_disjoint(protected_half_edges) {
                contiguous_fragment_retire_plan(&face_half_edges, &blocked)
            } else {
                None
            };
            let fragment_retire_refusal = if fragment_retire.is_some() {
                None
            } else if !self_owned_protected_pairs.is_empty() {
                Some(BridgeBoundaryFragmentRetireRefusal::SelfOwnedProtectedBridgeOperandPair)
            } else if !blocked.is_disjoint(protected_half_edges) {
                Some(BridgeBoundaryFragmentRetireRefusal::ProtectedBridgeOperandOverlap)
            } else {
                Some(BridgeBoundaryFragmentRetireRefusal::NonContiguousFragment)
            };
            Some(BridgeBoundaryPartialFacePlan {
                face,
                blocked_half_edges: blocked.into_iter().collect(),
                face_half_edges,
                fragment_retire,
                self_owned_protected_pairs,
                fragment_retire_refusal,
            })
        })
        .collect()
}

fn self_owned_protected_pairs(
    brep: &BRep,
    blocked: &BTreeSet<HalfEdgeHandle>,
    protected_half_edges: &BTreeSet<HalfEdgeHandle>,
) -> Vec<BridgeBoundarySelfOwnedPairPlan> {
    blocked
        .iter()
        .filter_map(|half_edge| {
            let Some(data) = brep.half_edges.get(*half_edge) else {
                return None;
            };
            let twin = data.twin;
            if *half_edge > twin {
                return None;
            }
            (protected_half_edges.contains(half_edge)
                && protected_half_edges.contains(&twin)
                && brep
                    .half_edges
                    .get(twin)
                    .is_some_and(|twin_data| twin_data.loop_ == data.loop_))
            .then_some(BridgeBoundarySelfOwnedPairPlan {
                half_edge: *half_edge,
                twin,
            })
        })
        .collect()
}

fn bridge_cycle_half_edges(cycles: &[FaceFamilyBridgeBoundaryCycle]) -> BTreeSet<HalfEdgeHandle> {
    cycles
        .iter()
        .flat_map(|cycle| {
            cycle
                .support_interval()
                .iter()
                .chain(cycle.outside_interval().iter())
        })
        .copied()
        .collect()
}

fn contiguous_fragment_retire_plan(
    face_half_edges: &[HalfEdgeHandle],
    blocked: &BTreeSet<HalfEdgeHandle>,
) -> Option<BridgeBoundaryFragmentRetirePlan> {
    if blocked.is_empty() || blocked.len() == face_half_edges.len() {
        return None;
    }
    let is_blocked = face_half_edges
        .iter()
        .map(|half_edge| blocked.contains(half_edge))
        .collect::<Vec<_>>();
    let mut starts = Vec::new();
    for i in 0..is_blocked.len() {
        let prev = if i == 0 { is_blocked.len() - 1 } else { i - 1 };
        if is_blocked[i] && !is_blocked[prev] {
            starts.push(i);
        }
    }
    if starts.len() != 1 {
        return None;
    }
    let first_index = starts[0];
    let mut after_last_index = (first_index + 1) % is_blocked.len();
    while is_blocked[after_last_index] {
        after_last_index = (after_last_index + 1) % is_blocked.len();
        if after_last_index == first_index {
            return None;
        }
    }
    Some(BridgeBoundaryFragmentRetirePlan {
        first_half_edge: face_half_edges[first_index],
        after_last_half_edge: face_half_edges[after_last_index],
    })
}

fn face_outer_half_edges(brep: &BRep, face: FaceHandle) -> Option<Vec<HalfEdgeHandle>> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    let mut half_edges = Vec::new();
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.loop_ != face_data.outer_loop {
            return None;
        }
        half_edges.push(cursor);
        cursor = half_edge.next;
        if cursor == start {
            return Some(half_edges);
        }
    }
    None
}
