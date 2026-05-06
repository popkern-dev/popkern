//! Candidate enumeration for repeated-apex side-support region planning.

use popkern_core::topology::{FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::repeated_apex_region_graph::{RegionGraph, RegionTrimEndpointPair};
use super::repeated_apex_side_support_facts::EndpointRegionMembership;

#[derive(Clone, Debug)]
pub(super) struct RegionOperationCandidate {
    pub(super) fillet_index: usize,
    pub(super) kind: RegionOperationKind,
    pub(super) endpoints: RegionCandidateEndpoints,
    pub(super) status: CandidateStatus,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum RegionOperationKind {
    SplitFaceWithEdge,
    BridgeSupportRegion,
    MergeAdjacentSupportRegion,
}

#[derive(Clone, Debug)]
pub(super) struct RegionCandidateEndpoints {
    pub(super) cylinder_face: FaceHandle,
    pub(super) support_face: FaceHandle,
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
    pub(super) start_membership: EndpointRegionMembership,
    pub(super) end_membership: EndpointRegionMembership,
}

#[derive(Clone, Debug)]
pub(super) enum CandidateStatus {
    Accepted,
    Rejected { reason: CandidateRejection },
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum CandidateRejection {
    EndpointsNotInOneSupportLoop,
    EndpointMembershipUnknown,
    EndpointMembershipAmbiguous,
}

#[derive(Clone, Debug)]
pub(super) struct CandidateSetSnapshot {
    candidates: Vec<RegionOperationCandidate>,
}

impl CandidateSetSnapshot {
    pub(super) fn enumerate(region_graph: &RegionGraph) -> Self {
        let mut candidates = Vec::new();
        for pair in region_graph.trim_endpoint_pairs() {
            match (&pair.start_membership, &pair.end_membership) {
                (
                    EndpointRegionMembership::SameSupportLoop {
                        loop_handle: start_loop,
                        ..
                    },
                    EndpointRegionMembership::SameSupportLoop {
                        loop_handle: end_loop,
                        ..
                    },
                ) if start_loop == end_loop => {
                    candidates.push(RegionOperationCandidate {
                        fillet_index: pair.fillet_index,
                        kind: RegionOperationKind::SplitFaceWithEdge,
                        endpoints: RegionCandidateEndpoints::from_pair(pair),
                        status: CandidateStatus::Accepted,
                    });
                }
                (EndpointRegionMembership::Unknown, _) | (_, EndpointRegionMembership::Unknown) => {
                    candidates.push(RegionOperationCandidate {
                        fillet_index: pair.fillet_index,
                        kind: RegionOperationKind::SplitFaceWithEdge,
                        endpoints: RegionCandidateEndpoints::from_pair(pair),
                        status: CandidateStatus::Rejected {
                            reason: CandidateRejection::EndpointMembershipUnknown,
                        },
                    });
                }
                (EndpointRegionMembership::Ambiguous { .. }, _)
                | (_, EndpointRegionMembership::Ambiguous { .. }) => {
                    candidates.push(RegionOperationCandidate {
                        fillet_index: pair.fillet_index,
                        kind: RegionOperationKind::SplitFaceWithEdge,
                        endpoints: RegionCandidateEndpoints::from_pair(pair),
                        status: CandidateStatus::Rejected {
                            reason: CandidateRejection::EndpointMembershipAmbiguous,
                        },
                    });
                }
                _ => {
                    candidates.push(RegionOperationCandidate {
                        fillet_index: pair.fillet_index,
                        kind: RegionOperationKind::SplitFaceWithEdge,
                        endpoints: RegionCandidateEndpoints::from_pair(pair),
                        status: CandidateStatus::Rejected {
                            reason: CandidateRejection::EndpointsNotInOneSupportLoop,
                        },
                    });
                    candidates.push(RegionOperationCandidate {
                        fillet_index: pair.fillet_index,
                        kind: RegionOperationKind::BridgeSupportRegion,
                        endpoints: RegionCandidateEndpoints::from_pair(pair),
                        status: CandidateStatus::Accepted,
                    });
                    candidates.push(RegionOperationCandidate {
                        fillet_index: pair.fillet_index,
                        kind: RegionOperationKind::MergeAdjacentSupportRegion,
                        endpoints: RegionCandidateEndpoints::from_pair(pair),
                        status: CandidateStatus::Accepted,
                    });
                }
            }
        }
        candidates.sort_by_key(|candidate| {
            (
                candidate.fillet_index,
                candidate.kind.diagnostic_label(),
                candidate.status.diagnostic_label(),
            )
        });
        Self { candidates }
    }

    #[cfg(test)]
    pub(super) fn len(&self) -> usize {
        self.candidates.len()
    }

    pub(super) fn iter(&self) -> impl Iterator<Item = &RegionOperationCandidate> {
        self.candidates.iter()
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "candidates=[{}]",
            self.candidates
                .iter()
                .map(RegionOperationCandidate::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }
}

impl RegionOperationCandidate {
    fn diagnostic_label(&self) -> String {
        format!(
            "FilletEnd[{}] {} {} {}",
            self.fillet_index,
            self.kind.diagnostic_label(),
            self.endpoints.diagnostic_label(),
            self.status.diagnostic_label()
        )
    }
}

impl RegionOperationKind {
    fn diagnostic_label(self) -> &'static str {
        match self {
            RegionOperationKind::SplitFaceWithEdge => "split-face-with-edge",
            RegionOperationKind::BridgeSupportRegion => "bridge-support-region",
            RegionOperationKind::MergeAdjacentSupportRegion => "merge-adjacent-support-region",
        }
    }
}

impl CandidateStatus {
    fn diagnostic_label(&self) -> &'static str {
        match self {
            CandidateStatus::Accepted => "accepted",
            CandidateStatus::Rejected { reason } => reason.diagnostic_label(),
        }
    }
}

impl CandidateRejection {
    fn diagnostic_label(self) -> &'static str {
        match self {
            CandidateRejection::EndpointsNotInOneSupportLoop => {
                "rejected:endpoints-not-in-one-support-loop"
            }
            CandidateRejection::EndpointMembershipUnknown => "rejected:endpoint-membership-unknown",
            CandidateRejection::EndpointMembershipAmbiguous => {
                "rejected:endpoint-membership-ambiguous"
            }
        }
    }
}

pub(super) struct BridgeLoopMembership {
    pub(super) cylinder_face: FaceHandle,
    pub(super) support_loop: LoopHandle,
    pub(super) support_half_edge: HalfEdgeHandle,
    pub(super) outside_face: FaceHandle,
    pub(super) outside_loop: LoopHandle,
    pub(super) outside_half_edge: HalfEdgeHandle,
}

impl RegionCandidateEndpoints {
    fn from_pair(pair: &RegionTrimEndpointPair) -> Self {
        Self {
            cylinder_face: pair.cylinder_face,
            support_face: pair.support_face,
            start: pair.start,
            end: pair.end,
            start_membership: pair.start_membership.clone(),
            end_membership: pair.end_membership.clone(),
        }
    }

    pub(super) fn bridge_loop_membership(&self) -> Option<BridgeLoopMembership> {
        match (&self.start_membership, &self.end_membership) {
            (
                EndpointRegionMembership::SameSupportLoop {
                    face: _,
                    loop_handle,
                    half_edge,
                },
                EndpointRegionMembership::OutsideSupportFamily {
                    face: outside_face,
                    loop_handle: outside_loop,
                    half_edge: outside_half_edge,
                },
            ) => Some(BridgeLoopMembership {
                cylinder_face: self.cylinder_face,
                support_loop: *loop_handle,
                support_half_edge: *half_edge,
                outside_face: *outside_face,
                outside_loop: *outside_loop,
                outside_half_edge: *outside_half_edge,
            }),
            (
                EndpointRegionMembership::OutsideSupportFamily {
                    face: outside_face,
                    loop_handle: outside_loop,
                    half_edge: outside_half_edge,
                },
                EndpointRegionMembership::SameSupportLoop {
                    face: _,
                    loop_handle,
                    half_edge,
                },
            ) => Some(BridgeLoopMembership {
                cylinder_face: self.cylinder_face,
                support_loop: *loop_handle,
                support_half_edge: *half_edge,
                outside_face: *outside_face,
                outside_loop: *outside_loop,
                outside_half_edge: *outside_half_edge,
            }),
            _ => None,
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "cylinder={:?} support={:?} {:?}->{:?} start={} end={}",
            self.cylinder_face,
            self.support_face,
            self.start,
            self.end,
            self.start_membership.diagnostic_label(),
            self.end_membership.diagnostic_label()
        )
    }
}
