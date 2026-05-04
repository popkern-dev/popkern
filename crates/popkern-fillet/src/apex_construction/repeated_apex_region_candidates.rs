//! Candidate enumeration for repeated-apex side-support region planning.

use popkern_core::topology::{FaceHandle, LoopHandle, VertexHandle};

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
}

#[derive(Clone, Debug)]
pub(super) struct RegionCandidateEndpoints {
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
}

#[derive(Clone, Debug)]
pub(super) struct CandidateSetSnapshot {
    pub(super) candidates: Vec<RegionOperationCandidate>,
}

impl CandidateSetSnapshot {
    pub(super) fn enumerate(region_graph: &RegionGraph) -> Self {
        let mut candidates = Vec::new();
        for pair in &region_graph.trim_endpoint_pairs {
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
        }
    }
}

impl RegionCandidateEndpoints {
    fn from_pair(pair: &RegionTrimEndpointPair) -> Self {
        Self {
            support_face: pair.support_face,
            start: pair.start,
            end: pair.end,
            start_membership: pair.start_membership.clone(),
            end_membership: pair.end_membership.clone(),
        }
    }

    pub(super) fn support_loop_membership(&self) -> Option<(FaceHandle, LoopHandle)> {
        match (&self.start_membership, &self.end_membership) {
            (
                EndpointRegionMembership::SameSupportLoop { face, loop_handle },
                EndpointRegionMembership::OutsideSupportFamily,
            )
            | (
                EndpointRegionMembership::OutsideSupportFamily,
                EndpointRegionMembership::SameSupportLoop { face, loop_handle },
            ) => Some((*face, *loop_handle)),
            _ => None,
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "support={:?} {:?}->{:?} start={} end={}",
            self.support_face,
            self.start,
            self.end,
            self.start_membership.diagnostic_label(),
            self.end_membership.diagnostic_label()
        )
    }
}
