//! Candidate enumeration for repeated-apex side-support region planning.

use super::repeated_apex_region_graph::{EndpointRegionMembership, RegionGraph};

#[derive(Clone, Debug)]
pub(super) struct RegionOperationCandidate {
    pub(super) fillet_index: usize,
    pub(super) kind: RegionOperationKind,
    pub(super) status: CandidateStatus,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum RegionOperationKind {
    SplitFaceWithEdge,
    RequiresSupportRegionBridge,
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
    MissingValidatedBridgePrimitive,
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
                        status: CandidateStatus::Accepted,
                    });
                }
                (EndpointRegionMembership::Unknown, _) | (_, EndpointRegionMembership::Unknown) => {
                    candidates.push(RegionOperationCandidate {
                        fillet_index: pair.fillet_index,
                        kind: RegionOperationKind::SplitFaceWithEdge,
                        status: CandidateStatus::Rejected {
                            reason: CandidateRejection::EndpointMembershipUnknown,
                        },
                    });
                }
                _ => {
                    candidates.push(RegionOperationCandidate {
                        fillet_index: pair.fillet_index,
                        kind: RegionOperationKind::SplitFaceWithEdge,
                        status: CandidateStatus::Rejected {
                            reason: CandidateRejection::EndpointsNotInOneSupportLoop,
                        },
                    });
                    candidates.push(RegionOperationCandidate {
                        fillet_index: pair.fillet_index,
                        kind: RegionOperationKind::RequiresSupportRegionBridge,
                        status: CandidateStatus::Rejected {
                            reason: CandidateRejection::MissingValidatedBridgePrimitive,
                        },
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
            "FilletEnd[{}] {} {}",
            self.fillet_index,
            self.kind.diagnostic_label(),
            self.status.diagnostic_label()
        )
    }
}

impl RegionOperationKind {
    fn diagnostic_label(self) -> &'static str {
        match self {
            RegionOperationKind::SplitFaceWithEdge => "split-face-with-edge",
            RegionOperationKind::RequiresSupportRegionBridge => "requires-support-region-bridge",
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
            CandidateRejection::MissingValidatedBridgePrimitive => {
                "rejected:missing-validated-bridge-primitive"
            }
        }
    }
}
