use super::*;

impl FaceFamilyBridgeBoundaryExposurePlan {
    pub(in crate::apex_construction) fn diagnostic_label(&self) -> String {
        format!(
            "exposure_status={} blockers=[{}] retire_faces=[{}] partial_faces=[{}]",
            self.status.diagnostic_label(),
            self.blockers
                .iter()
                .map(BridgeBoundaryOwnershipBlocker::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.retire_faces
                .iter()
                .map(BridgeBoundaryRetireFacePlan::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.partial_faces
                .iter()
                .map(BridgeBoundaryPartialFacePlan::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; ")
        )
    }
}

impl FaceFamilyBridgeBoundaryReadiness {
    pub(in crate::apex_construction) fn is_boundary_ready(&self) -> bool {
        matches!(self, FaceFamilyBridgeBoundaryReadiness::BoundaryReady)
    }

    pub(in crate::apex_construction) fn diagnostic_label(&self) -> String {
        match self {
            FaceFamilyBridgeBoundaryReadiness::BoundaryReady => "boundary-ready".to_string(),
            FaceFamilyBridgeBoundaryReadiness::InteriorStillOwned { blocking_half_edge } => {
                format!("interior-still-owned blocking_half_edge={blocking_half_edge:?}")
            }
        }
    }
}

impl BridgeBoundaryOwnershipBlocker {
    fn diagnostic_label(&self) -> String {
        format!(
            "source={:?}/{:?} edge={:?} blocked_by={:?}/{:?}",
            self.source_face,
            self.source_half_edge,
            self.edge,
            self.blocking_face,
            self.blocking_half_edge
        )
    }
}

impl BridgeBoundaryRetireFacePlan {
    fn diagnostic_label(&self) -> String {
        format!("{:?} half_edges={:?}", self.face, self.half_edges)
    }
}

impl BridgeBoundaryPartialFacePlan {
    fn diagnostic_label(&self) -> String {
        format!(
            "{:?} blocked={}/{} blocked_half_edges={:?} face_half_edges={:?} self_owned_pairs=[{}] fragment_retire={} fragment_retire_refusal={}",
            self.face,
            self.blocked_half_edges.len(),
            self.face_half_edges.len(),
            self.blocked_half_edges,
            self.face_half_edges,
            self.self_owned_protected_pairs
                .iter()
                .map(BridgeBoundarySelfOwnedPairPlan::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.fragment_retire
                .as_ref()
                .map(BridgeBoundaryFragmentRetirePlan::diagnostic_label)
                .unwrap_or_else(|| "none".to_string()),
            self.fragment_retire_refusal
                .map(BridgeBoundaryFragmentRetireRefusal::diagnostic_label)
                .unwrap_or("none")
        )
    }
}

impl BridgeBoundarySelfOwnedPairPlan {
    fn diagnostic_label(&self) -> String {
        format!(
            "{:?}<->{:?} delta=dV=0 dE=1 dF=0 dL=0",
            self.half_edge, self.twin
        )
    }
}

impl BridgeBoundaryFragmentRetirePlan {
    fn diagnostic_label(&self) -> String {
        format!(
            "first={:?} after_last={:?} delta=dV=0 dE=1 dF=0 dL=0",
            self.first_half_edge, self.after_last_half_edge
        )
    }
}

impl BridgeBoundaryFragmentRetireRefusal {
    fn diagnostic_label(self) -> &'static str {
        match self {
            BridgeBoundaryFragmentRetireRefusal::ProtectedBridgeOperandOverlap => {
                "protected-bridge-operand-overlap"
            }
            BridgeBoundaryFragmentRetireRefusal::SelfOwnedProtectedBridgeOperandPair => {
                "self-owned-protected-bridge-operand-pair"
            }
            BridgeBoundaryFragmentRetireRefusal::NonContiguousFragment => "non-contiguous-fragment",
        }
    }
}

impl BridgeBoundaryExposureStatus {
    fn diagnostic_label(self) -> &'static str {
        match self {
            BridgeBoundaryExposureStatus::AlreadyBoundaryReady => "already-boundary-ready",
            BridgeBoundaryExposureStatus::RetireWholeFaces => "retire-whole-faces",
            BridgeBoundaryExposureStatus::PartialFaceOwnership => "partial-face-ownership",
        }
    }
}
