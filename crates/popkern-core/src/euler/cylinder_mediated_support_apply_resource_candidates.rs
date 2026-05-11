//! Typed construction candidates for cylinder-mediated resource reconciliation.

use super::apply_cylinder_resources::BridgeResourceReconciliationStatus;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct CylinderResourceCandidateSummary {
    pub(super) status: BridgeResourceReconciliationStatus,
    pub(super) exact_pair_actions: usize,
    pub(super) shared_outside_actions: usize,
    pub(super) single_covered_actions: usize,
    pub(super) unsupported_actions: usize,
}

impl CylinderResourceCandidateSummary {
    pub(super) fn diagnostic_label(self) -> String {
        let candidate = CylinderResourceConstructionCandidate::from_summary(self);
        format!(
            "cylinder_resource_construction_candidate=family={} feasibility={} reason={} action_candidates=exact_pair:{} shared_outside:{} single_covered:{} unsupported:{}",
            candidate.family.diagnostic_label(),
            candidate.feasibility.diagnostic_label(),
            candidate.reason,
            self.exact_pair_actions,
            self.shared_outside_actions,
            self.single_covered_actions,
            self.unsupported_actions
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CylinderResourceConstructionCandidate {
    family: CylinderResourceConstructionFamily,
    feasibility: CylinderResourceConstructionFeasibility,
    reason: &'static str,
}

impl CylinderResourceConstructionCandidate {
    fn from_summary(summary: CylinderResourceCandidateSummary) -> Self {
        match summary.status {
            BridgeResourceReconciliationStatus::ExactResourcePairing => Self {
                family: CylinderResourceConstructionFamily::ExactResourcePairing,
                feasibility: CylinderResourceConstructionFeasibility::ExecutableByExistingResources,
                reason: "exact-cylinder-resource-pairing-ready",
            },
            BridgeResourceReconciliationStatus::SharedOutsideTransitionNeedsResource => Self {
                family: CylinderResourceConstructionFamily::SharedOutsideTransitionResource,
                feasibility: CylinderResourceConstructionFeasibility::RequiresResourceAuthoring,
                reason: "shared-outside-transition-needs-cylinder-resource",
            },
            BridgeResourceReconciliationStatus::SingleSupportCoveredBySharedCylinderResource => {
                Self {
                    family:
                        CylinderResourceConstructionFamily::SingleSupportCoveredBySharedCylinder,
                    feasibility: CylinderResourceConstructionFeasibility::RequiresResourceSplit,
                    reason: "single-support-covered-by-shared-cylinder-resource",
                }
            }
            BridgeResourceReconciliationStatus::SharedStartSplitCylinderResourcesReady => Self {
                family: CylinderResourceConstructionFamily::SharedStartSplitCylinderResources,
                feasibility: CylinderResourceConstructionFeasibility::ExecutableAfterNormalization,
                reason: "shared-start-split-cylinder-resources-ready",
            },
            BridgeResourceReconciliationStatus::SplitCylinderResourcesNeedMerge => Self {
                family: CylinderResourceConstructionFamily::SplitCylinderResourceMerge,
                feasibility: CylinderResourceConstructionFeasibility::RequiresResourceMerge,
                reason: "split-cylinder-resource-merge-required",
            },
            BridgeResourceReconciliationStatus::IncompatibleResourceSets => Self {
                family: CylinderResourceConstructionFamily::IncompatibleResourceSets,
                feasibility: CylinderResourceConstructionFeasibility::InvalidResourceSets,
                reason: "incompatible-cylinder-resource-sets",
            },
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CylinderResourceConstructionFamily {
    ExactResourcePairing,
    SharedOutsideTransitionResource,
    SingleSupportCoveredBySharedCylinder,
    SharedStartSplitCylinderResources,
    SplitCylinderResourceMerge,
    IncompatibleResourceSets,
}

impl CylinderResourceConstructionFamily {
    fn diagnostic_label(self) -> &'static str {
        match self {
            Self::ExactResourcePairing => "exact-resource-pairing",
            Self::SharedOutsideTransitionResource => "shared-outside-transition-resource",
            Self::SingleSupportCoveredBySharedCylinder => {
                "single-support-covered-by-shared-cylinder-resource"
            }
            Self::SharedStartSplitCylinderResources => "shared-start-split-cylinder-resources",
            Self::SplitCylinderResourceMerge => "split-cylinder-resource-merge",
            Self::IncompatibleResourceSets => "incompatible-resource-sets",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CylinderResourceConstructionFeasibility {
    ExecutableByExistingResources,
    ExecutableAfterNormalization,
    RequiresResourceAuthoring,
    RequiresResourceSplit,
    RequiresResourceMerge,
    InvalidResourceSets,
}

impl CylinderResourceConstructionFeasibility {
    fn diagnostic_label(self) -> &'static str {
        match self {
            Self::ExecutableByExistingResources => "executable-by-existing-resources",
            Self::ExecutableAfterNormalization => "executable-after-normalization",
            Self::RequiresResourceAuthoring => "requires-resource-authoring",
            Self::RequiresResourceSplit => "requires-resource-split",
            Self::RequiresResourceMerge => "requires-resource-merge",
            Self::InvalidResourceSets => "invalid-resource-sets",
        }
    }
}
