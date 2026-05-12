//! Patch candidate diagnostics for high-valence and vertex-patch strategies.
#![allow(missing_docs)]

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PatchCandidateApplicability {
    Applicable,
    NotApplicable,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PatchCandidateFeasibility {
    Ready,
    MissingBoundaryPrimitive(&'static str),
    MissingPatchPrimitive(&'static str),
    Invalid(&'static str),
    NotApplicable,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PatchCandidateQuality {
    Acceptable,
    ManagedRisk(&'static str),
    InvalidResidual(&'static str),
    InvalidTrim(&'static str),
    NotApplicable,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PatchCandidateDowngradePolicy {
    PreserveRequestedContinuity,
    DowngradeAllowed(&'static str),
    NotApplicable,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PatchCandidateFact {
    pub applicability: PatchCandidateApplicability,
    pub n_sided_feasibility: PatchCandidateFeasibility,
    pub vertex_patch_feasibility: PatchCandidateFeasibility,
    pub boundary_compatibility: PatchCandidateFeasibility,
    pub continuity_target: &'static str,
    pub quality: PatchCandidateQuality,
    pub downgrade_policy: PatchCandidateDowngradePolicy,
}

impl PatchCandidateFact {
    pub const NOT_APPLICABLE: Self = Self {
        applicability: PatchCandidateApplicability::NotApplicable,
        n_sided_feasibility: PatchCandidateFeasibility::NotApplicable,
        vertex_patch_feasibility: PatchCandidateFeasibility::NotApplicable,
        boundary_compatibility: PatchCandidateFeasibility::NotApplicable,
        continuity_target: "not-applicable",
        quality: PatchCandidateQuality::NotApplicable,
        downgrade_policy: PatchCandidateDowngradePolicy::NotApplicable,
    };
}
