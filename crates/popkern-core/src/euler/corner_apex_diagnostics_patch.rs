//! Patch candidate classifiers for corner/apex diagnostics.
#![allow(missing_docs)]

use super::patch_types::{
    PatchCandidateApplicability, PatchCandidateDowngradePolicy, PatchCandidateFact,
    PatchCandidateFeasibility, PatchCandidateQuality,
};
use super::types::CornerApexStrategyFamily;

pub(super) fn patch_candidate_fact(family: CornerApexStrategyFamily) -> PatchCandidateFact {
    match family {
        CornerApexStrategyFamily::HighValencePatch => PatchCandidateFact {
            applicability: PatchCandidateApplicability::Applicable,
            n_sided_feasibility: PatchCandidateFeasibility::MissingPatchPrimitive(
                "n-sided-corner-patch",
            ),
            vertex_patch_feasibility: PatchCandidateFeasibility::MissingPatchPrimitive(
                "vertex-patch-upgrade",
            ),
            boundary_compatibility: PatchCandidateFeasibility::MissingBoundaryPrimitive(
                "patch-boundary-compatibility",
            ),
            continuity_target: "G1",
            quality: PatchCandidateQuality::ManagedRisk("high-valence-patch-quality"),
            downgrade_policy: PatchCandidateDowngradePolicy::DowngradeAllowed("G1-continuity"),
        },
        CornerApexStrategyFamily::VertexPatch => PatchCandidateFact {
            applicability: PatchCandidateApplicability::Applicable,
            n_sided_feasibility: PatchCandidateFeasibility::NotApplicable,
            vertex_patch_feasibility: PatchCandidateFeasibility::MissingPatchPrimitive(
                "vertex-patch-construction",
            ),
            boundary_compatibility: PatchCandidateFeasibility::MissingBoundaryPrimitive(
                "patch-boundary-compatibility",
            ),
            continuity_target: "G1",
            quality: PatchCandidateQuality::ManagedRisk("vertex-patch-quality"),
            downgrade_policy: PatchCandidateDowngradePolicy::DowngradeAllowed("G1-continuity"),
        },
        _ => PatchCandidateFact::NOT_APPLICABLE,
    }
}
