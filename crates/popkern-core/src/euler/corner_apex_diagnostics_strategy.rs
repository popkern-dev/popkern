//! Strategy metadata classifiers for corner/apex diagnostics.
#![allow(missing_docs)]

use super::strategy_types::{
    AdvancedStrategyMode, CapRunoutMode, CliffTerminationMode, OverflowHandlingMode,
    StrategyMetadataFact, StrategyOptionReadiness,
};
use super::types::{BlendTypeCompatibility, CornerApexStrategyFamily, VertexClassificationFact};

pub(super) fn strategy_metadata_fact(
    family: CornerApexStrategyFamily,
    classification: VertexClassificationFact,
) -> StrategyMetadataFact {
    StrategyMetadataFact {
        cap_runout: option_state(family, CornerApexStrategyFamily::CapRunout, "not-free-end"),
        cap_runout_mode: if family == CornerApexStrategyFamily::CapRunout {
            CapRunoutMode::EndCap
        } else {
            CapRunoutMode::NotApplicable
        },
        cliff_termination: option_state(
            family,
            CornerApexStrategyFamily::CliffTermination,
            "support-not-cliff",
        ),
        cliff_termination_mode: if family == CornerApexStrategyFamily::CliffTermination {
            CliffTerminationMode::ManagedCliff
        } else {
            CliffTerminationMode::NotApplicable
        },
        mutual_trim: option_state(
            family,
            CornerApexStrategyFamily::MutualTrim,
            "not-through-vertex",
        ),
        three_face_blend: option_state(
            family,
            CornerApexStrategyFamily::ThreeFaceBlend,
            "not-standard-corner",
        ),
        vertex_patch: vertex_patch_state(family, classification),
        overflow_handling: option_state(
            family,
            CornerApexStrategyFamily::OverflowSmoothing,
            "support-not-overflow",
        ),
        overflow_mode: if family == CornerApexStrategyFamily::OverflowSmoothing {
            OverflowHandlingMode::SmoothOverflow
        } else {
            OverflowHandlingMode::NotApplicable
        },
        advanced_mode: advanced_mode(family, classification),
        setback: option_state(
            family,
            CornerApexStrategyFamily::Setback,
            "not-setback-case",
        ),
        deferred_construction: option_state(
            family,
            CornerApexStrategyFamily::DeferredVertexPatch,
            "not-deferred-patch",
        ),
        geometry_mode: "analytic",
        suboptimal_but_valid: family == CornerApexStrategyFamily::CliffTermination,
    }
}

fn option_state(
    family: CornerApexStrategyFamily,
    target: CornerApexStrategyFamily,
    reason: &'static str,
) -> StrategyOptionReadiness {
    if family == target {
        StrategyOptionReadiness::Selected
    } else {
        StrategyOptionReadiness::Rejected(reason)
    }
}

fn vertex_patch_state(
    family: CornerApexStrategyFamily,
    classification: VertexClassificationFact,
) -> StrategyOptionReadiness {
    if matches!(
        family,
        CornerApexStrategyFamily::VertexPatch | CornerApexStrategyFamily::HighValencePatch
    ) {
        StrategyOptionReadiness::Selected
    } else if classification.blend_types == BlendTypeCompatibility::Heterogeneous {
        StrategyOptionReadiness::ReadyAlternative
    } else {
        StrategyOptionReadiness::Rejected("not-vertex-patch-case")
    }
}

fn advanced_mode(
    family: CornerApexStrategyFamily,
    classification: VertexClassificationFact,
) -> AdvancedStrategyMode {
    if family == CornerApexStrategyFamily::Setback {
        AdvancedStrategyMode::Setback
    } else if family == CornerApexStrategyFamily::DeferredVertexPatch {
        AdvancedStrategyMode::DeferredPatch
    } else if classification.blend_types == BlendTypeCompatibility::Heterogeneous {
        AdvancedStrategyMode::MixedChamferRound
    } else {
        AdvancedStrategyMode::NotApplicable
    }
}
