//! Selection-intent extraction for corner/apex diagnostics.

use super::selection_types::{
    CornerSelectionIntentFact, CornerSelectionOptionState, CornerSelectionPropagation,
    CornerSelectionSource,
};
use super::types::{BlendKind, CornerApexDiagnosticRequest};

pub(super) fn selection_intent(request: &CornerApexDiagnosticRequest) -> CornerSelectionIntentFact {
    let intent = request.selection_intent;
    CornerSelectionIntentFact {
        source: selection_source(intent),
        propagation: selection_propagation(intent),
        requested_setback: option_state(intent.contains("setback")),
        requested_cap: option_state(intent.contains("cap")),
        per_edge_law_assignments: request
            .selected_edges
            .iter()
            .filter(|edge| edge.blend_kind == BlendKind::FilletVariableRadius)
            .count(),
    }
}

fn selection_source(intent: &str) -> CornerSelectionSource {
    if intent.contains("explicit") || intent == "unit-test" {
        CornerSelectionSource::ExplicitEdgeSet
    } else if intent.contains("propagat") || intent.contains("inferred") {
        CornerSelectionSource::InferredPropagation
    } else if intent.is_empty() {
        CornerSelectionSource::NotProvided
    } else {
        CornerSelectionSource::Unknown
    }
}

fn selection_propagation(intent: &str) -> CornerSelectionPropagation {
    if intent.contains("branch-stop") {
        CornerSelectionPropagation::BranchStop
    } else if intent.contains("partial-chain") {
        CornerSelectionPropagation::PartialChain
    } else if intent.contains("mixed-convexity-refusal") {
        CornerSelectionPropagation::Refused("mixed-convexity")
    } else if intent.contains("tangent-chain") || intent.contains("propagat") {
        CornerSelectionPropagation::TangentChain
    } else if intent.is_empty() {
        CornerSelectionPropagation::Unknown
    } else {
        CornerSelectionPropagation::NotRequested
    }
}

fn option_state(requested: bool) -> CornerSelectionOptionState {
    if requested {
        CornerSelectionOptionState::Requested
    } else {
        CornerSelectionOptionState::NotProvided
    }
}
