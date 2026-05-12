//! Selection-intent facts for corner/apex diagnostics.

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CornerSelectionIntentFact {
    pub source: CornerSelectionSource,
    pub propagation: CornerSelectionPropagation,
    pub requested_setback: CornerSelectionOptionState,
    pub requested_cap: CornerSelectionOptionState,
    pub per_edge_law_assignments: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CornerSelectionSource {
    ExplicitEdgeSet,
    InferredPropagation,
    NotProvided,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CornerSelectionPropagation {
    NotRequested,
    TangentChain,
    BranchStop,
    PartialChain,
    Refused(&'static str),
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CornerSelectionOptionState {
    Requested,
    NotProvided,
}
