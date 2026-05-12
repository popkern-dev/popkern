//! Lifecycle action diagnostics.
#![allow(missing_docs)]

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum LifecycleState {
    ScopedDerivedDiscardable,
    StaleGraph,
    TopologyFingerprintMismatch,
    ReusedAfterMutation,
    ArtifactNeedsCleanup(&'static str),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DiagnosticArtifactKind {
    AcceptedGraphFact,
    DiagnosticArtifact,
    DeferredExploration,
    MissingArtifactRecord,
    LeakedPlanAnnotation,
    UnknownClosure,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum LifecycleAction {
    Consumable,
    RebuildRequired(&'static str),
    RefuseStalePlan(&'static str),
    CleanupRequired(&'static str),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LifecycleFact {
    pub state: LifecycleState,
    pub artifact: DiagnosticArtifactKind,
    pub action: LifecycleAction,
}
