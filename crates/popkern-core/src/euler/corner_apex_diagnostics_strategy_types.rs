//! Strategy metadata diagnostics for corner/apex planning.
#![allow(missing_docs)]

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StrategyOptionReadiness {
    Selected,
    ReadyAlternative,
    Rejected(&'static str),
    NotApplicable,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CapRunoutMode {
    NotApplicable,
    EndCap,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CliffTerminationMode {
    NotApplicable,
    ManagedCliff,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum OverflowHandlingMode {
    NotApplicable,
    SmoothOverflow,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AdvancedStrategyMode {
    NotApplicable,
    ChamferOnly,
    MixedChamferRound,
    Setback,
    DeferredPatch,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StrategyMetadataFact {
    pub cap_runout: StrategyOptionReadiness,
    pub cap_runout_mode: CapRunoutMode,
    pub cliff_termination: StrategyOptionReadiness,
    pub cliff_termination_mode: CliffTerminationMode,
    pub mutual_trim: StrategyOptionReadiness,
    pub three_face_blend: StrategyOptionReadiness,
    pub vertex_patch: StrategyOptionReadiness,
    pub overflow_handling: StrategyOptionReadiness,
    pub overflow_mode: OverflowHandlingMode,
    pub advanced_mode: AdvancedStrategyMode,
    pub setback: StrategyOptionReadiness,
    pub deferred_construction: StrategyOptionReadiness,
    pub geometry_mode: &'static str,
    pub suboptimal_but_valid: bool,
}
