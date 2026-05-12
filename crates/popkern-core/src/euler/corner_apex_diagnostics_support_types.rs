//! Support adequacy sub-facts for corner/apex diagnostics.
#![allow(missing_docs)]

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CenterLocusFeasibilityState {
    Clear,
    Penetrates,
    MissingPrimitive(&'static str),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CenterLocusFeasibilityFact {
    pub state: CenterLocusFeasibilityState,
}

impl CenterLocusFeasibilityFact {
    pub fn clear() -> Self {
        Self {
            state: CenterLocusFeasibilityState::Clear,
        }
    }

    pub fn penetrates() -> Self {
        Self {
            state: CenterLocusFeasibilityState::Penetrates,
        }
    }

    pub fn missing_primitive(reason: &'static str) -> Self {
        Self {
            state: CenterLocusFeasibilityState::MissingPrimitive(reason),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SupportFootprintFeasibilityState {
    Clear,
    Overflow,
    Cliff,
    BlendTerritoryCollision,
    BoundaryUnreachable,
    MissingPrimitive(&'static str),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SupportFootprintFeasibilityFact {
    pub overflow: SupportFootprintFeasibilityState,
    pub cliff: SupportFootprintFeasibilityState,
    pub blend_collision: SupportFootprintFeasibilityState,
    pub boundary_reachability: SupportFootprintFeasibilityState,
}

impl SupportFootprintFeasibilityFact {
    pub fn clear() -> Self {
        Self {
            overflow: SupportFootprintFeasibilityState::Clear,
            cliff: SupportFootprintFeasibilityState::Clear,
            blend_collision: SupportFootprintFeasibilityState::Clear,
            boundary_reachability: SupportFootprintFeasibilityState::Clear,
        }
    }
}
