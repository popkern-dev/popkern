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
