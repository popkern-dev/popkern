//! Construction-risk action diagnostics.
#![allow(missing_docs)]

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ConstructionRiskAction {
    None,
    Retry(&'static str),
    Reclassify(&'static str),
    Invalid(&'static str),
    MissingPrimitive(&'static str),
    RollbackRequired(&'static str),
    UnknownInputFact,
}
