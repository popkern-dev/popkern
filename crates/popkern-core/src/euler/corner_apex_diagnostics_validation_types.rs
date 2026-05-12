//! Integration validation action diagnostics.
#![allow(missing_docs)]

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IntegrationValidationAction {
    None,
    Repair(&'static str),
    Downgrade(&'static str),
    Invalid(&'static str),
    RelaxedAcceptance(&'static str),
    SemanticCheck(&'static str),
    UnknownInputFact,
}
