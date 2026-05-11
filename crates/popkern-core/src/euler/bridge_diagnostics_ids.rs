//! Local ids for bridge diagnostic graph entries.

/// Local graph id for a diagnostic resource.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DiagnosticResourceId(pub usize);

/// Local graph id for a diagnostic operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DiagnosticOperationId(pub usize);

/// Local graph id for a diagnostic conflict component.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ResourceConflictId(pub usize);
