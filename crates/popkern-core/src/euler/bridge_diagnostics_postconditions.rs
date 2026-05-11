//! Bridge diagnostic postcondition preview data types.

use super::DiagnosticOperationId;

/// Local graph id for a postcondition preview.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PostconditionPredictionId(pub usize);

/// Predicted topology delta for a diagnostic construction candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DiagnosticTopologyDelta {
    /// Vertex count delta.
    pub vertices: i64,
    /// Edge count delta.
    pub edges: i64,
    /// Face count delta.
    pub faces: i64,
    /// Loop count delta.
    pub loops: i64,
}

/// Predicted postcondition for a construction candidate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PostconditionPrediction {
    /// Operation being predicted.
    pub operation: DiagnosticOperationId,
    /// Expected local topology delta.
    pub expected_delta: DiagnosticTopologyDelta,
    /// Stable prediction label.
    pub label: &'static str,
}
