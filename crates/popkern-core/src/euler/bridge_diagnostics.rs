//! Read-only bridge resource diagnostics.
//!
//! This module starts as a narrow scaffold for bridge-cycle diagnostics. It
//! intentionally does not walk intervals, classify conflicts, or authorize
//! mutation yet.

use crate::error::CoreError;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

/// Local graph id for a diagnostic resource.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DiagnosticResourceId(pub usize);

/// Local graph id for a diagnostic operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DiagnosticOperationId(pub usize);

/// Local graph id for a diagnostic conflict component.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ResourceConflictId(pub usize);

/// Read-only bridge diagnostic request.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BridgeDiagnosticRequest {
    /// Planned bridge operations in caller-selected semantic order.
    pub operations: Vec<BridgeDiagnosticOperationRequest>,
}

/// One planned bridge operation to inspect.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BridgeDiagnosticOperationRequest {
    /// Semantic purpose of this operation.
    pub purpose: DiagnosticOperationPurpose,
    /// First source loop.
    pub first_loop: LoopHandle,
    /// Directed interval on the first loop, including the after-last half-edge.
    pub first_interval: Vec<HalfEdgeHandle>,
    /// Second source loop.
    pub second_loop: LoopHandle,
    /// Directed interval on the second loop, including the after-last half-edge.
    pub second_interval: Vec<HalfEdgeHandle>,
}

/// Read-only bridge diagnostic report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BridgeDiagnosticReport {
    /// Structured diagnostic graph.
    pub graph: TopologyDiagnosticGraph,
    /// Stable compact summary derived from the graph.
    pub summary: BridgeDiagnosticSummary,
}

/// Structured topology diagnostic graph.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TopologyDiagnosticGraph {
    /// Scoped topology fingerprint captured when the graph was built.
    pub version: TopologyGraphFingerprint,
    /// Diagnostic resources.
    pub resources: Vec<DiagnosticResource>,
    /// Diagnostic operations.
    pub operations: Vec<DiagnosticOperation>,
    /// Resource conflict components.
    pub conflicts: Vec<ResourceConflict>,
    /// Operator readiness entries.
    pub readiness: Vec<OperatorReadiness>,
    /// Predicted postconditions.
    pub postconditions: Vec<PostconditionPrediction>,
}

/// Scoped topology fingerprint for stale graph detection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TopologyGraphFingerprint {
    /// Vertex count at graph construction.
    pub vertices: usize,
    /// Edge count at graph construction.
    pub edges: usize,
    /// Half-edge count at graph construction.
    pub half_edges: usize,
    /// Loop count at graph construction.
    pub loops: usize,
    /// Face count at graph construction.
    pub faces: usize,
    /// Deterministic checksum over scoped loop and interval membership.
    pub scoped_membership_checksum: u64,
}

/// One diagnostic resource.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiagnosticResource {
    /// Local resource id.
    pub id: DiagnosticResourceId,
    /// Resource kind.
    pub kind: DiagnosticResourceKind,
    /// BRep handle payload for debugging and later revalidation.
    pub handle: DiagnosticHandle,
    /// Semantic/topology owner if known.
    pub owner: DiagnosticOwner,
    /// Resource availability status.
    pub status: DiagnosticResourceStatus,
}

/// One diagnostic operation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiagnosticOperation {
    /// Local operation id.
    pub id: DiagnosticOperationId,
    /// Operation kind.
    pub kind: DiagnosticOperationKind,
    /// Semantic purpose.
    pub purpose: DiagnosticOperationPurpose,
    /// Resource ids consumed by the operation.
    pub consumes: Vec<DiagnosticResourceId>,
    /// Resource kinds created by the operation.
    pub creates: Vec<DiagnosticResourceKind>,
}

/// One resource conflict component.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResourceConflict {
    /// Local conflict id.
    pub id: ResourceConflictId,
    /// Conflict kind.
    pub kind: ResourceConflictKind,
    /// Participating resource ids.
    pub resources: Vec<DiagnosticResourceId>,
    /// Participating operation ids.
    pub operations: Vec<DiagnosticOperationId>,
    /// Conflict status.
    pub status: ResourceConflictStatus,
}

/// One operator readiness entry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OperatorReadiness {
    /// Operation being evaluated.
    pub operation: DiagnosticOperationId,
    /// Candidate operator family.
    pub operator: DiagnosticOperatorFamily,
    /// Readiness status.
    pub status: OperatorReadinessStatus,
    /// Stable reason label.
    pub reason: &'static str,
}

/// Predicted postcondition placeholder.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PostconditionPrediction {
    /// Operation being predicted.
    pub operation: DiagnosticOperationId,
    /// Stable prediction label.
    pub label: &'static str,
}

/// Compact bridge diagnostic summary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BridgeDiagnosticSummary {
    /// Number of operations in the graph.
    pub operations: usize,
    /// Number of resources in the graph.
    pub resources: usize,
    /// Number of conflict components in the graph.
    pub conflicts: usize,
    /// Number of readiness entries in the graph.
    pub readiness: usize,
    /// Number of postcondition predictions in the graph.
    pub postconditions: usize,
}

/// Diagnostic resource kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DiagnosticResourceKind {
    /// Half-edge side.
    HalfEdgeSide,
    /// Edge identity.
    EdgeIdentity,
    /// Vertex.
    Vertex,
    /// Loop interval.
    LoopInterval,
    /// Face region.
    FaceRegion,
    /// Support-family region.
    SupportFamily,
    /// Cylinder segment.
    CylinderSegment,
}

/// BRep handle payload carried by a diagnostic resource.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DiagnosticHandle {
    /// No handle has been assigned yet.
    None,
    /// Half-edge handle.
    HalfEdge(HalfEdgeHandle),
    /// Loop handle.
    Loop(LoopHandle),
    /// Face handle.
    Face(FaceHandle),
}

/// Diagnostic owner payload.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DiagnosticOwner {
    /// Owner has not been classified yet.
    Unknown,
    /// Resource is owned by a loop.
    Loop(LoopHandle),
    /// Resource is owned by a face.
    Face(FaceHandle),
}

/// Resource availability status.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DiagnosticResourceStatus {
    /// Status has not been classified yet.
    Unknown,
    /// Resource is live.
    Live,
    /// Resource is stale or invalid.
    Invalid,
}

/// Diagnostic operation kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DiagnosticOperationKind {
    /// Bridge operation.
    Bridge,
}

/// Semantic purpose of a diagnostic operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DiagnosticOperationPurpose {
    /// Purpose has not been specialized yet.
    BridgeCycle,
}

/// Conflict kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ResourceConflictKind {
    /// Duplicate half-edge side consumption.
    DuplicateHalfEdgeConsumption,
    /// Duplicate boundary-edge consumption.
    DuplicateBoundaryEdgeConsumption,
    /// Shared start overlap.
    SharedStart,
    /// Shared end overlap.
    SharedEnd,
    /// Crossed endpoint overlap.
    CrossedEndpoint,
    /// Aligned overlap.
    AlignedOverlap,
    /// Nested overlap.
    NestedOverlap,
    /// Interior overlap.
    InteriorOverlap,
    /// Already-boundary overlap.
    AlreadyBoundaryOverlap,
    /// Paired side needs exposure.
    PairedSideNeedsExposure,
    /// Unsupported owner.
    UnsupportedOwner,
    /// Stale handle.
    StaleHandle,
    /// Unknown conflict.
    Unknown,
}

/// Conflict status.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ResourceConflictStatus {
    /// Conflict is valid as-is.
    Valid,
    /// Conflict can be normalized by a known operator.
    Normalizable,
    /// Conflict is valid topology but lacks an operator.
    MissingOperator,
    /// Conflict is invalid input.
    Invalid,
    /// Conflict is not classified yet.
    Unknown,
}

/// Candidate operator family.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DiagnosticOperatorFamily {
    /// Aggregate face-family bridge.
    FaceFamilyBridge,
    /// Same-loop face bridge.
    SameLoopFaceBridge,
    /// Same-loop boundary exposure.
    SameLoopBoundaryExposure,
    /// Boundary-side allocation or separator operator.
    BoundarySideAllocationOrSeparator,
}

/// Operator readiness status.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum OperatorReadinessStatus {
    /// Operator can consume the operation as-is.
    Ready,
    /// Operator can consume the operation after a known normalization.
    ReadyAfterNormalization,
    /// Required operator is missing.
    MissingOperator,
    /// Input is invalid for this operator.
    InvalidInput,
    /// Readiness has not been classified yet.
    Unknown,
}

impl BRep {
    /// Build a read-only bridge diagnostic report.
    ///
    /// Slice 1 intentionally returns an empty graph shell with a scoped
    /// fingerprint. Later slices populate inventory, conflicts, and readiness.
    pub fn bridge_diagnostic_report(
        &self,
        request: BridgeDiagnosticRequest,
    ) -> Result<BridgeDiagnosticReport, CoreError> {
        let graph =
            TopologyDiagnosticGraph::empty(TopologyGraphFingerprint::capture(self, &request));
        Ok(BridgeDiagnosticReport {
            summary: BridgeDiagnosticSummary::from_graph(&graph),
            graph,
        })
    }
}

impl TopologyDiagnosticGraph {
    fn empty(version: TopologyGraphFingerprint) -> Self {
        Self {
            version,
            resources: Vec::new(),
            operations: Vec::new(),
            conflicts: Vec::new(),
            readiness: Vec::new(),
            postconditions: Vec::new(),
        }
    }
}

impl TopologyGraphFingerprint {
    fn capture(brep: &BRep, request: &BridgeDiagnosticRequest) -> Self {
        Self {
            vertices: brep.vertices.len(),
            edges: brep.edges.len(),
            half_edges: brep.half_edges.len(),
            loops: brep.loops.len(),
            faces: brep.faces.len(),
            scoped_membership_checksum: scoped_membership_checksum(request),
        }
    }
}

impl BridgeDiagnosticSummary {
    fn from_graph(graph: &TopologyDiagnosticGraph) -> Self {
        Self {
            operations: graph.operations.len(),
            resources: graph.resources.len(),
            conflicts: graph.conflicts.len(),
            readiness: graph.readiness.len(),
            postconditions: graph.postconditions.len(),
        }
    }

    /// Stable compact diagnostic label.
    pub fn diagnostic_label(&self) -> String {
        format!(
            "bridge_diagnostics=operations={} resources={} conflicts={} readiness={} postconditions={}",
            self.operations, self.resources, self.conflicts, self.readiness, self.postconditions
        )
    }
}

fn scoped_membership_checksum(request: &BridgeDiagnosticRequest) -> u64 {
    let mut checksum = 0xcbf2_9ce4_8422_2325u64;
    for operation in &request.operations {
        checksum = checksum_debug(checksum, operation.purpose);
        checksum = checksum_debug(checksum, operation.first_loop);
        checksum = checksum_debug(checksum, &operation.first_interval);
        checksum = checksum_debug(checksum, operation.second_loop);
        checksum = checksum_debug(checksum, &operation.second_interval);
    }
    checksum
}

fn checksum_debug<T: core::fmt::Debug>(mut checksum: u64, value: T) -> u64 {
    for byte in format!("{value:?}").bytes() {
        checksum ^= u64::from(byte);
        checksum = checksum.wrapping_mul(0x0000_0100_0000_01b3);
    }
    checksum
}

#[cfg(test)]
mod tests {
    use crate::euler::{BridgeDiagnosticRequest, DiagnosticOperationPurpose};
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn empty_bridge_diagnostic_report_is_deterministic() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let request = BridgeDiagnosticRequest { operations: vec![] };

        let first = brep
            .bridge_diagnostic_report(request.clone())
            .expect("empty report builds");
        let second = brep
            .bridge_diagnostic_report(request)
            .expect("empty report builds");

        assert_eq!(first, second);
        assert_eq!(
            first.summary.diagnostic_label(),
            "bridge_diagnostics=operations=0 resources=0 conflicts=0 readiness=0 postconditions=0"
        );
    }

    #[test]
    fn bridge_diagnostic_fingerprint_tracks_scoped_request_membership() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let loop_ = brep.loops.keys().next().expect("cube has a loop");
        let half_edge = brep.loops[loop_].half_edge;
        let next = brep.half_edges[half_edge].next;

        let first = brep
            .bridge_diagnostic_report(BridgeDiagnosticRequest {
                operations: vec![crate::euler::BridgeDiagnosticOperationRequest {
                    purpose: DiagnosticOperationPurpose::BridgeCycle,
                    first_loop: loop_,
                    first_interval: vec![half_edge, next],
                    second_loop: loop_,
                    second_interval: vec![next, half_edge],
                }],
            })
            .expect("report builds");
        let second = brep
            .bridge_diagnostic_report(BridgeDiagnosticRequest {
                operations: vec![crate::euler::BridgeDiagnosticOperationRequest {
                    purpose: DiagnosticOperationPurpose::BridgeCycle,
                    first_loop: loop_,
                    first_interval: vec![next, half_edge],
                    second_loop: loop_,
                    second_interval: vec![half_edge, next],
                }],
            })
            .expect("report builds");

        assert_ne!(
            first.graph.version.scoped_membership_checksum,
            second.graph.version.scoped_membership_checksum
        );
    }
}
