//! Shared read-only BRep repair knowledge types.
//!
//! This module is intentionally data-only in the first migration slice. It
//! defines the stable vocabulary that apex diagnostics and generic BRep repair
//! operators will migrate toward, without collecting facts or changing repair
//! behavior.

use crate::naming::NameId;
use crate::topology::brep_repair_knowledge_authority::BrepRepairAuthorityFact;
use crate::topology::{
    BRep, BrepEndpointMembershipFacts, BrepLoopAdjacencyFacts, BrepLoopFacts,
    BrepTrimEndpointAgreementFacts, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle,
    ShellHandle, SolidHandle, VertexHandle,
};

/// Current schema version for exported BRep repair knowledge artifacts.
pub const BREP_REPAIR_KNOWLEDGE_SCHEMA_VERSION: BrepRepairKnowledgeSchema =
    BrepRepairKnowledgeSchema { major: 0, minor: 1 };

/// Versioned schema identifier for repair knowledge graphs and JSON artifacts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BrepRepairKnowledgeSchema {
    /// Breaking schema version.
    pub major: u16,
    /// Additive schema version.
    pub minor: u16,
}

/// Read-only scoped knowledge graph for BRep repair diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepRepairKnowledgeGraph {
    /// Schema version used by this graph.
    pub schema: BrepRepairKnowledgeSchema,
    /// Coarse BRep freshness fingerprint.
    pub fingerprint: BrepKnowledgeFingerprint,
    /// Scope requested by the caller.
    pub scope: BrepRepairKnowledgeScope,
    /// Loop topology facts collected in the scope.
    pub loops: Vec<BrepLoopFacts>,
    /// Loop adjacency facts collected in the scope.
    pub loop_adjacencies: Vec<BrepLoopAdjacencyFacts>,
    /// Endpoint membership facts collected in the scope.
    pub endpoint_memberships: Vec<BrepEndpointMembershipFacts>,
    /// Trim endpoint agreement facts collected in the scope.
    pub trim_endpoint_agreements: Vec<BrepTrimEndpointAgreementFacts>,
    /// Facts collected or derived inside the scope.
    pub facts: Vec<BrepRepairFact>,
    /// Diagnostics rendered from facts or scope failures.
    pub diagnostics: Vec<BrepRepairDiagnostic>,
}

/// Coarse BRep freshness fingerprint used to reject stale knowledge snapshots.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BrepKnowledgeFingerprint {
    /// Vertex count at capture.
    pub vertices: usize,
    /// Edge count at capture.
    pub edges: usize,
    /// Half-edge count at capture.
    pub half_edges: usize,
    /// Loop count at capture.
    pub loops: usize,
    /// Face count at capture.
    pub faces: usize,
    /// Shell count at capture.
    pub shells: usize,
}

impl BrepKnowledgeFingerprint {
    /// Capture a coarse freshness fingerprint from a BRep.
    pub fn from_brep(brep: &BRep) -> Self {
        Self {
            vertices: brep.vertices.len(),
            edges: brep.edges.len(),
            half_edges: brep.half_edges.len(),
            loops: brep.loops.len(),
            faces: brep.faces.len(),
            shells: brep.shells.len(),
        }
    }

    /// Compact migration-compatible diagnostic label.
    pub fn diagnostic_label(self) -> String {
        format!(
            "V={} E={} HE={} L={} F={} S={}",
            self.vertices, self.edges, self.half_edges, self.loops, self.faces, self.shells
        )
    }
}

impl BrepRepairKnowledgeGraph {
    /// Create an empty graph for additive migration and tests.
    pub fn empty(scope: BrepRepairKnowledgeScope) -> Self {
        Self {
            schema: BREP_REPAIR_KNOWLEDGE_SCHEMA_VERSION,
            fingerprint: BrepKnowledgeFingerprint {
                vertices: 0,
                edges: 0,
                half_edges: 0,
                loops: 0,
                faces: 0,
                shells: 0,
            },
            scope,
            loops: Vec::new(),
            loop_adjacencies: Vec::new(),
            endpoint_memberships: Vec::new(),
            trim_endpoint_agreements: Vec::new(),
            facts: Vec::new(),
            diagnostics: Vec::new(),
        }
    }
}

/// Caller-selected local repair knowledge scope.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BrepRepairKnowledgeScope {
    /// Initial topology seeds requested by the caller.
    pub seeds: Vec<BrepRepairScopeSeed>,
    /// Explicit expansion policy.
    pub expansion: BrepRepairScopeExpansion,
    /// Caller provenance labels attached to this scope.
    pub provenance: Vec<BrepRepairProvenanceLabel>,
    /// Final expanded handles once a collector has run.
    pub expanded: BrepRepairScopeHandles,
}

impl BrepRepairKnowledgeScope {
    /// Build a scope that only includes the supplied seeds.
    pub fn seeds_only(seeds: Vec<BrepRepairScopeSeed>) -> Self {
        Self {
            seeds,
            expansion: BrepRepairScopeExpansion::SeedsOnly,
            provenance: Vec::new(),
            expanded: BrepRepairScopeHandles::default(),
        }
    }
}

/// Topology seed for a scoped repair knowledge graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BrepRepairScopeSeed {
    /// Face seed.
    Face(FaceHandle),
    /// Edge seed.
    Edge(EdgeHandle),
    /// Vertex seed.
    Vertex(VertexHandle),
    /// Loop seed.
    Loop(LoopHandle),
    /// Half-edge seed.
    HalfEdge(HalfEdgeHandle),
}

/// Explicit expansion policy for a repair knowledge scope.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub enum BrepRepairScopeExpansion {
    /// Do not expand beyond seeds.
    #[default]
    SeedsOnly,
    /// Include loops incident to seed handles.
    IncidentLoops,
    /// Include faces adjacent to seed handles.
    AdjacentFaces,
    /// Include the radial cycle around seed edges or half-edges.
    RadialCycle,
    /// Include the full shell containing seed handles.
    Shell,
    /// Include handles owned by a semantic feature.
    FeatureOwnedSubcomplex(SemanticFeatureId),
}

/// Stable semantic-feature id used to request feature-owned scopes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SemanticFeatureId(pub u64);

/// Expanded topology handle sets for a collected scope.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BrepRepairScopeHandles {
    /// Scoped vertices.
    pub vertices: Vec<VertexHandle>,
    /// Scoped half-edges.
    pub half_edges: Vec<HalfEdgeHandle>,
    /// Scoped edges.
    pub edges: Vec<EdgeHandle>,
    /// Scoped loops.
    pub loops: Vec<LoopHandle>,
    /// Scoped faces.
    pub faces: Vec<FaceHandle>,
    /// Scoped shells.
    pub shells: Vec<ShellHandle>,
    /// Scoped solids.
    pub solids: Vec<SolidHandle>,
}

/// Caller provenance label for a graph, fact, or diagnostic.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BrepRepairProvenanceLabel(pub String);

/// Options controlling read-only knowledge collection and export.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BrepRepairKnowledgeOptions {
    /// Geometry evaluation depth.
    pub evaluation_level: BrepRepairEvaluationLevel,
    /// Maximum samples per scoped edge/trim/surface relation.
    pub sample_budget: usize,
    /// Whether preview points and parameter samples are retained.
    pub include_debug_geometry: bool,
    /// Diagnostic artifact export mode.
    pub diagnostic_export: BrepRepairDiagnosticExport,
    /// Whether whole-model expansion is permitted.
    pub whole_model_allowed: bool,
    /// Whether exported diagnostics must include all stable UX fields.
    pub strict_schema: bool,
}

impl Default for BrepRepairKnowledgeOptions {
    fn default() -> Self {
        Self {
            evaluation_level: BrepRepairEvaluationLevel::EndpointGeometry,
            sample_budget: 0,
            include_debug_geometry: false,
            diagnostic_export: BrepRepairDiagnosticExport::None,
            whole_model_allowed: false,
            strict_schema: true,
        }
    }
}

/// Geometry evaluation depth for knowledge collection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub enum BrepRepairEvaluationLevel {
    /// Only topology facts are collected.
    TopologyOnly,
    /// Endpoint geometry is evaluated.
    #[default]
    EndpointGeometry,
    /// Sampled edge/trim/surface geometry is evaluated.
    SampledGeometry,
    /// Debug geometry is retained for artifacts.
    FullDebugArtifact,
}

/// Diagnostic export mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub enum BrepRepairDiagnosticExport {
    /// Do not export a diagnostic artifact.
    #[default]
    None,
    /// Export stable compact JSON.
    CompactJson,
    /// Export full debug JSON.
    FullJson,
}

/// Local fact id within a graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BrepRepairFactId(pub usize);

/// Local diagnostic id within a graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BrepRepairDiagnosticId(pub usize);

/// One repair knowledge fact.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepRepairFact {
    /// Fact header and lifecycle metadata.
    pub header: BrepRepairFactHeader,
    /// Topology handles directly described by the fact.
    pub handles: Vec<BrepRepairTopologyHandle>,
    /// Human-readable migration label retained for compatibility.
    pub migration_label: Option<String>,
    /// Authority-specific payload when this fact classifies trusted evidence.
    pub authority: Option<BrepRepairAuthorityFact>,
}

/// Stable metadata required on every shared repair fact.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepRepairFactHeader {
    /// Local fact id.
    pub id: BrepRepairFactId,
    /// Fact family.
    pub kind: BrepRepairFactKind,
    /// Fact lifecycle state.
    pub state: BrepRepairFactState,
    /// Dependency handles and fact ids.
    pub dependencies: BrepRepairFactDependencies,
    /// Provenance labels.
    pub provenance: Vec<BrepRepairProvenanceLabel>,
    /// Local tolerance evidence, when a numeric judgment was made.
    pub tolerance: Option<LocalToleranceEvidence>,
}

/// Shared repair fact family.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BrepRepairFactKind {
    /// BRep freshness and scope fingerprint.
    Fingerprint,
    /// Loop-local topology fact.
    Loop,
    /// Loop adjacency fact.
    LoopAdjacency,
    /// Endpoint membership fact.
    EndpointMembership,
    /// Trim endpoint agreement fact.
    TrimEndpointAgreement,
    /// Surface-boundary fact.
    SurfaceBoundary,
    /// Edge/surface candidate fact.
    EdgeSurfaceCandidate,
    /// Authority classification fact.
    Authority,
    /// Operator-local policy decision fact.
    OperatorPolicy,
}

/// Lifecycle state for a repair fact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BrepRepairFactState {
    /// Read directly from topology, geometry, trims, or semantic names.
    Observed,
    /// Computed from valid facts.
    Derived,
    /// Admissible input to a repair policy.
    Candidate,
    /// Accepted or declined without mutating the BRep.
    Decision,
    /// Dependency fingerprint no longer matches the current BRep.
    Stale,
    /// Facts disagree inside the local tolerance contract.
    Contradictory,
    /// Expected evidence is missing, unevaluable, or unsupported.
    Unavailable,
}

/// Fact dependencies.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BrepRepairFactDependencies {
    /// Required topology handles.
    pub topology: Vec<BrepRepairTopologyHandle>,
    /// Required semantic names.
    pub semantic_names: Vec<NameId>,
    /// Required prior facts.
    pub facts: Vec<BrepRepairFactId>,
}

/// Topology handle carried by a shared repair fact or diagnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BrepRepairTopologyHandle {
    /// Vertex handle.
    Vertex(VertexHandle),
    /// Half-edge handle.
    HalfEdge(HalfEdgeHandle),
    /// Edge handle.
    Edge(EdgeHandle),
    /// Loop handle.
    Loop(LoopHandle),
    /// Face handle.
    Face(FaceHandle),
    /// Shell handle.
    Shell(ShellHandle),
    /// Solid handle.
    Solid(SolidHandle),
}

/// Tolerance evidence for a numeric classification.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocalToleranceEvidence {
    /// Local tolerance used by the fact.
    pub tolerance: f64,
    /// Measured residual compared against the tolerance.
    pub residual: f64,
    /// Source of the tolerance.
    pub source: LocalToleranceSource,
}

/// Source of a local tolerance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum LocalToleranceSource {
    /// BRep/model default tolerance.
    Model,
    /// Topology entity-local tolerance.
    Entity,
    /// Caller-supplied repair policy tolerance.
    RepairPolicy,
    /// Imported file or entity tolerance.
    Import,
    /// Test fixture tolerance.
    TestFixture,
}

/// Structured repair diagnostic derived from facts.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepRepairDiagnostic {
    /// Local diagnostic id.
    pub id: BrepRepairDiagnosticId,
    /// Stable public code.
    pub code: BrepRepairDiagnosticCode,
    /// Severity for API and GUI consumers.
    pub severity: BrepRepairDiagnosticSeverity,
    /// Facts summarized by this diagnostic.
    pub facts: Vec<BrepRepairFactId>,
    /// Topology handles referenced by this diagnostic.
    pub handles: Vec<BrepRepairTopologyHandle>,
    /// Suggested repair family, if applicable.
    pub suggested_repair: BrepRepairSuggestedRepair,
    /// Render key for user-facing message catalogs.
    pub render_key: String,
    /// Pre-render-key compact label retained while migration to
    /// `render_key` catalogs is in progress.
    pub migration_label: Option<String>,
}

/// Stable public diagnostic code.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BrepRepairDiagnosticCode(pub String);

/// Diagnostic severity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BrepRepairDiagnosticSeverity {
    /// Informational fact.
    Info,
    /// Warning that does not block the operation.
    Warning,
    /// Repairable error.
    RepairableError,
    /// Hard error that blocks the operation.
    HardError,
    /// Internal invariant failure.
    InternalError,
}

/// Suggested repair family for a diagnostic.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BrepRepairSuggestedRepair {
    /// No repair is suggested.
    None,
    /// Classification only.
    Classification,
    /// Topological sewing.
    Sewing,
    /// Trim reconstruction.
    TrimReconstruction,
    /// Edge-curve completion.
    EdgeCurveCompletion,
    /// Degenerate-edge normalization.
    DegenerateEdgeNormalization,
    /// Topology rebuild.
    TopologyRebuild,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_graph_uses_current_schema_and_preserves_scope() {
        let scope = BrepRepairKnowledgeScope::seeds_only(vec![]);
        let graph = BrepRepairKnowledgeGraph::empty(scope.clone());

        assert_eq!(graph.schema, BREP_REPAIR_KNOWLEDGE_SCHEMA_VERSION);
        assert_eq!(graph.scope, scope);
        assert!(graph.facts.is_empty());
        assert!(graph.diagnostics.is_empty());
    }

    #[test]
    fn fingerprint_label_matches_compact_apex_format() {
        let brep = BRep::new(crate::tolerance::Tolerance::default());
        let fingerprint = BrepKnowledgeFingerprint::from_brep(&brep);

        assert_eq!(fingerprint.diagnostic_label(), "V=0 E=0 HE=0 L=0 F=0 S=0");
    }

    #[test]
    fn default_options_are_endpoint_only_and_transition_safe() {
        let options = BrepRepairKnowledgeOptions::default();

        assert_eq!(
            options.evaluation_level,
            BrepRepairEvaluationLevel::EndpointGeometry
        );
        assert_eq!(options.sample_budget, 0);
        assert!(!options.include_debug_geometry);
        assert!(!options.whole_model_allowed);
        assert!(options.strict_schema);
    }
}
