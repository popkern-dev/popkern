//! JSON export schema for shared BRep repair knowledge.

use serde::Serialize;

use crate::topology::{
    BrepRepairDiagnostic, BrepRepairDiagnosticSeverity, BrepRepairFact, BrepRepairFactKind,
    BrepRepairFactState, BrepRepairKnowledgeGraph, BrepRepairKnowledgeSchema,
    BrepRepairScopeExpansion, BrepRepairScopeSeed, BrepRepairSuggestedRepair,
    BrepRepairTopologyHandle, LocalToleranceEvidence, LocalToleranceSource,
};

/// Stable JSON artifact for API, GUI, and debug consumers.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct BrepRepairKnowledgeArtifact {
    /// Artifact schema version.
    pub schema_version: BrepRepairKnowledgeSchemaRecord,
    /// Requested and expanded graph scope.
    pub scope: BrepRepairKnowledgeScopeRecord,
    /// Stable facts in graph-local id order.
    pub facts: Vec<BrepRepairFactRecord>,
    /// Stable diagnostics in graph-local id order.
    pub diagnostics: Vec<BrepRepairDiagnosticRecord>,
    /// Numeric measurements extracted from facts.
    pub measurements: Vec<BrepRepairMeasurementRecord>,
    /// Pre-render-key compact labels retained only for
    /// migration/debugging.
    pub migration_labels: Vec<BrepRepairMigrationLabelRecord>,
}

/// Serializable schema version.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct BrepRepairKnowledgeSchemaRecord {
    /// Breaking schema version.
    pub major: u16,
    /// Additive schema version.
    pub minor: u16,
}

/// Serializable scope record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct BrepRepairKnowledgeScopeRecord {
    /// Initial requested seeds.
    pub seeds: Vec<String>,
    /// Expansion policy.
    pub expansion: String,
    /// Expanded handles after collection.
    pub expanded_handles: Vec<String>,
}

/// Serializable fact record.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct BrepRepairFactRecord {
    /// Local fact id.
    pub id: usize,
    /// Stable fact family label.
    pub kind: String,
    /// Fact lifecycle state.
    pub state: String,
    /// Topology handles directly described by the fact.
    pub handles: Vec<String>,
    /// Local tolerance evidence, when present.
    pub tolerance: Option<BrepRepairToleranceRecord>,
}

/// Serializable diagnostic record.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct BrepRepairDiagnosticRecord {
    /// Local diagnostic id.
    pub id: usize,
    /// Stable public diagnostic code.
    pub code: String,
    /// Severity for API and GUI consumers.
    pub severity: String,
    /// Referenced fact ids.
    pub fact_ids: Vec<usize>,
    /// Referenced topology handles.
    pub topology_handles: Vec<String>,
    /// Recorded semantic names. Empty until semantic-name facts are collected.
    pub semantic_names: Vec<String>,
    /// Suggested repair family.
    pub suggested_repair: String,
    /// Render key for GUI message catalogs.
    pub render_key: String,
}

/// Serializable measurement record.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct BrepRepairMeasurementRecord {
    /// Fact that owns this measurement.
    pub fact_id: usize,
    /// Local tolerance used by the measurement.
    pub tolerance: f64,
    /// Measured residual.
    pub residual: f64,
    /// Source of the tolerance.
    pub tolerance_source: String,
    /// Sample count when known by this schema.
    pub sample_count: Option<usize>,
}

/// Serializable tolerance record.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct BrepRepairToleranceRecord {
    /// Local tolerance.
    pub tolerance: f64,
    /// Measured residual.
    pub residual: f64,
    /// Source of the tolerance.
    pub source: String,
}

/// Serializable migration-label record (the pre-render-key compact
/// diagnostic label, retained for migration/debugging).
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct BrepRepairMigrationLabelRecord {
    /// Owning fact id.
    pub fact_id: usize,
    /// Human-readable migration label.
    pub label: String,
}

/// Build a stable JSON-ready repair knowledge artifact.
pub fn build_brep_repair_knowledge_artifact(
    graph: &BrepRepairKnowledgeGraph,
) -> BrepRepairKnowledgeArtifact {
    build_brep_repair_knowledge_artifact_with_samples(graph, None)
}

/// Build a stable JSON-ready repair knowledge artifact with known sample count.
pub fn build_brep_repair_knowledge_artifact_with_samples(
    graph: &BrepRepairKnowledgeGraph,
    sample_count: Option<usize>,
) -> BrepRepairKnowledgeArtifact {
    BrepRepairKnowledgeArtifact {
        schema_version: schema_record(graph.schema),
        scope: scope_record(graph),
        facts: graph.facts.iter().map(fact_record).collect(),
        diagnostics: graph.diagnostics.iter().map(diagnostic_record).collect(),
        measurements: graph
            .facts
            .iter()
            .filter_map(|fact| measurement_record(fact, sample_count))
            .collect(),
        migration_labels: graph
            .facts
            .iter()
            .filter_map(migration_label_record)
            .collect(),
    }
}

/// Render a diagnostic from structured fields, not from legacy labels.
pub fn render_brep_repair_diagnostic(diagnostic: &BrepRepairDiagnostic) -> String {
    let subject = diagnostic
        .handles
        .first()
        .map(|handle| handle_label(*handle))
        .unwrap_or_else(|| "unscoped BRep".to_string());
    format!(
        "{} at {} via {}",
        diagnostic.code.0,
        subject,
        suggested_repair_label(diagnostic.suggested_repair)
    )
}

fn schema_record(schema: BrepRepairKnowledgeSchema) -> BrepRepairKnowledgeSchemaRecord {
    BrepRepairKnowledgeSchemaRecord {
        major: schema.major,
        minor: schema.minor,
    }
}

fn scope_record(graph: &BrepRepairKnowledgeGraph) -> BrepRepairKnowledgeScopeRecord {
    let mut expanded = Vec::new();
    expanded.extend(
        graph
            .scope
            .expanded
            .vertices
            .iter()
            .map(|handle| handle_label(BrepRepairTopologyHandle::Vertex(*handle))),
    );
    expanded.extend(
        graph
            .scope
            .expanded
            .half_edges
            .iter()
            .map(|handle| handle_label(BrepRepairTopologyHandle::HalfEdge(*handle))),
    );
    expanded.extend(
        graph
            .scope
            .expanded
            .loops
            .iter()
            .map(|handle| handle_label(BrepRepairTopologyHandle::Loop(*handle))),
    );
    expanded.extend(
        graph
            .scope
            .expanded
            .faces
            .iter()
            .map(|handle| handle_label(BrepRepairTopologyHandle::Face(*handle))),
    );
    BrepRepairKnowledgeScopeRecord {
        seeds: graph.scope.seeds.iter().map(seed_label).collect(),
        expansion: expansion_label(graph.scope.expansion).to_string(),
        expanded_handles: expanded,
    }
}

fn fact_record(fact: &BrepRepairFact) -> BrepRepairFactRecord {
    BrepRepairFactRecord {
        id: fact.header.id.0,
        kind: fact_kind_label(fact.header.kind).to_string(),
        state: fact_state_label(fact.header.state).to_string(),
        handles: fact
            .handles
            .iter()
            .map(|handle| handle_label(*handle))
            .collect(),
        tolerance: fact.header.tolerance.map(tolerance_record),
    }
}

fn diagnostic_record(diagnostic: &BrepRepairDiagnostic) -> BrepRepairDiagnosticRecord {
    BrepRepairDiagnosticRecord {
        id: diagnostic.id.0,
        code: diagnostic.code.0.clone(),
        severity: severity_label(diagnostic.severity).to_string(),
        fact_ids: diagnostic.facts.iter().map(|id| id.0).collect(),
        topology_handles: diagnostic
            .handles
            .iter()
            .map(|handle| handle_label(*handle))
            .collect(),
        semantic_names: Vec::new(),
        suggested_repair: suggested_repair_label(diagnostic.suggested_repair).to_string(),
        render_key: diagnostic.render_key.clone(),
    }
}

fn measurement_record(
    fact: &BrepRepairFact,
    sample_count: Option<usize>,
) -> Option<BrepRepairMeasurementRecord> {
    fact.header
        .tolerance
        .map(|tolerance| BrepRepairMeasurementRecord {
            fact_id: fact.header.id.0,
            tolerance: tolerance.tolerance,
            residual: tolerance.residual,
            tolerance_source: tolerance_source_label(tolerance.source).to_string(),
            sample_count,
        })
}

fn migration_label_record(fact: &BrepRepairFact) -> Option<BrepRepairMigrationLabelRecord> {
    fact.migration_label
        .as_ref()
        .map(|label| BrepRepairMigrationLabelRecord {
            fact_id: fact.header.id.0,
            label: label.clone(),
        })
}

fn tolerance_record(tolerance: LocalToleranceEvidence) -> BrepRepairToleranceRecord {
    BrepRepairToleranceRecord {
        tolerance: tolerance.tolerance,
        residual: tolerance.residual,
        source: tolerance_source_label(tolerance.source).to_string(),
    }
}

fn seed_label(seed: &BrepRepairScopeSeed) -> String {
    match *seed {
        BrepRepairScopeSeed::Face(handle) => handle_label(BrepRepairTopologyHandle::Face(handle)),
        BrepRepairScopeSeed::Edge(handle) => handle_label(BrepRepairTopologyHandle::Edge(handle)),
        BrepRepairScopeSeed::Vertex(handle) => {
            handle_label(BrepRepairTopologyHandle::Vertex(handle))
        }
        BrepRepairScopeSeed::Loop(handle) => handle_label(BrepRepairTopologyHandle::Loop(handle)),
        BrepRepairScopeSeed::HalfEdge(handle) => {
            handle_label(BrepRepairTopologyHandle::HalfEdge(handle))
        }
    }
}

fn expansion_label(expansion: BrepRepairScopeExpansion) -> &'static str {
    match expansion {
        BrepRepairScopeExpansion::SeedsOnly => "seeds_only",
        BrepRepairScopeExpansion::IncidentLoops => "incident_loops",
        BrepRepairScopeExpansion::AdjacentFaces => "adjacent_faces",
        BrepRepairScopeExpansion::RadialCycle => "radial_cycle",
        BrepRepairScopeExpansion::Shell => "shell",
        BrepRepairScopeExpansion::FeatureOwnedSubcomplex(_) => "feature_owned_subcomplex",
    }
}

fn handle_label(handle: BrepRepairTopologyHandle) -> String {
    match handle {
        BrepRepairTopologyHandle::Vertex(handle) => format!("vertex:{handle:?}"),
        BrepRepairTopologyHandle::HalfEdge(handle) => format!("half_edge:{handle:?}"),
        BrepRepairTopologyHandle::Edge(handle) => format!("edge:{handle:?}"),
        BrepRepairTopologyHandle::Loop(handle) => format!("loop:{handle:?}"),
        BrepRepairTopologyHandle::Face(handle) => format!("face:{handle:?}"),
        BrepRepairTopologyHandle::Shell(handle) => format!("shell:{handle:?}"),
        BrepRepairTopologyHandle::Solid(handle) => format!("solid:{handle:?}"),
    }
}

fn fact_kind_label(kind: BrepRepairFactKind) -> &'static str {
    match kind {
        BrepRepairFactKind::Fingerprint => "fingerprint",
        BrepRepairFactKind::Loop => "loop",
        BrepRepairFactKind::LoopAdjacency => "loop_adjacency",
        BrepRepairFactKind::EndpointMembership => "endpoint_membership",
        BrepRepairFactKind::TrimEndpointAgreement => "trim_endpoint_agreement",
        BrepRepairFactKind::SurfaceBoundary => "surface_boundary",
        BrepRepairFactKind::EdgeSurfaceCandidate => "edge_surface_candidate",
        BrepRepairFactKind::Authority => "authority",
        BrepRepairFactKind::OperatorPolicy => "operator_policy",
    }
}

fn fact_state_label(state: BrepRepairFactState) -> &'static str {
    match state {
        BrepRepairFactState::Observed => "observed",
        BrepRepairFactState::Derived => "derived",
        BrepRepairFactState::Candidate => "candidate",
        BrepRepairFactState::Decision => "decision",
        BrepRepairFactState::Stale => "stale",
        BrepRepairFactState::Contradictory => "contradictory",
        BrepRepairFactState::Unavailable => "unavailable",
    }
}

fn severity_label(severity: BrepRepairDiagnosticSeverity) -> &'static str {
    match severity {
        BrepRepairDiagnosticSeverity::Info => "info",
        BrepRepairDiagnosticSeverity::Warning => "warning",
        BrepRepairDiagnosticSeverity::RepairableError => "repairable_error",
        BrepRepairDiagnosticSeverity::HardError => "hard_error",
        BrepRepairDiagnosticSeverity::InternalError => "internal_error",
    }
}

fn suggested_repair_label(repair: BrepRepairSuggestedRepair) -> &'static str {
    match repair {
        BrepRepairSuggestedRepair::None => "none",
        BrepRepairSuggestedRepair::Classification => "classification",
        BrepRepairSuggestedRepair::Sewing => "sewing",
        BrepRepairSuggestedRepair::TrimReconstruction => "trim_reconstruction",
        BrepRepairSuggestedRepair::EdgeCurveCompletion => "edge_curve_completion",
        BrepRepairSuggestedRepair::DegenerateEdgeNormalization => "degenerate_edge_normalization",
        BrepRepairSuggestedRepair::TopologyRebuild => "topology_rebuild",
    }
}

fn tolerance_source_label(source: LocalToleranceSource) -> &'static str {
    match source {
        LocalToleranceSource::Model => "model",
        LocalToleranceSource::Entity => "entity",
        LocalToleranceSource::RepairPolicy => "repair_policy",
        LocalToleranceSource::Import => "import",
        LocalToleranceSource::TestFixture => "test_fixture",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::{
        BRep, BrepRepairKnowledgeOptions, BrepRepairKnowledgeScope, BrepRepairScopeSeed,
        analyze_brep_repair_knowledge,
    };

    #[test]
    fn compact_json_artifact_exports_non_apex_fixture_schema() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face = brep.faces.keys().next().unwrap();
        let graph = analyze_brep_repair_knowledge(
            &brep,
            BrepRepairKnowledgeScope::seeds_only(vec![BrepRepairScopeSeed::Face(face)]),
            BrepRepairKnowledgeOptions {
                sample_budget: 5,
                ..Default::default()
            },
        );
        let artifact = build_brep_repair_knowledge_artifact(&graph);
        let json = serde_json::to_value(&artifact).unwrap();

        assert_eq!(
            artifact.schema_version.major,
            crate::topology::BREP_REPAIR_KNOWLEDGE_SCHEMA_VERSION.major
        );
        assert!(json.get("schema_version").is_some());
        assert!(json.get("scope").is_some());
        assert!(json.get("facts").is_some());
        assert!(json.get("diagnostics").is_some());
        assert!(json.get("measurements").is_some());
        assert!(json.get("migration_labels").is_some());
        assert!(
            artifact
                .diagnostics
                .iter()
                .any(|diagnostic| diagnostic.code == "brep.loop.facts")
        );
    }

    #[test]
    fn renderer_uses_structured_diagnostic_fields() {
        let diagnostic = BrepRepairDiagnostic {
            id: crate::topology::BrepRepairDiagnosticId(0),
            code: crate::topology::BrepRepairDiagnosticCode("brep.trim.pcurve_missing".to_string()),
            severity: BrepRepairDiagnosticSeverity::RepairableError,
            facts: vec![crate::topology::BrepRepairFactId(7)],
            handles: Vec::new(),
            suggested_repair: BrepRepairSuggestedRepair::TrimReconstruction,
            render_key: "brep.trim.pcurve_missing".to_string(),
            migration_label: Some("legacy text that must not be rendered".to_string()),
        };

        let rendered = render_brep_repair_diagnostic(&diagnostic);

        assert!(rendered.contains("brep.trim.pcurve_missing"));
        assert!(rendered.contains("trim_reconstruction"));
        assert!(!rendered.contains("legacy text"));
    }
}
