use std::path::Path;

use popkern_core::topology::{
    BRep, BrepRepairActionArtifact, BrepRepairDiagnosticExport, BrepRepairEvaluationLevel,
    BrepRepairKnowledgeArtifact, BrepRepairKnowledgeOptions, BrepRepairKnowledgeScope,
    BrepRepairScopeExpansion, BrepRepairScopeSeed, analyze_brep_repair_knowledge,
    build_brep_repair_knowledge_artifact_with_samples,
};

const REPAIR_ACTION_SCHEMA: &str = "popkern.repair_actions";
const REPAIR_ACTION_SCHEMA_VERSION: u32 = 1;
const REPAIR_KNOWLEDGE_SCHEMA: &str = "popkern.brep_repair_knowledge";
const REPAIR_KNOWLEDGE_SCHEMA_VERSION: u32 = 1;

pub fn write_repair_action_artifact(
    path: &Path,
    label: &str,
    artifact: &BrepRepairActionArtifact,
) -> Result<(), String> {
    let payload = serde_json::json!({
        "schema": REPAIR_ACTION_SCHEMA,
        "schema_version": REPAIR_ACTION_SCHEMA_VERSION,
        "label": label,
        "artifact": artifact,
    });
    let json = serde_json::to_string_pretty(&payload).map_err(|err| err.to_string())?;
    std::fs::write(path, json).map_err(|err| err.to_string())
}

pub fn write_repair_knowledge_artifact(
    path: &Path,
    label: &str,
    export_mode: &str,
    artifact: &BrepRepairKnowledgeArtifact,
) -> Result<(), String> {
    let payload = serde_json::json!({
        "schema": REPAIR_KNOWLEDGE_SCHEMA,
        "schema_version": REPAIR_KNOWLEDGE_SCHEMA_VERSION,
        "label": label,
        "export_mode": export_mode,
        "artifact": artifact,
    });
    let json = serde_json::to_string_pretty(&payload).map_err(|err| err.to_string())?;
    std::fs::write(path, json).map_err(|err| err.to_string())
}

pub fn write_whole_brep_repair_knowledge_artifact(
    path: &Path,
    label: &str,
    brep: &BRep,
    samples: usize,
    export: BrepRepairDiagnosticExport,
) -> Result<(), String> {
    let graph = analyze_brep_repair_knowledge(
        brep,
        whole_brep_repair_knowledge_scope(brep),
        BrepRepairKnowledgeOptions {
            evaluation_level: match export {
                BrepRepairDiagnosticExport::FullJson => {
                    BrepRepairEvaluationLevel::FullDebugArtifact
                }
                _ => BrepRepairEvaluationLevel::SampledGeometry,
            },
            sample_budget: samples,
            include_debug_geometry: matches!(export, BrepRepairDiagnosticExport::FullJson),
            diagnostic_export: export,
            whole_model_allowed: true,
            strict_schema: true,
        },
    );
    let artifact = build_brep_repair_knowledge_artifact_with_samples(&graph, Some(samples));
    write_repair_knowledge_artifact(path, label, export_mode_label(export), &artifact)
}

fn whole_brep_repair_knowledge_scope(brep: &BRep) -> BrepRepairKnowledgeScope {
    let mut scope = BrepRepairKnowledgeScope::seeds_only(
        brep.faces.keys().map(BrepRepairScopeSeed::Face).collect(),
    );
    scope.expansion = BrepRepairScopeExpansion::AdjacentFaces;
    scope
}

fn export_mode_label(export: BrepRepairDiagnosticExport) -> &'static str {
    match export {
        BrepRepairDiagnosticExport::None => "none",
        BrepRepairDiagnosticExport::CompactJson => "compact",
        BrepRepairDiagnosticExport::FullJson => "full",
    }
}
