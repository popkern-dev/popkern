use std::path::Path;

use popkern_core::topology::BrepRepairActionArtifact;

const REPAIR_ACTION_SCHEMA: &str = "popkern.repair_actions";
const REPAIR_ACTION_SCHEMA_VERSION: u32 = 1;

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
