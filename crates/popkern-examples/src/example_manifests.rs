use std::path::{Path, PathBuf};

/// Stable schema label for example/export run manifests.
pub const EXAMPLE_RUN_MANIFEST_SCHEMA: &str = "popkern.example_export_run";
/// Stable schema version for example/export run manifests.
pub const EXAMPLE_RUN_MANIFEST_SCHEMA_VERSION: u32 = 1;

/// Input data for writing one example/export run manifest.
pub struct ExampleRunManifestInput<'a> {
    /// Stable runner identifier, for example `pns_toms_step`.
    pub runner: &'a str,
    /// User-facing run preset or mode.
    pub preset: &'a str,
    /// Fixture or case ids selected for this run.
    pub selected_fixtures: &'a [&'a str],
    /// Per-fixture success payloads owned by the caller's domain.
    pub fixtures: Vec<serde_json::Value>,
    /// Per-fixture failure payloads owned by the caller's domain.
    pub failures: Vec<serde_json::Value>,
    /// Optional runner-specific metadata.
    pub metadata: serde_json::Value,
}

/// Write a versioned example/export run manifest.
pub fn write_example_run_manifest(
    output_dir: &Path,
    file_name: &str,
    input: ExampleRunManifestInput<'_>,
) -> Result<PathBuf, String> {
    let path = output_dir.join(file_name);
    let payload = serde_json::json!({
        "schema": EXAMPLE_RUN_MANIFEST_SCHEMA,
        "schema_version": EXAMPLE_RUN_MANIFEST_SCHEMA_VERSION,
        "runner": input.runner,
        "preset": input.preset,
        "selected_fixtures": input.selected_fixtures,
        "summary": {
            "selected": input.selected_fixtures.len(),
            "succeeded": input.fixtures.len(),
            "failed": input.failures.len(),
        },
        "metadata": input.metadata,
        "fixtures": input.fixtures,
        "failures": input.failures,
    });
    let json = serde_json::to_string_pretty(&payload).map_err(|err| err.to_string())?;
    std::fs::write(&path, json).map_err(|err| err.to_string())?;
    Ok(path)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn manifest_writer_uses_shared_schema_and_summary() {
        let path = std::env::temp_dir().join(format!(
            "popkern_example_manifest_{}.json",
            std::process::id()
        ));
        let output_dir = path.parent().expect("temp path has parent");
        let file_name = path.file_name().unwrap().to_str().unwrap();
        let selected = ["tee", "cube"];
        write_example_run_manifest(
            output_dir,
            file_name,
            ExampleRunManifestInput {
                runner: "test_runner",
                preset: "smoke",
                selected_fixtures: &selected,
                fixtures: vec![serde_json::json!({"stem": "tee"})],
                failures: vec![serde_json::json!({"stem": "cube", "reason": "blocked"})],
                metadata: serde_json::json!({"kind": "unit"}),
            },
        )
        .unwrap();
        let payload: serde_json::Value =
            serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
        let _ = std::fs::remove_file(path);
        assert_eq!(payload["schema"], EXAMPLE_RUN_MANIFEST_SCHEMA);
        assert_eq!(
            payload["schema_version"],
            EXAMPLE_RUN_MANIFEST_SCHEMA_VERSION
        );
        assert_eq!(payload["runner"], "test_runner");
        assert_eq!(payload["summary"]["selected"], 2);
        assert_eq!(payload["summary"]["succeeded"], 1);
        assert_eq!(payload["summary"]["failed"], 1);
    }
}
