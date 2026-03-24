//! Model file loading and saving.
//!
//! Currently supports BRep JSON (`.json`). OBJ, STL, and STEP formats
//! are recognized but return stub errors until their respective crates
//! are integrated.

use std::path::Path;

use popkern_core::topology::BRep;

use crate::error::VizError;

/// Load a model file into a [`BRep`].
///
/// Dispatches on the file extension:
/// - `.json` — BRep JSON (via serde)
/// - `.obj` — stub (not yet implemented)
/// - `.stl` — stub (not yet implemented)
/// - `.step` / `.stp` — stub (requires popkern-io)
///
/// Returns an error for unrecognized extensions.
pub fn load_model(path: &Path) -> Result<BRep, VizError> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    match ext.as_str() {
        "json" => load_brep_json(path),
        "obj" => Err(VizError::Load("OBJ loading not yet implemented".into())),
        "stl" => Err(VizError::Load("STL loading not yet implemented".into())),
        "step" | "stp" => Err(VizError::Load(
            "STEP loading not yet implemented (requires popkern-io)".into(),
        )),
        _ => Err(VizError::Load(format!("unknown file extension: {ext}"))),
    }
}

/// Load a BRep from a JSON file.
fn load_brep_json(path: &Path) -> Result<BRep, VizError> {
    let data = std::fs::read_to_string(path)?;
    let brep: BRep = serde_json::from_str(&data)?;
    Ok(brep)
}

/// Save a [`BRep`] to a JSON file.
///
/// Serializes the model with pretty-printing for human readability.
pub fn save_brep_json(brep: &BRep, path: &Path) -> Result<(), VizError> {
    let data = serde_json::to_string_pretty(brep)?;
    std::fs::write(path, data)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    #[test]
    fn load_unsupported_extension_is_error() {
        let result = load_model(Path::new("/tmp/model.xyz"));
        assert!(result.is_err());
        let msg = result.unwrap_err().to_string();
        assert!(msg.contains("unknown file extension"));
    }

    #[test]
    fn load_obj_is_stub_error() {
        let result = load_model(Path::new("/tmp/model.obj"));
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("not yet implemented")
        );
    }

    #[test]
    fn load_stl_is_stub_error() {
        let result = load_model(Path::new("/tmp/model.stl"));
        assert!(result.is_err());
    }

    #[test]
    fn load_step_is_stub_error() {
        let result = load_model(Path::new("/tmp/model.step"));
        assert!(result.is_err());
        let result = load_model(Path::new("/tmp/model.stp"));
        assert!(result.is_err());
    }

    #[test]
    fn load_missing_json_file_is_error() {
        let result = load_model(Path::new("/tmp/nonexistent_popkern_test.json"));
        assert!(result.is_err());
    }

    #[test]
    fn save_and_load_roundtrip() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");

        let dir = std::env::temp_dir();
        let path = dir.join("popkern_loader_test.json");

        save_brep_json(&brep, &path).expect("save should succeed");

        let loaded = load_brep_json(&path).expect("load should succeed");
        assert_eq!(loaded.faces.len(), brep.faces.len());
        assert_eq!(loaded.edges.len(), brep.edges.len());
        assert_eq!(loaded.vertices.len(), brep.vertices.len());

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn load_model_json_roundtrip() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let dir = std::env::temp_dir();
        let path = dir.join("popkern_loader_model_test.json");

        save_brep_json(&brep, &path).expect("save should succeed");

        let loaded = load_model(&path).expect("load_model should succeed");
        assert_eq!(loaded.faces.len(), brep.faces.len());

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn save_brep_json_creates_file() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let dir = std::env::temp_dir();
        let path = dir.join("popkern_save_test.json");

        save_brep_json(&brep, &path).expect("save should succeed");
        assert!(path.exists());
        let meta = std::fs::metadata(&path).expect("metadata");
        assert!(meta.len() > 0);

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn load_invalid_json_is_error() {
        let dir = std::env::temp_dir();
        let path = dir.join("popkern_invalid_json_test.json");

        std::fs::write(&path, "this is not valid json").expect("write");
        let result = load_model(&path);
        assert!(result.is_err());

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn extension_case_insensitive() {
        // .JSON should work the same as .json (but file may not exist)
        let result = load_model(Path::new("/tmp/nonexistent.JSON"));
        // Should fail with IO error (file not found), not "unknown extension"
        let msg = result.unwrap_err().to_string();
        assert!(!msg.contains("unknown file extension"), "got: {msg}");
    }

    #[test]
    fn empty_extension_is_error() {
        let result = load_model(Path::new("/tmp/noext"));
        assert!(result.is_err());
        assert!(
            result
                .unwrap_err()
                .to_string()
                .contains("unknown file extension")
        );
    }
}
