//! Structural test: dependency direction enforcement.
//!
//! Ensures no crate imports from a higher layer.
//! Layer order: core (leaf) → bool → io → mesh → ffi → viz
//!
//! REMEDIATION: Move shared types to popkern-core or restructure the dependency.

use std::fs;
use std::path::Path;

/// Workspace root relative to this crate's manifest directory.
fn workspace_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("Could not resolve workspace root")
}

#[test]
fn core_has_no_internal_crate_dependencies() {
    let cargo_toml_path = workspace_root().join("crates/popkern-core/Cargo.toml");
    let cargo_toml =
        fs::read_to_string(&cargo_toml_path).expect("Could not read popkern-core/Cargo.toml");

    let forbidden = [
        "popkern-bool",
        "popkern-io",
        "popkern-mesh",
        "popkern-ffi",
        "popkern-viz",
    ];

    for dep in &forbidden {
        assert!(
            !cargo_toml.contains(dep),
            "popkern-core must not depend on {}. \
             Move shared types to popkern-core or restructure the dependency.",
            dep
        );
    }
}

#[test]
fn bool_depends_only_on_core() {
    let cargo_toml_path = workspace_root().join("crates/popkern-bool/Cargo.toml");
    let cargo_toml =
        fs::read_to_string(&cargo_toml_path).expect("Could not read popkern-bool/Cargo.toml");

    // popkern-bool must not depend on higher-layer crates
    let forbidden = ["popkern-io", "popkern-mesh", "popkern-ffi", "popkern-viz"];

    for dep in &forbidden {
        assert!(
            !cargo_toml.contains(dep),
            "popkern-bool must not depend on {}. \
             Move shared types to popkern-core or restructure the dependency.",
            dep
        );
    }
}

#[test]
fn mesh_depends_only_on_core() {
    let cargo_toml_path = workspace_root().join("crates/popkern-mesh/Cargo.toml");
    let cargo_toml =
        fs::read_to_string(&cargo_toml_path).expect("Could not read popkern-mesh/Cargo.toml");

    // popkern-mesh must not depend on higher-layer crates in [dependencies]
    // (popkern-bool is allowed only in [dev-dependencies])
    let forbidden = ["popkern-io", "popkern-ffi", "popkern-viz"];

    for dep in &forbidden {
        assert!(
            !cargo_toml.contains(dep),
            "popkern-mesh must not depend on {}. \
             Move shared types to popkern-core or restructure the dependency.",
            dep
        );
    }

    // Verify popkern-core is a dependency
    assert!(
        cargo_toml.contains("popkern-core"),
        "popkern-mesh must depend on popkern-core"
    );
}

#[test]
fn viz_depends_only_on_core_and_mesh() {
    let cargo_toml_path = workspace_root().join("crates/popkern-viz/Cargo.toml");
    let cargo_toml =
        fs::read_to_string(&cargo_toml_path).expect("Could not read popkern-viz/Cargo.toml");

    // popkern-viz must not depend on higher-layer or sibling crates it shouldn't use
    let forbidden = ["popkern-bool", "popkern-io", "popkern-ffi"];

    for dep in &forbidden {
        assert!(
            !cargo_toml.contains(dep),
            "popkern-viz must not depend on {}. \
             Move shared types to popkern-core or restructure the dependency.",
            dep
        );
    }

    // Verify required dependencies
    assert!(
        cargo_toml.contains("popkern-core"),
        "popkern-viz must depend on popkern-core"
    );
    assert!(
        cargo_toml.contains("popkern-mesh"),
        "popkern-viz must depend on popkern-mesh"
    );
}
