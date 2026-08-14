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

/// True if `cargo_toml` declares a dependency on `dep` — as an inline
/// key (`popkern-x = { ... }`) or a table header
/// (`[dependencies.popkern-x]`, `[dev-dependencies.popkern-x]`).
///
/// A plain `contains(dep)` scan is too blunt: it fires on crate names
/// mentioned in comments (a Cargo.toml comment naming a sibling crate
/// once failed the bool gate, PK-272). Comments are stripped before
/// matching.
fn declares_dependency(cargo_toml: &str, dep: &str) -> bool {
    cargo_toml.lines().any(|line| {
        let code = line.split('#').next().unwrap_or("").trim();
        let key_form = code
            .strip_prefix(dep)
            .is_some_and(|rest| rest.trim_start().starts_with('='));
        let table_form = code.starts_with('[') && code.contains(&format!("dependencies.{dep}"));
        key_form || table_form
    })
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
            !declares_dependency(&cargo_toml, dep),
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
            !declares_dependency(&cargo_toml, dep),
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
            !declares_dependency(&cargo_toml, dep),
            "popkern-mesh must not depend on {}. \
             Move shared types to popkern-core or restructure the dependency.",
            dep
        );
    }

    // Verify popkern-core is a dependency
    assert!(
        declares_dependency(&cargo_toml, "popkern-core"),
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
            !declares_dependency(&cargo_toml, dep),
            "popkern-viz must not depend on {}. \
             Move shared types to popkern-core or restructure the dependency.",
            dep
        );
    }

    // Verify required dependencies
    assert!(
        declares_dependency(&cargo_toml, "popkern-core"),
        "popkern-viz must depend on popkern-core"
    );
    assert!(
        declares_dependency(&cargo_toml, "popkern-mesh"),
        "popkern-viz must depend on popkern-mesh"
    );
}

#[test]
fn measure_depends_only_on_core_and_mesh() {
    let cargo_toml_path = workspace_root().join("crates/popkern-measure/Cargo.toml");
    let cargo_toml =
        fs::read_to_string(&cargo_toml_path).expect("Could not read popkern-measure/Cargo.toml");

    let forbidden = [
        "popkern-history",
        "popkern-io",
        "popkern-ffi",
        "popkern-viz",
    ];

    for dep in &forbidden {
        assert!(
            !declares_dependency(&cargo_toml, dep),
            "popkern-measure must not depend on {}. \
             Move shared types to popkern-core or restructure the dependency.",
            dep
        );
    }

    assert!(
        declares_dependency(&cargo_toml, "popkern-core"),
        "popkern-measure must depend on popkern-core"
    );
    assert!(
        declares_dependency(&cargo_toml, "popkern-mesh"),
        "popkern-measure must depend on popkern-mesh"
    );
}

#[test]
fn part_is_a_top_tier_facade() {
    let cargo_toml_path = workspace_root().join("crates/popkern-part/Cargo.toml");
    let cargo_toml =
        fs::read_to_string(&cargo_toml_path).expect("Could not read popkern-part/Cargo.toml");

    // popkern-part sits at the top of the dependency graph (peer of
    // popkern-ffi): it may depend on any kernel crate but nothing may
    // sit above it except examples/ffi, and it must not pull in the
    // presentation layer.
    let forbidden = ["popkern-ffi", "popkern-viz"];
    for dep in forbidden {
        assert!(
            !declares_dependency(&cargo_toml, dep),
            "popkern-part must not depend on {dep}"
        );
    }
    assert!(
        declares_dependency(&cargo_toml, "popkern-core"),
        "popkern-part must depend on popkern-core"
    );

    // And no kernel crate may depend back on the facade.
    for below in [
        "popkern-core",
        "popkern-bool",
        "popkern-draft",
        "popkern-fillet",
        "popkern-sweep",
        "popkern-shell",
        "popkern-io",
        "popkern-mesh",
    ] {
        let toml = fs::read_to_string(workspace_root().join(format!("crates/{below}/Cargo.toml")))
            .unwrap_or_default();
        assert!(
            !declares_dependency(&toml, "popkern-part"),
            "{below} must not depend on popkern-part (facade is top-tier)"
        );
    }
}
