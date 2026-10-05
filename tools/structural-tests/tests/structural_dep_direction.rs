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
        "popkern-blend",
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

/// `popkern-script` is the recorder tier above Part and History: it
/// may depend on both facades, and no kernel crate may depend on it.
#[test]
fn script_sits_above_part_and_history() {
    let cargo_toml_path = workspace_root().join("crates/popkern-script/Cargo.toml");
    let cargo_toml =
        fs::read_to_string(&cargo_toml_path).expect("Could not read popkern-script/Cargo.toml");
    for dep in ["popkern-history", "popkern-part", "popkern-core"] {
        assert!(
            declares_dependency(&cargo_toml, dep),
            "popkern-script must depend on {dep}"
        );
    }
    assert!(
        !declares_dependency(&cargo_toml, "popkern-viz"),
        "popkern-script must not depend on popkern-viz (browser viewers are JS; native viewer is separate)"
    );
    for below in [
        "popkern-core",
        "popkern-constraint",
        "popkern-bool",
        "popkern-draft",
        "popkern-blend",
        "popkern-sweep",
        "popkern-shell",
        "popkern-io",
        "popkern-mesh",
        "popkern-sketch",
        "popkern-measure",
        "popkern-viz",
        "popkern-part",
        "popkern-history",
    ] {
        let toml = fs::read_to_string(workspace_root().join(format!("crates/{below}/Cargo.toml")))
            .unwrap_or_default();
        assert!(
            !declares_dependency(&toml, "popkern-script"),
            "{below} must not depend on popkern-script (recorder is above the facades)"
        );
    }
}

/// `popkern-assembly` is the product-structure tier above the recorder
/// (design `docs/architecture/assemblies.md` §3): a part definition is
/// a recorded `Model`, so the tier depends on `popkern-script`; nothing
/// at or below the recorder may depend on it, and it never reaches up
/// to a skin or host.
#[test]
fn assembly_sits_above_script() {
    let cargo_toml =
        fs::read_to_string(workspace_root().join("crates/popkern-assembly/Cargo.toml"))
            .expect("Could not read popkern-assembly/Cargo.toml");
    assert!(
        declares_dependency(&cargo_toml, "popkern-script"),
        "popkern-assembly must depend on popkern-script"
    );
    for above in ["popkern-rhai", "popkern-cli", "popkern-live", "popkern-viz"] {
        assert!(
            !declares_dependency(&cargo_toml, above),
            "popkern-assembly must not depend on {above}"
        );
    }
    for below in [
        "popkern-core",
        "popkern-constraint",
        "popkern-bool",
        "popkern-draft",
        "popkern-blend",
        "popkern-sweep",
        "popkern-shell",
        "popkern-io",
        "popkern-mesh",
        "popkern-sketch",
        "popkern-measure",
        "popkern-viz",
        "popkern-part",
        "popkern-history",
        "popkern-script",
    ] {
        let toml = fs::read_to_string(workspace_root().join(format!("crates/{below}/Cargo.toml")))
            .unwrap_or_default();
        assert!(
            !declares_dependency(&toml, "popkern-assembly"),
            "{below} must not depend on popkern-assembly (assemblies sit above the recorder)"
        );
    }
}

/// `popkern-rhai` is a skin over `popkern-script`: it may depend on the
/// recorder and the facades it re-exports, never on the native viewer,
/// and nothing below the recorder may depend on it.
#[test]
fn rhai_sits_above_script() {
    let cargo_toml_path = workspace_root().join("crates/popkern-rhai/Cargo.toml");
    let cargo_toml =
        fs::read_to_string(&cargo_toml_path).expect("Could not read popkern-rhai/Cargo.toml");
    assert!(
        declares_dependency(&cargo_toml, "popkern-script"),
        "popkern-rhai must depend on popkern-script"
    );
    assert!(
        !declares_dependency(&cargo_toml, "popkern-viz"),
        "popkern-rhai must not depend on popkern-viz"
    );
    for below in [
        "popkern-core",
        "popkern-constraint",
        "popkern-bool",
        "popkern-draft",
        "popkern-blend",
        "popkern-sweep",
        "popkern-shell",
        "popkern-io",
        "popkern-mesh",
        "popkern-sketch",
        "popkern-measure",
        "popkern-viz",
        "popkern-part",
        "popkern-history",
        "popkern-script",
        "popkern-assembly",
    ] {
        let toml = fs::read_to_string(workspace_root().join(format!("crates/{below}/Cargo.toml")))
            .unwrap_or_default();
        assert!(
            !declares_dependency(&toml, "popkern-rhai"),
            "{below} must not depend on popkern-rhai (skins are above the recorder)"
        );
    }
}

/// `popkern-cli` is the command over the Rhai skin: it may depend on
/// the skin, the recorder, the facades and the CPU renderer; nothing
/// depends on it.
#[test]
fn cli_sits_above_rhai() {
    let cargo_toml_path = workspace_root().join("crates/popkern-cli/Cargo.toml");
    let cargo_toml =
        fs::read_to_string(&cargo_toml_path).expect("Could not read popkern-cli/Cargo.toml");
    for dep in ["popkern-rhai", "popkern-script", "popkern-viz"] {
        assert!(
            declares_dependency(&cargo_toml, dep),
            "popkern-cli must depend on {dep}"
        );
    }
    for below in [
        "popkern-core",
        "popkern-constraint",
        "popkern-bool",
        "popkern-draft",
        "popkern-blend",
        "popkern-sweep",
        "popkern-shell",
        "popkern-io",
        "popkern-mesh",
        "popkern-sketch",
        "popkern-measure",
        "popkern-viz",
        "popkern-part",
        "popkern-history",
        "popkern-script",
        "popkern-assembly",
        "popkern-rhai",
        "popkern-live",
        "popkern-examples",
    ] {
        let toml = fs::read_to_string(workspace_root().join(format!("crates/{below}/Cargo.toml")))
            .unwrap_or_default();
        assert!(
            !declares_dependency(&toml, "popkern-cli"),
            "{below} must not depend on popkern-cli (the command is the top of the graph)"
        );
    }
}

/// `popkern-project` runs project migrations over the Rhai skin: it
/// depends on the skin and the assembly document, never on a host
/// (command, browser engine, viewer), and nothing at or below the skin
/// depends on it.
#[test]
fn project_sits_above_rhai() {
    let cargo_toml = fs::read_to_string(workspace_root().join("crates/popkern-project/Cargo.toml"))
        .expect("Could not read popkern-project/Cargo.toml");
    for dep in ["popkern-rhai", "popkern-assembly"] {
        assert!(
            declares_dependency(&cargo_toml, dep),
            "popkern-project must depend on {dep}"
        );
    }
    for above in [
        "popkern-cli",
        "popkern-live",
        "popkern-viz",
        "popkern-examples",
    ] {
        assert!(
            !declares_dependency(&cargo_toml, above),
            "popkern-project must not depend on {above}"
        );
    }
    for below in [
        "popkern-core",
        "popkern-mesh",
        "popkern-part",
        "popkern-history",
        "popkern-script",
        "popkern-assembly",
        "popkern-rhai",
    ] {
        let toml = fs::read_to_string(workspace_root().join(format!("crates/{below}/Cargo.toml")))
            .unwrap_or_default();
        assert!(
            !declares_dependency(&toml, "popkern-project"),
            "{below} must not depend on popkern-project (projects sit above the skin)"
        );
    }
}

/// `popkern-server` hosts the page and native sessions over the
/// browser engine: it depends on popkern-live and popkern-project and on
/// no viewer, and only the command depends on it.
#[test]
fn server_sits_above_live() {
    let cargo_toml = fs::read_to_string(workspace_root().join("crates/popkern-server/Cargo.toml"))
        .expect("Could not read popkern-server/Cargo.toml");
    for dep in ["popkern-live", "popkern-project"] {
        assert!(
            declares_dependency(&cargo_toml, dep),
            "popkern-server must depend on {dep}"
        );
    }
    assert!(
        !declares_dependency(&cargo_toml, "popkern-viz"),
        "popkern-server must not depend on popkern-viz"
    );
    for below in [
        "popkern-core",
        "popkern-history",
        "popkern-script",
        "popkern-assembly",
        "popkern-rhai",
        "popkern-project",
        "popkern-live",
    ] {
        let toml = fs::read_to_string(workspace_root().join(format!("crates/{below}/Cargo.toml")))
            .unwrap_or_default();
        assert!(
            !declares_dependency(&toml, "popkern-server"),
            "{below} must not depend on popkern-server (the server hosts them)"
        );
    }
}

/// `popkern-live` is the browser engine over the Rhai skin: it
/// may depend on the skin, the recorder and the mesher; never on the
/// native viewer (browsers render kernel triangles in JS) and nothing
/// below the command tier depends on it.
#[test]
fn live_sits_above_rhai() {
    let cargo_toml_path = workspace_root().join("crates/popkern-live/Cargo.toml");
    let cargo_toml =
        fs::read_to_string(&cargo_toml_path).expect("Could not read popkern-live/Cargo.toml");
    for dep in ["popkern-rhai", "popkern-script", "popkern-mesh"] {
        assert!(
            declares_dependency(&cargo_toml, dep),
            "popkern-live must depend on {dep}"
        );
    }
    assert!(
        !declares_dependency(&cargo_toml, "popkern-viz"),
        "popkern-live must not depend on popkern-viz (browsers render in JS)"
    );
    for below in [
        "popkern-core",
        "popkern-constraint",
        "popkern-bool",
        "popkern-draft",
        "popkern-blend",
        "popkern-sweep",
        "popkern-shell",
        "popkern-io",
        "popkern-mesh",
        "popkern-sketch",
        "popkern-measure",
        "popkern-viz",
        "popkern-part",
        "popkern-history",
        "popkern-script",
        "popkern-assembly",
        "popkern-rhai",
    ] {
        let toml = fs::read_to_string(workspace_root().join(format!("crates/{below}/Cargo.toml")))
            .unwrap_or_default();
        assert!(
            !declares_dependency(&toml, "popkern-live"),
            "{below} must not depend on popkern-live (hosts are the top of the graph)"
        );
    }
}
