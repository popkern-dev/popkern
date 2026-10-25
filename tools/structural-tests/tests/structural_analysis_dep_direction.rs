//! Structural test: dependency direction for the analysis crates
//! (PK-948, docs/architecture/analysis-integration.md).
//!
//! Rules:
//! - `popkern-analysis` (the solver) depends on `popkern-math`,
//!   `thiserror` and `serde` only: it sits at the Math tier beside
//!   `popkern-constraint`;
//! - `popkern-analysis-brep` (the adapter) sees the solver and the
//!   core, mesh and measure crates only: it sits at the Measure tier;
//! - among the analysis crates only `popkern-analysis-cli` and
//!   `popkern-analysis-fixtures` see `popkern-io` (files and STEP);
//! - no kernel crate depends on an analysis crate except the hosts
//!   above the adapter: `popkern-live`, `popkern-server`, `popkern-cli`
//!   and `popkern-examples`.
//!
//! The check reads `[dependencies]` (normal dependencies). Dev
//! dependencies are test-only and do not cross the boundary.
//!
//! REMEDIATION: move the code to the crate the rule allows, or route
//! the data through `AnalysisModel`.

use std::path::{Path, PathBuf};

const SOLVER: &str = "popkern-analysis";
const SOLVER_ALLOWED: [&str; 3] = ["popkern-math", "thiserror", "serde"];
const ADAPTER: &str = "popkern-analysis-brep";
const ADAPTER_POPKERN_ALLOWED: [&str; 5] = [
    "popkern-analysis",
    "popkern-math",
    "popkern-core",
    "popkern-mesh",
    "popkern-measure",
];
const IO_ALLOWED: [&str; 2] = ["popkern-analysis-cli", "popkern-analysis-fixtures"];
/// The kernel crates allowed to depend on an analysis crate.
const ANALYSIS_HOSTS: [&str; 4] = [
    "popkern-live",
    "popkern-server",
    "popkern-cli",
    "popkern-examples",
];

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("Could not resolve workspace root")
}

/// Names of the normal dependencies declared in a manifest's
/// `[dependencies]` table, plus any `[dependencies.<name>]` tables.
/// Comments stripped.
fn normal_dependencies(manifest: &str) -> Vec<String> {
    let mut deps = Vec::new();
    let mut in_deps = false;
    for raw in manifest.lines() {
        let line = raw.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('[') {
            let header = line.trim_matches(|c| c == '[' || c == ']');
            in_deps = header == "dependencies";
            if let Some(name) = header.strip_prefix("dependencies.") {
                deps.push(name.to_string());
            }
            continue;
        }
        if !in_deps {
            continue;
        }
        let Some((key, _)) = line.split_once('=') else {
            continue;
        };
        let key = key.trim();
        deps.push(key.strip_suffix(".workspace").unwrap_or(key).to_string());
    }
    deps
}

fn is_analysis(krate: &str) -> bool {
    krate.starts_with("popkern-analysis")
}

/// Rule violations for crate `krate` with manifest text `manifest`.
fn violations(krate: &str, manifest: &str) -> Vec<String> {
    let mut found = Vec::new();
    for name in normal_dependencies(manifest) {
        let name = name.as_str();
        if krate == SOLVER && !SOLVER_ALLOWED.contains(&name) {
            found.push(format!(
                "{krate} depends on {name}; the solver allows only {SOLVER_ALLOWED:?}"
            ));
        }
        if krate == ADAPTER
            && name.starts_with("popkern-")
            && !ADAPTER_POPKERN_ALLOWED.contains(&name)
        {
            found.push(format!(
                "{krate} depends on {name}; the adapter allows only {ADAPTER_POPKERN_ALLOWED:?}"
            ));
        }
        if is_analysis(krate) && name == "popkern-io" && !IO_ALLOWED.contains(&krate) {
            found.push(format!(
                "{krate} depends on popkern-io; among the analysis crates only {IO_ALLOWED:?} may"
            ));
        }
        if !is_analysis(krate) && is_analysis(name) && !ANALYSIS_HOSTS.contains(&krate) {
            found.push(format!(
                "{krate} depends on {name}; only the hosts {ANALYSIS_HOSTS:?} may"
            ));
        }
    }
    found
}

/// Every crate manifest under `crates/` and `tools/`, by package
/// directory name (the analysis wasm bench is `tools/analysis-wasm-bench`,
/// package `popkern-analysis-wasm-bench`).
fn manifests() -> Vec<(String, String)> {
    let root = workspace_root();
    let mut found = Vec::new();
    for dir in ["crates", "tools"] {
        let Ok(entries) = std::fs::read_dir(root.join(dir)) else {
            continue;
        };
        for entry in entries.flatten() {
            let manifest = entry.path().join("Cargo.toml");
            let Ok(text) = std::fs::read_to_string(&manifest) else {
                continue;
            };
            let name = text
                .lines()
                .find_map(|l| l.trim().strip_prefix("name = \""))
                .and_then(|rest| rest.split('"').next())
                .unwrap_or("")
                .to_string();
            found.push((name, text));
        }
    }
    found.sort();
    found
}

#[test]
fn analysis_crates_respect_the_dependency_direction() {
    let all = manifests();
    let analysis: Vec<&str> = all
        .iter()
        .map(|(n, _)| n.as_str())
        .filter(|n| is_analysis(n))
        .collect();
    assert!(
        analysis.len() >= 5,
        "expected the five analysis crates, found {analysis:?}"
    );
    let found: Vec<String> = all
        .iter()
        .flat_map(|(name, text)| violations(name, text))
        .collect();
    assert!(
        found.is_empty(),
        "dependency direction violations:\n{}",
        found.join("\n")
    );
}

#[test]
fn the_parser_reads_every_declaration_form() {
    let manifest = "[package]\nname = \"x\"\n[dependencies]\n# popkern-io in a comment\npopkern-math.workspace = true\nthiserror = \"2\"\nserde = { workspace = true }\n[dependencies.popkern-core]\npath = \"../popkern-core\"\n[dev-dependencies]\npopkern-io = { path = \"../popkern-io\" }\n";
    assert_eq!(
        normal_dependencies(manifest),
        ["popkern-math", "thiserror", "serde", "popkern-core"]
    );
}

#[test]
fn the_solver_rule_fails_on_a_geometry_dependency() {
    let manifest = "[dependencies]\npopkern-math = { path = \"../popkern-math\" }\npopkern-core = { path = \"../popkern-core\" }\n";
    let found = violations(SOLVER, manifest);
    assert_eq!(found.len(), 1, "{found:?}");
}

#[test]
fn the_adapter_rule_fails_on_io_and_on_an_operations_crate() {
    let manifest = "[dependencies]\npopkern-io = { path = \"../popkern-io\" }\npopkern-bool = { path = \"../popkern-bool\" }\n";
    // popkern-io breaks both the adapter allow-list and the io rule.
    assert_eq!(violations(ADAPTER, manifest).len(), 3);
    assert!(violations("popkern-analysis-cli", manifest).is_empty());
}

#[test]
fn the_host_rule_fails_below_the_hosts() {
    let manifest =
        "[dependencies]\npopkern-analysis-brep = { path = \"../popkern-analysis-brep\" }\n";
    assert_eq!(violations("popkern-mesh", manifest).len(), 1);
    assert!(violations("popkern-live", manifest).is_empty());
}
