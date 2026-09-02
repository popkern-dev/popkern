//! PK-413: the scripting guide cannot drift from the surface.
//!
//! - Every stable op, every selector predicate and every §10 error
//!   class in `docs/reference/scripting-surface-checklist.json` is
//!   named in `docs/reference/scripting-guide.md` by its exact
//!   spelling.
//! - Every option key the guide's feature table names for an op is
//!   a real key of that op's Rhai binding (`RHAI_BINDINGS`), and
//!   every binding key is in the guide.
//! - Every ```rhai block in the guide runs through `ScriptHost` and
//!   serves, or raises exactly the class its `// raises:` comment
//!   names.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use popkern_rhai::conformance_audit::RHAI_BINDINGS;
use popkern_rhai::{HostOptions, ScriptHost};
use serde_json::Value;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

fn guide() -> String {
    std::fs::read_to_string(root().join("docs/reference/scripting-guide.md")).expect("guide")
}

fn checklist() -> Value {
    serde_json::from_str(
        &std::fs::read_to_string(root().join("docs/reference/scripting-surface-checklist.json"))
            .expect("checklist"),
    )
    .expect("checklist json")
}

/// The op as the Rhai skin spells it (the binding's function name).
fn rhai_name(op: &str) -> String {
    RHAI_BINDINGS
        .iter()
        .find(|b| b.op == op)
        .map(|b| b.function.to_string())
        .unwrap_or_else(|| op.to_string())
}

#[test]
fn every_stable_op_predicate_and_error_class_is_named_in_the_guide() {
    let guide = guide();
    let checklist = checklist();
    let mut missing = Vec::new();
    for op in checklist["ops"].as_array().expect("ops") {
        if op["mark"] != "stable" {
            continue;
        }
        let name = rhai_name(op["name"].as_str().expect("name"));
        if !guide.contains(&format!("{name}(")) && !guide.contains(&format!("`{name}`")) {
            missing.push(format!("op {name}"));
        }
    }
    for predicate in checklist["predicates"].as_array().expect("predicates") {
        let name = predicate["name"].as_str().expect("name");
        if !guide.contains(&format!("{name}(")) {
            missing.push(format!("predicate {name}"));
        }
    }
    for class in checklist["error_classes"].as_array().expect("classes") {
        let name = class["name"].as_str().expect("name");
        if !guide.contains(&format!("`{name}`")) {
            missing.push(format!("error class {name}"));
        }
    }
    assert!(missing.is_empty(), "the guide does not name: {missing:?}");
}

#[test]
fn the_guides_option_keys_are_the_bindings_option_keys() {
    let guide = guide();
    // The feature table: one row per verb, options in the third cell.
    let mut named: BTreeSet<(String, String)> = BTreeSet::new();
    for line in guide.lines().filter(|l| l.starts_with("| `")) {
        let cells: Vec<&str> = line.split('|').collect();
        if cells.len() < 5 {
            continue;
        }
        // Every backticked verb in the first cell (`union(tool)`,
        // `subtract(tool)`, … share a row).
        let verbs: Vec<String> = backticked(cells[1])
            .into_iter()
            .map(|v| v.split('(').next().unwrap_or(&v).to_string())
            .collect();
        for key in backticked(cells[3]) {
            if key.chars().all(|c| c.is_ascii_lowercase() || c == '_') {
                for verb in &verbs {
                    named.insert((verb.clone(), key.clone()));
                }
            }
        }
    }
    let feature_verbs = [
        "extrude",
        "revolve",
        "union",
        "subtract",
        "intersect",
        "fillet",
        "chamfer",
        "shell",
        "translate",
        "rotate",
    ];
    let mut wrong = Vec::new();
    for (verb, key) in &named {
        let Some(binding) = RHAI_BINDINGS.iter().find(|b| b.function == verb) else {
            continue;
        };
        if !binding.options.contains(&key.as_str()) {
            wrong.push(format!("{verb}: `{key}` is not an option of the binding"));
        }
    }
    for verb in feature_verbs {
        let binding = RHAI_BINDINGS
            .iter()
            .find(|b| b.function == verb)
            .expect("binding");
        for key in binding.options {
            if !named.contains(&(verb.to_string(), key.to_string())) {
                wrong.push(format!(
                    "{verb}: option `{key}` is not in the guide's table"
                ));
            }
        }
    }
    assert!(wrong.is_empty(), "{wrong:?}");
}

fn backticked(cell: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut rest = cell;
    while let Some(start) = rest.find('`') {
        let after = &rest[start + 1..];
        let Some(end) = after.find('`') else { break };
        out.push(after[..end].to_string());
        rest = &after[end + 1..];
    }
    out
}

/// The ```rhai blocks of the guide with the class each must raise
/// (`None` = must serve).
fn blocks() -> Vec<(usize, String, Option<String>)> {
    let guide = guide();
    let mut out = Vec::new();
    let mut lines = guide.lines().enumerate().peekable();
    while let Some((n, line)) = lines.next() {
        if line.trim() != "```rhai" {
            continue;
        }
        let mut body = String::new();
        for (_, l) in lines.by_ref() {
            if l.trim() == "```" {
                break;
            }
            body.push_str(l);
            body.push('\n');
        }
        let raises = body
            .lines()
            .next()
            .and_then(|l| l.trim().strip_prefix("// raises:"))
            .map(|c| c.trim().to_string());
        out.push((n + 1, body, raises));
    }
    out
}

#[test]
fn every_code_block_in_the_guide_runs_as_the_prose_says() {
    let blocks = blocks();
    assert!(blocks.len() >= 10, "the guide carries runnable blocks");
    for (line, body, raises) in blocks {
        let mut host = ScriptHost::with_options(HostOptions {
            max_operations: Some(50_000_000),
            timeout: None,
        });
        let result = host.run(&body, &format!("guide:{line}"));
        match (result, raises) {
            (Ok(outcome), None) => {
                assert!(
                    outcome.model.is_some(),
                    "guide block at line {line} produced no model"
                );
            }
            (Err(error), Some(class)) => {
                assert_eq!(
                    error.kind,
                    class,
                    "guide block at line {line}: {}",
                    error.render()
                );
            }
            (Err(error), None) => {
                panic!("guide block at line {line} failed: {}", error.render())
            }
            (Ok(_), Some(class)) => {
                panic!("guide block at line {line} served but the guide says it raises {class}")
            }
        }
    }
}
