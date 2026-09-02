//! Audit a skin against the generated checklist: every stable op has a
//! binding whose function is registered and whose parameters cover
//! the spec's, every preview op is bound or explicitly allow-listed
//! with a reason, reserved ops are not registered, every predicate is
//! callable, and every §10 error class (with its fields) is produced
//! by the shared payload projection.

use rhai::Engine;

use crate::conformance::{Checklist, OpRow};
use crate::engine::registered_function_names;

pub use crate::conformance_bindings::{RHAI_BINDINGS, RHAI_UNBOUND_PREVIEW};
pub use crate::conformance_rust_bindings::RUST_BINDINGS;

/// How a spec op is spelled in a skin.
#[derive(Debug, Clone, Copy)]
pub struct Binding {
    /// Spec op name.
    pub op: &'static str,
    /// Registered function (Rhai) or `pub fn` (Rust) name.
    pub function: &'static str,
    /// Positional parameter names, receiver excluded.
    pub positional: &'static [&'static str],
    /// Keys of the trailing options map.
    pub options: &'static [&'static str],
    /// Spec parameters the skin fixes by convention (documented).
    pub defaulted: &'static [&'static str],
}

/// One conformance failure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Finding {
    /// The op, predicate, or error class.
    pub subject: String,
    /// What is missing.
    pub detail: String,
}

fn finding(subject: &str, detail: impl Into<String>) -> Finding {
    Finding {
        subject: subject.to_string(),
        detail: detail.into(),
    }
}

fn check_binding(row: &OpRow, binding: &Binding, registered: &[String]) -> Vec<Finding> {
    let mut out = Vec::new();
    if !registered.iter().any(|n| n == binding.function) {
        out.push(finding(
            &row.name,
            format!(
                "bound to '{}' which the engine does not register",
                binding.function
            ),
        ));
    }
    for parameter in &row.parameters {
        let covered = binding.positional.contains(&parameter.name.as_str())
            || binding.options.contains(&parameter.name.as_str())
            || binding.defaulted.contains(&parameter.name.as_str());
        if !covered {
            out.push(finding(
                &row.name,
                format!(
                    "spec parameter '{}' is neither positional, an option key, nor documented as defaulted",
                    parameter.name
                ),
            ));
        }
    }
    out
}

/// Audit the Rhai engine against the checklist.
pub fn audit_rhai(engine: &Engine, checklist: &Checklist) -> Vec<Finding> {
    let registered = registered_function_names(engine);
    let mut findings = Vec::new();
    for row in &checklist.ops {
        let binding = RHAI_BINDINGS.iter().find(|b| b.op == row.name);
        match (row.mark.as_str(), binding) {
            ("reserved", _) => {
                if registered.iter().any(|n| n == &row.name) {
                    findings.push(finding(&row.name, "reserved op must not be registered"));
                }
            }
            ("stable", Some(binding)) => findings.extend(check_binding(row, binding, &registered)),
            ("stable", None) => findings.push(finding(&row.name, "stable op has no Rhai binding")),
            ("preview", binding) => {
                let allowed = RHAI_UNBOUND_PREVIEW.iter().any(|(op, section, _)| {
                    *op == row.name && (section.is_empty() || *section == row.section)
                });
                match binding {
                    _ if allowed => {}
                    Some(binding) => findings.extend(check_binding(row, binding, &registered)),
                    None => findings.push(finding(
                        &row.name,
                        "preview op is neither bound nor allow-listed with a reason",
                    )),
                }
            }
            (other, _) => findings.push(finding(&row.name, format!("unknown mark '{other}'"))),
        }
    }
    for predicate in &checklist.predicates {
        if !registered.iter().any(|n| n == &predicate.name) {
            findings.push(finding(
                &predicate.name,
                format!("{} predicate is not registered", predicate.kind),
            ));
        }
    }
    findings
}

/// Audit the shared payload projection sources (`error_payload.rs`
/// plus the history failure types) for every §10 class and field.
pub fn audit_error_classes(checklist: &Checklist, sources: &str) -> Vec<Finding> {
    let mut findings = Vec::new();
    for class in &checklist.error_classes {
        let produced = sources.contains(&format!("\"kind\": \"{}\"", class.name))
            || sources.contains(&format!("{}(", class.name))
            || sources.contains(&format!("{} {{", class.name));
        if !produced {
            findings.push(finding(&class.name, "error class is never produced"));
            continue;
        }
        for field in &class.fields {
            let present = sources.contains(&format!("\"{field}\""))
                || sources.contains(&format!("{field}:"))
                || sources.contains(&format!("{field},"));
            if !present {
                findings.push(finding(
                    &class.name,
                    format!("payload field '{field}' does not appear in the projection"),
                ));
            }
        }
    }
    findings
}

/// Audit the Rust skin's sources: every stable op maps to a `pub fn`.
pub fn audit_rust(checklist: &Checklist, sources: &str) -> Vec<Finding> {
    let mut findings = Vec::new();
    for row in checklist.ops.iter().filter(|r| r.mark == "stable") {
        match RUST_BINDINGS.iter().find(|(op, _)| *op == row.name) {
            Some((_, needle)) => {
                if !sources.contains(&format!("pub {needle}("))
                    && !sources.contains(&format!("pub {needle}<"))
                {
                    findings.push(finding(
                        &row.name,
                        format!("Rust binding '{needle}' is not a pub fn in popkern-script"),
                    ));
                }
            }
            None => findings.push(finding(&row.name, "stable op has no Rust binding")),
        }
    }
    findings
}
