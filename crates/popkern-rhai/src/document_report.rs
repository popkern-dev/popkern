//! Reports and outcomes over an authored document of either kind
//! (design `docs/architecture/assemblies.md` §6.5): a part document
//! reports exactly as before; an assembly reports every part
//! definition's rows in definition order, each tagged with its `part`,
//! with the summary counts (`is_clean`, `failures`,
//! `downgraded_blends`) totalled, and a `parts` table.

use popkern_assembly::Document;
use serde_json::{Value, json};

use crate::host::RunOutcome;
use crate::outcome::{Outcome, first_failure, report_value};

impl RunOutcome {
    /// Take the authored document (assembly first, then model).
    pub fn take_document(&mut self) -> Option<Document> {
        if let Some(assembly) = self.assembly.take() {
            return Some(Document::Assembly(Box::new(assembly)));
        }
        self.model.take().map(|m| Document::Part(Box::new(m)))
    }
}

/// The report value of `doc` (spec §11 shape; aggregated for an
/// assembly).
pub fn document_report_value(doc: &mut Document) -> Value {
    match doc {
        Document::Part(model) => report_value(model),
        Document::Assembly(assembly) => {
            let name = assembly.name().to_string();
            let mut features = Vec::new();
            let mut parts = Vec::new();
            let (mut failures, mut downgraded, mut clean) = (0u64, 0u64, true);
            for (part, model) in assembly.part_models_mut() {
                let part = part.unwrap_or_default();
                let report = report_value(model);
                failures += report["failures"].as_u64().unwrap_or(0);
                downgraded += report["downgraded_blends"].as_u64().unwrap_or(0);
                clean &= report["is_clean"].as_bool().unwrap_or(false);
                let rows = report["features"].as_array().cloned().unwrap_or_default();
                parts.push(json!({
                    "part": part,
                    "features": rows.len(),
                    "is_clean": report["is_clean"],
                }));
                for mut row in rows {
                    if let Value::Object(map) = &mut row {
                        map.insert("part".into(), json!(part));
                    }
                    features.push(row);
                }
            }
            json!({
                "assembly": name,
                "features": features,
                "parts": parts,
                "is_clean": clean,
                "failures": failures,
                "downgraded_blends": downgraded,
            })
        }
    }
}

/// The first failed feature of `doc` as a kernel decline (an assembly's
/// names its part).
pub fn document_first_failure(doc: &mut Document) -> Option<Outcome> {
    for (part, model) in doc.models_mut() {
        if let Some(outcome) = first_failure(model) {
            return Some(match (part, outcome) {
                (
                    Some(part),
                    Outcome::KernelDecline {
                        kind,
                        rendered,
                        mut payload,
                    },
                ) => {
                    if let Value::Object(map) = &mut payload {
                        map.insert("part".into(), json!(part));
                    }
                    Outcome::KernelDecline {
                        kind,
                        rendered: format!("in part '{part}': {rendered}"),
                        payload,
                    }
                }
                (_, outcome) => outcome,
            });
        }
    }
    None
}

/// Feature count of `doc` (summed over an assembly's parts).
pub fn document_feature_count(doc: &mut Document) -> usize {
    doc.models_mut()
        .into_iter()
        .map(|(_, m)| m.features().len())
        .sum()
}

/// A document-level failure as a run outcome: a part's recorder error
/// keeps its own class (tagged with the part); anything else is an
/// `AssemblyError` decline carrying the structured payload.
pub fn from_assembly_error(error: &popkern_assembly::AssemblyError, step: &str) -> Outcome {
    match error {
        popkern_assembly::AssemblyError::Script { source, definition } => {
            match crate::outcome::from_script_error(source, step) {
                Outcome::KernelDecline {
                    kind,
                    rendered,
                    mut payload,
                } if definition != "<document>" => {
                    if let Value::Object(map) = &mut payload {
                        map.insert("part".into(), json!(definition));
                    }
                    Outcome::KernelDecline {
                        kind,
                        rendered: format!("in part '{definition}': {rendered}"),
                        payload,
                    }
                }
                other => other,
            }
        }
        other => Outcome::KernelDecline {
            kind: "AssemblyError".into(),
            rendered: other.to_string(),
            payload: other.payload(),
        },
    }
}
