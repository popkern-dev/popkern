//! Shared by the run tests: repository paths, the report schema check
//! and a `.vtu` reader.
#![allow(dead_code)]

use std::path::{Path, PathBuf};

use serde_json::Value;

pub fn repo(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

pub fn schema() -> Value {
    let text = std::fs::read_to_string(repo("docs/schema/analysis-report-v1.schema.json"))
        .unwrap_or_else(|e| panic!("{e}"));
    serde_json::from_str(&text).unwrap_or_else(|e| panic!("{e}"))
}

fn type_ok(v: &Value, t: &str) -> bool {
    match t {
        "object" => v.is_object(),
        "array" => v.is_array(),
        "string" => v.is_string(),
        "number" => v.is_number(),
        "integer" => v.is_u64() || v.is_i64(),
        "boolean" => v.is_boolean(),
        "null" => v.is_null(),
        _ => false,
    }
}

/// Validates `v` against `s` for the keywords the schema uses: type,
/// required, properties, items, enum, const and local $ref.
pub fn validate(root: &Value, s: &Value, v: &Value, path: &str, errors: &mut Vec<String>) {
    if let Some(r) = s.get("$ref").and_then(Value::as_str) {
        let target = r
            .trim_start_matches("#/")
            .split('/')
            .fold(root, |node, key| &node[key]);
        return validate(root, target, v, path, errors);
    }
    if let Some(t) = s.get("type") {
        let ok = match t {
            Value::String(t) => type_ok(v, t),
            Value::Array(ts) => ts.iter().filter_map(Value::as_str).any(|t| type_ok(v, t)),
            _ => false,
        };
        if !ok {
            errors.push(format!("{path}: {v} is not {t}"));
            return;
        }
    }
    if let Some(c) = s.get("const")
        && c != v
    {
        errors.push(format!("{path}: {v} is not {c}"));
    }
    if let Some(e) = s.get("enum").and_then(Value::as_array)
        && !e.contains(v)
    {
        errors.push(format!("{path}: {v} not in {e:?}"));
    }
    if let (Some(req), true) = (s.get("required").and_then(Value::as_array), v.is_object()) {
        for k in req.iter().filter_map(Value::as_str) {
            if v.get(k).is_none() {
                errors.push(format!("{path}: missing {k}"));
            }
        }
    }
    if let (Some(props), Some(obj)) = (
        s.get("properties").and_then(Value::as_object),
        v.as_object(),
    ) {
        for (k, sub) in props {
            if let Some(x) = obj.get(k) {
                validate(root, sub, x, &format!("{path}.{k}"), errors);
            }
        }
    }
    if let (Some(items), Some(arr)) = (s.get("items"), v.as_array()) {
        for (i, x) in arr.iter().enumerate() {
            validate(root, items, x, &format!("{path}[{i}]"), errors);
        }
    }
}

pub fn assert_valid(report: &Value) {
    let s = schema();
    let mut errors = Vec::new();
    validate(&s, &s, report, "$", &mut errors);
    assert!(errors.is_empty(), "{errors:#?}");
    match report["status"].as_str() {
        Some("solved") => assert!(report.get("solve").is_some() && report.get("result").is_some()),
        Some("declined") => assert!(report.get("error").is_some()),
        other => panic!("status {other:?}"),
    }
}

/// The numbers of the `<DataArray Name="name">` in `vtu`.
pub fn data_array(vtu: &str, name: &str) -> Vec<f64> {
    let tag = format!("Name=\"{name}\"");
    let start = vtu.find(&tag).unwrap_or_else(|| panic!("no {name}"));
    let body = &vtu[start..];
    let open = body.find('>').unwrap_or_else(|| panic!("{name}")) + 1;
    let close = body
        .find("</DataArray>")
        .unwrap_or_else(|| panic!("{name}"));
    body[open..close]
        .split_whitespace()
        .map(|t| {
            t.parse::<f64>()
                .unwrap_or_else(|e| panic!("{name}: {t}: {e}"))
        })
        .collect()
}
