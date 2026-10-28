//! Structural test: the cell complex's data-structure half depends on
//! `popkern_core` only (PK-964, plan §4.13).
//!
//! The complex starts inside `popkern-bool`. Its model, ids, operand
//! sets, validator and artifact dump must stay free of boolean-pipeline
//! types (`crate::…`), so they can move to `popkern_core` unchanged once
//! a second crate (STEP import, sewing, analysis) needs them. The
//! pipeline half (lift, radial order, regions, selection, serving) is
//! where the boolean's own types belong.
//!
//! REMEDIATION: move the code that needs a boolean type into the
//! pipeline half, or pass the data in as plain values.

use std::path::Path;

const DATA_HALF: &[&str] = &[
    "model.rs",
    "ids.rs",
    "operand_set.rs",
    "validate.rs",
    "artifact.rs",
];

#[test]
fn complex_data_half_uses_core_only() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../crates/popkern-bool/src/complex")
        .canonicalize()
        .expect("complex module exists");
    let mut violations = Vec::new();
    for name in DATA_HALF {
        let path = dir.join(name);
        let text = std::fs::read_to_string(&path).expect("readable");
        for (i, line) in text.lines().enumerate() {
            let code = line.split("//").next().unwrap_or("");
            if code.contains("crate::") {
                violations.push(format!("{name}:{}: {}", i + 1, line.trim()));
            }
        }
    }
    assert!(
        violations.is_empty(),
        "the complex data half must not use boolean-pipeline types:\n{}",
        violations.join("\n")
    );
}
