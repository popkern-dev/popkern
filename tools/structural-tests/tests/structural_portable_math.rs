//! Structural test: no platform-libm transcendentals in library code.
//!
//! `f64::sin` and friends call the platform C math library (glibc
//! natively, the Rust runtime's `libm` port in wasm), which differ in
//! the last ulp — measured under PK-674: the propeller and bored_block
//! gallery scripts built different BRep coordinates natively and in
//! the browser, and the propeller's blades tessellated differently.
//! Library code must use `popkern_math::PortableMath` (`x.psin()`,
//! `y.patan2(x)`, …), which compiles the same algorithm on every
//! target. Test code is exempt (expected values do not reach a part),
//! as are `popkern-math` itself and `popkern-examples` (native-only
//! showcase binaries and probes, never part of a wasm build).
//!
//! The scan is textual: `.sin(` on a non-float type with its own `sin`
//! (the certified interval `Iv`) is spelled `Iv::sin(x)` in library
//! code so the method form stays reserved for floats.
//!
//! REMEDIATION: `use popkern_math::PortableMath;` and call the `p`
//! method (`theta.cos()` → `theta.pcos()`).

use std::path::Path;

const TRANSCENDENTALS: &[&str] = &[
    "sin", "cos", "tan", "asin", "acos", "atan", "atan2", "sin_cos", "sinh", "cosh", "tanh",
    "asinh", "acosh", "atanh", "exp", "exp_m1", "ln", "ln_1p", "log2", "log10", "cbrt", "hypot",
    "powf",
];

fn workspace_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("Could not resolve workspace root")
}

fn is_test_file(file: &Path) -> bool {
    if file.components().any(|c| {
        c.as_os_str() == "tests" || c.as_os_str() == "benches" || c.as_os_str() == "examples"
    }) {
        return true;
    }
    let name = file
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    name == "tests.rs"
        || name.ends_with("_tests.rs")
        || name.contains("_tests_part")
        || name.starts_with("test_")
        || name.starts_with("tests_")
}

/// The transcendental a line calls through the inherent method or
/// `f64::`/`f32::` path, if any.
fn platform_call(line: &str) -> Option<&'static str> {
    TRANSCENDENTALS.iter().copied().find(|name| {
        line.contains(&format!(".{name}("))
            || line.contains(&format!("f64::{name}("))
            || line.contains(&format!("f32::{name}("))
    })
}

/// Which lines belong to a `#[cfg(test)]` item. The item runs from
/// the attribute to the brace that closes its body; a braceless item
/// (`mod tests;`, a gated `use`) covers only its own lines — it must
/// not swallow the library code after it.
fn cfg_test_lines(lines: &[&str]) -> Vec<bool> {
    let mut test_only = vec![false; lines.len()];
    let mut i = 0;
    while i < lines.len() {
        if !lines[i].trim().starts_with("#[cfg(test)]") {
            i += 1;
            continue;
        }
        let mut k = i + 1;
        while k < lines.len() && !lines[k].contains('{') && !lines[k].trim_end().ends_with(';') {
            k += 1;
        }
        if k < lines.len() && lines[k].contains('{') {
            let mut depth = 0i64;
            while k < lines.len() {
                depth += lines[k].matches('{').count() as i64;
                depth -= lines[k].matches('}').count() as i64;
                if depth <= 0 {
                    break;
                }
                k += 1;
            }
        }
        let end = k.min(lines.len() - 1);
        for flag in &mut test_only[i..=end] {
            *flag = true;
        }
        i = end + 1;
    }
    test_only
}

#[test]
fn no_platform_transcendentals_in_library_code() {
    let files = popkern_lints::find_rs_files(&workspace_root().join("crates"));
    let mut violations = Vec::new();
    for file in &files {
        let path_str = file.to_string_lossy();
        if path_str.contains("popkern-math")
            || path_str.contains("popkern-examples")
            || is_test_file(file)
        {
            continue;
        }
        let Ok(content) = std::fs::read_to_string(file) else {
            continue;
        };
        let lines: Vec<&str> = content.lines().collect();
        let test_only = cfg_test_lines(&lines);
        for (idx, line) in lines.iter().enumerate() {
            if test_only[idx] || line.trim_start().starts_with("//") {
                continue;
            }
            if let Some(name) = platform_call(line) {
                violations.push(format!("{path_str}:{}: {name}: {}", idx + 1, line.trim()));
            }
        }
    }
    assert!(
        violations.is_empty(),
        "platform-libm transcendentals in library code (results differ natively and in \
         wasm; use popkern_math::PortableMath):\n{}",
        violations.join("\n")
    );
}

#[test]
fn the_scan_sees_each_spelling() {
    assert_eq!(platform_call("let s = theta.sin();"), Some("sin"));
    assert_eq!(platform_call("let a = f64::atan2(y, x);"), Some("atan2"));
    assert_eq!(platform_call("let s = theta.psin();"), None);
    assert_eq!(platform_call("let a = y.patan2(x);"), None);
}

#[test]
fn a_braceless_test_item_does_not_hide_library_code() {
    let lines = [
        "#[cfg(test)]",
        "mod tests;",
        "fn f(x: f64) -> f64 {",
        "    x.sin()",
        "}",
        "#[cfg(test)]",
        "fn helper(",
        "    x: f64,",
        ") -> f64 {",
        "    x.cos()",
        "}",
    ];
    let test_only = cfg_test_lines(&lines);
    assert_eq!(
        test_only,
        [
            true, true, false, false, false, true, true, true, true, true, true
        ]
    );
}
