//! Structural test: no `std::time::Instant` in library code.
//!
//! `Instant::now()` panics at runtime on `wasm32-unknown-unknown`
//! (the browser target has no ambient monotonic clock), so every
//! timed kernel path that used it was a runtime trap in the browser —
//! measured under PK-369/PK-374: the first `extrude` aborted with
//! `unreachable`. Library code must use
//! `popkern_core::time::Stopwatch`, which degrades to honestly-absent
//! timing instead of aborting. Test code, benches, and
//! `popkern-examples` (native-only binaries) are exempt, as is the
//! shim itself (`popkern-core/src/time.rs`), which owns the one
//! platform-guarded `Instant` use.
//!
//! REMEDIATION: replace `std::time::Instant::now()` with
//! `popkern_core::time::Stopwatch::start()` and consume `elapsed()`
//! as `Option<Duration>` (record absence, never a fabricated zero).

use std::path::Path;

/// Workspace root relative to this crate's manifest directory.
fn workspace_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("Could not resolve workspace root")
}

#[test]
fn no_std_instant_in_library_code() {
    let crates_dir = workspace_root().join("crates");
    let files = popkern_lints::find_rs_files(&crates_dir);

    let mut violations = Vec::new();

    for file in &files {
        let path_str = file.to_string_lossy();

        // The shim itself owns the one platform-guarded Instant use.
        if path_str.ends_with("popkern-core/src/time.rs") {
            continue;
        }
        // popkern-examples is native-only (showcase binaries and
        // perf probes); it is not part of any wasm build.
        if path_str.contains("popkern-examples") {
            continue;
        }
        // Test and bench code is exempt (native-only by definition).
        let is_test_file = file.components().any(|c| {
            c.as_os_str() == "tests"
                || c.as_os_str() == "integration_tests"
                || c.as_os_str() == "benches"
        });
        if is_test_file {
            continue;
        }
        let name = file
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        if name == "tests.rs"
            || name.ends_with("_tests.rs")
            || name.starts_with("test_")
            || name.starts_with("tests_")
        {
            continue;
        }

        let Ok(content) = std::fs::read_to_string(file) else {
            continue;
        };
        let mut in_cfg_test = false;
        let mut depth_at_cfg_test = 0usize;
        let mut brace_depth = 0usize;
        for (idx, line) in content.lines().enumerate() {
            if line.contains("#[cfg(test)]") {
                in_cfg_test = true;
                depth_at_cfg_test = brace_depth;
            }
            brace_depth += line.matches('{').count();
            brace_depth = brace_depth.saturating_sub(line.matches('}').count());
            if in_cfg_test && brace_depth <= depth_at_cfg_test && line.contains('}') {
                in_cfg_test = false;
            }
            if in_cfg_test {
                continue;
            }
            if line.trim_start().starts_with("//") {
                continue;
            }
            if line.contains("time::Instant") || line.contains("Instant::now") {
                violations.push(format!("{}:{}: {}", path_str, idx + 1, line.trim()));
            }
        }
    }

    assert!(
        violations.is_empty(),
        "std::time::Instant in library code (traps on the wasm browser target; \
         use popkern_core::time::Stopwatch):\n{}",
        violations.join("\n")
    );
}
