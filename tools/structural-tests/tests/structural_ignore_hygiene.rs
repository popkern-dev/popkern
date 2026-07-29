//! Structural test: `#[ignore]` hygiene (PK-183).
//!
//! Two contracts:
//!
//! 1. **No bare `#[ignore]`** — every ignore attribute must carry a
//!    reason string naming its class, so a reader can triage the
//!    site without reading the whole file. Acceptable classes:
//!    - `env: <why>` — needs an environment the default gate lacks
//!      (e.g. a GPU adapter); the test must skip gracefully or be
//!      unrunnable without it.
//!    - `heavy: <why + where it runs>` — cost documented and the
//!      tier that runs it named.
//!    - `known gap: PK-nnn <what>` — a real kernel gap with the
//!      owning ledger task id.
//!    Diagnostic probes are NOT an acceptable class: per the
//!    project's no-print-trace rule they are removed or promoted to
//!    structured diagnostics before commit.
//!
//! 2. **Ignore-count ratchet** — the total number of ignore
//!    attributes in the workspace may only go DOWN. The baseline
//!    lives in `ignore_count_baseline.txt`; lower it when you retire
//!    ignores. Raising it requires deliberately editing the baseline
//!    file, which is the review hook.
//!
//! REMEDIATION: fix the test and un-ignore it, delete it if it
//! asserts retired behavior, or annotate it with one of the classes
//! above (and a ledger task for `known gap`).

use std::fs;
use std::path::{Path, PathBuf};

/// Workspace root relative to this crate's manifest directory.
fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("Could not resolve workspace root")
}

/// One `#[ignore]`-attribute occurrence found in a source file.
#[derive(Debug, PartialEq)]
struct IgnoreSite {
    /// 1-based line number.
    line: usize,
    /// True when the attribute carries no reason string.
    bare: bool,
}

/// Scan Rust source text for `#[ignore]` attributes.
///
/// Line-based: a line whose first non-whitespace is a comment marker
/// is skipped, so prose mentions of `#[ignore]` in doc comments do
/// not count. `#[cfg_attr(..., ignore)]` (no reason possible in that
/// position) counts as bare.
fn scan_ignore_sites(source: &str) -> Vec<IgnoreSite> {
    let mut sites = Vec::new();
    for (idx, line) in source.lines().enumerate() {
        let trimmed = line.trim_start();
        if trimmed.starts_with("//") {
            continue;
        }
        let Some(pos) = trimmed.find("#[ignore") else {
            // cfg_attr form: `#[cfg_attr(<pred>, ignore)]`
            if trimmed.contains("#[cfg_attr") && trimmed.contains(", ignore)") {
                sites.push(IgnoreSite {
                    line: idx + 1,
                    bare: true,
                });
            }
            continue;
        };
        let rest = &trimmed[pos + "#[ignore".len()..];
        let bare = !rest.trim_start().starts_with('=');
        sites.push(IgnoreSite {
            line: idx + 1,
            bare,
        });
    }
    sites
}

/// Every `.rs` file under `crates/` (sources and tests alike).
fn workspace_rust_files() -> Vec<PathBuf> {
    let mut files = Vec::new();
    let mut stack = vec![workspace_root().join("crates")];
    while let Some(dir) = stack.pop() {
        let Ok(entries) = fs::read_dir(&dir) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                if path.file_name().is_some_and(|n| n == "target") {
                    continue;
                }
                stack.push(path);
            } else if path.extension().is_some_and(|e| e == "rs") {
                files.push(path);
            }
        }
    }
    files.sort();
    files
}

#[test]
fn no_bare_ignore_attributes() {
    let root = workspace_root();
    let mut offenders = Vec::new();
    for path in workspace_rust_files() {
        let Ok(source) = fs::read_to_string(&path) else {
            continue;
        };
        for site in scan_ignore_sites(&source) {
            if site.bare {
                offenders.push(format!(
                    "{}:{}",
                    path.strip_prefix(&root).unwrap_or(&path).display(),
                    site.line
                ));
            }
        }
    }
    assert!(
        offenders.is_empty(),
        "bare #[ignore] attributes found — annotate with a class \
         (`env: …`, `heavy: …`, `known gap: PK-nnn …`), fix and un-ignore, \
         or delete the retired test (PK-183):\n  {}",
        offenders.join("\n  ")
    );
}

#[test]
fn ignore_attribute_count_only_ratchets_down() {
    let root = workspace_root();
    let baseline_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("ignore_count_baseline.txt");
    let baseline: usize = fs::read_to_string(&baseline_path)
        .expect("read ignore_count_baseline.txt")
        .trim()
        .parse()
        .expect("ignore_count_baseline.txt must contain a single integer");

    let mut sites = Vec::new();
    for path in workspace_rust_files() {
        let Ok(source) = fs::read_to_string(&path) else {
            continue;
        };
        for site in scan_ignore_sites(&source) {
            sites.push(format!(
                "{}:{}",
                path.strip_prefix(&root).unwrap_or(&path).display(),
                site.line
            ));
        }
    }
    assert!(
        sites.len() <= baseline,
        "ignore-attribute count grew: {} > baseline {}. New ignores are not \
         accepted — fix or delete instead. If an `env:`/`heavy:`/`known gap:` \
         ignore is genuinely required, raise the baseline deliberately in \
         {} and say why in the commit. Sites:\n  {}",
        sites.len(),
        baseline,
        baseline_path.display(),
        sites.join("\n  ")
    );
    if sites.len() < baseline {
        eprintln!(
            "ignore count {} is below baseline {} — lower {} to lock in the progress",
            sites.len(),
            baseline,
            baseline_path.display()
        );
    }
}

// ── Scanner negative cases (the gate is proven, not assumed) ──

#[test]
fn scanner_flags_bare_ignore() {
    let src = "#[test]\n#[ignore]\nfn probe() {}\n";
    let sites = scan_ignore_sites(src);
    assert_eq!(
        sites,
        vec![IgnoreSite {
            line: 2,
            bare: true
        }]
    );
}

#[test]
fn scanner_accepts_annotated_ignore_and_skips_comments() {
    let src = "\
/// doc that mentions `#[ignore]` in prose
// plain comment with #[ignore]
#[test]
#[ignore = \"env: requires GPU adapter\"]
fn gpu_test() {}
#[cfg_attr(debug_assertions, ignore)]
fn cfg_gated() {}
";
    let sites = scan_ignore_sites(src);
    assert_eq!(
        sites,
        vec![
            IgnoreSite {
                line: 4,
                bare: false
            },
            IgnoreSite {
                line: 6,
                bare: true
            },
        ]
    );
}
