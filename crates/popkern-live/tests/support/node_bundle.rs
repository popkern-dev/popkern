//! The nodejs wasm bundle the node smokes drive (PK-410): located,
//! optionally built (`POPKERN_LIVE_BUILD=1`), and refused when stale.

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::process::Command;

/// The workspace root.
pub fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

/// The fresh nodejs bundle directory, or `None` (with the reason
/// printed) when it is missing, stale, or node is absent.
pub fn bundle() -> Option<PathBuf> {
    let dir = root().join("target/live/node");
    let module = dir.join("popkern_live.js");
    if std::env::var("POPKERN_LIVE_BUILD").is_ok_and(|v| v == "1") {
        let status = Command::new(root().join("tools/build-live.sh"))
            .arg("--node")
            .status()
            .expect("build script");
        assert!(status.success(), "tools/build-live.sh --node failed");
    }
    if !module.exists() {
        eprintln!(
            "[pk410] SKIP: no nodejs bundle at {} (run tools/build-live.sh --node, or set \
             POPKERN_LIVE_BUILD=1)",
            module.display()
        );
        return None;
    }
    // FRESHNESS gate (PK-447): the smoke tests whatever bundle sits
    // here, and a STALE bundle reports the OLD kernel's declines as
    // current results (measured 2026-09-05: a day-old bundle declined
    // the new grind gallery script and an hour went into phantom
    // wasm-vs-native divergence). A bundle older than any workspace
    // crate source or gallery script must not speak for the kernel —
    // skip loudly, naming the newer source.
    if let Ok(bundle_meta) = std::fs::metadata(&module)
        && let Ok(bundle_mtime) = bundle_meta.modified()
    {
        let mut newest: Option<(std::time::SystemTime, PathBuf)> = None;
        let mut stack = vec![root().join("crates")];
        while let Some(d) = stack.pop() {
            let Ok(entries) = std::fs::read_dir(&d) else {
                continue;
            };
            for entry in entries.flatten() {
                let path = entry.path();
                let name = entry.file_name();
                let name = name.to_string_lossy();
                if path.is_dir() {
                    // Everything except build outputs, test-only
                    // trees and hidden dirs — those don't ship in
                    // the wasm bundle.
                    if name != "target" && name != "tests" && !name.starts_with('.') {
                        stack.push(path);
                    }
                } else if (name.ends_with(".rs") || name.ends_with(".rhai") || name == "Cargo.toml")
                    && let Ok(m) = entry.metadata()
                    && let Ok(t) = m.modified()
                    && newest.as_ref().is_none_or(|(nt, _)| t > *nt)
                {
                    newest = Some((t, path));
                }
            }
        }
        if let Some((newest_mtime, newest_path)) = newest
            && newest_mtime > bundle_mtime
        {
            eprintln!(
                "[pk410] SKIP: nodejs bundle is STALE — {} is newer than {} (rebuild with \
                 tools/build-live.sh --node, or set POPKERN_LIVE_BUILD=1); a stale \
                 bundle must not report kernel results as current",
                newest_path.display(),
                module.display()
            );
            return None;
        }
    }
    if Command::new("node").arg("--version").output().is_err() {
        eprintln!("[pk410] SKIP: node is not installed");
        return None;
    }
    Some(dir)
}
