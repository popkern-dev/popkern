//! `popkern build` / `popkern new` (PK-660): pending migrations against
//! a saved document, the divergence and missing-migration table
//! (live-projects-and-history.md §6.2), exit codes, and outputs left
//! untouched by a failed build.

use std::path::{Path, PathBuf};
use std::process::Command;

fn popkern(args: &[&str]) -> (i32, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_popkern"))
        .args(args)
        .output()
        .expect("run popkern");
    (
        out.status.code().unwrap_or(-1),
        String::from_utf8_lossy(&out.stdout).into_owned(),
    )
}

fn copy_dir(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).expect("mkdir");
    for entry in std::fs::read_dir(from).expect("read fixture") {
        let path = entry.expect("entry").path();
        let target = to.join(path.file_name().expect("name"));
        if path.is_dir() {
            copy_dir(&path, &target);
        } else {
            std::fs::copy(&path, &target).expect("copy");
        }
    }
}

/// A fresh copy of fixture `name` in a per-test folder.
fn project(name: &str, test: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("project_build")
        .join(test);
    let _ = std::fs::remove_dir_all(&dir);
    copy_dir(
        &Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../popkern-project/tests/fixtures")
            .join(name),
        &dir,
    );
    dir
}

fn s(path: &Path) -> &str {
    path.to_str().expect("utf8 path")
}

#[test]
fn a_saved_document_applies_only_pending_migrations() {
    let dir = project("boss_plate", "pending");
    let pkn = dir.join("plate.pkn");
    let held = dir.join("rounds.hold");
    std::fs::rename(dir.join("migrations/003_rounds.rhai"), &held).expect("hold back");
    let (code, _) = popkern(&["build", s(&dir), "--pkn", s(&pkn)]);
    assert_eq!(code, 0);

    std::fs::rename(&held, dir.join("migrations/003_rounds.rhai")).expect("restore");
    let (code, out) = popkern(&["build", s(&dir), "--pkn", s(&pkn), "--json"]);
    assert_eq!(code, 0, "{out}");
    let report: serde_json::Value = serde_json::from_str(&out).expect("json");
    let ran: Vec<&str> = report["versions"]
        .as_array()
        .expect("versions")
        .iter()
        .map(|v| v["name"].as_str().expect("name"))
        .collect();
    assert_eq!(ran, ["003_rounds"], "only the pending migration ran");

    let (code, out) = popkern(&["build", s(&dir), "--pkn", s(&pkn)]);
    assert_eq!(code, 0);
    assert!(out.contains("0 pending"), "{out}");
}

#[test]
fn a_diverged_history_rebuilds_or_is_refused_under_strict() {
    let dir = project("boss_plate", "diverged");
    let pkn = dir.join("plate.pkn");
    assert_eq!(popkern(&["build", s(&dir), "--pkn", s(&pkn)]).0, 0);
    let plate = dir.join("migrations/001_plate.rhai");
    let text = std::fs::read_to_string(&plate).expect("read");
    std::fs::write(&plate, text.replace("blind(0.4)", "blind(0.5)")).expect("edit");
    let before = std::fs::read(&pkn).expect("saved");

    let (code, out) = popkern(&["build", s(&dir), "--pkn", s(&pkn), "--strict"]);
    assert_eq!(code, 6, "{out}");
    assert_eq!(std::fs::read(&pkn).expect("saved"), before, "untouched");

    let (code, out) = popkern(&["build", s(&dir), "--pkn", s(&pkn)]);
    assert_eq!(code, 0, "{out}");
    assert!(out.contains("rebuilding"), "{out}");
}

#[test]
fn a_deleted_applied_migration_is_refused_unless_rebuilding() {
    let dir = project("boss_plate", "missing");
    let pkn = dir.join("plate.pkn");
    assert_eq!(popkern(&["build", s(&dir), "--pkn", s(&pkn)]).0, 0);
    std::fs::remove_file(dir.join("migrations/003_rounds.rhai")).expect("delete");
    let (code, out) = popkern(&["build", s(&dir), "--pkn", s(&pkn)]);
    assert_eq!(code, 5, "{out}");
    assert!(out.contains("003_rounds"), "{out}");
    let (code, _) = popkern(&["build", s(&dir), "--pkn", s(&pkn), "--rebuild"]);
    assert_eq!(code, 0);
}

#[test]
fn a_failing_migration_reports_its_line_and_keeps_outputs() {
    let dir = project("labelled_rim", "failing");
    let pkn = dir.join("rim.pkn");
    assert_eq!(popkern(&["build", s(&dir), "--pkn", s(&pkn)]).0, 0);
    let before = std::fs::read(&pkn).expect("saved");
    let block = dir.join("migrations/001_block.rhai");
    let text = std::fs::read_to_string(&block).expect("read");
    let without: String = text
        .lines()
        .filter(|l| !l.contains("label_edge"))
        .map(|l| format!("{l}\n"))
        .collect();
    std::fs::write(&block, without).expect("edit");

    let (code, out) = popkern(&["build", s(&dir), "--pkn", s(&pkn), "--rebuild", "--json"]);
    assert_ne!(code, 0);
    let report: serde_json::Value = serde_json::from_str(&out).expect("json");
    assert_eq!(report["error"]["migration"], "002_round.rhai");
    assert_eq!(report["error"]["line"], 2);
    assert_eq!(std::fs::read(&pkn).expect("saved"), before, "untouched");
}

#[test]
fn a_malformed_project_exits_five() {
    let dir = project("misnamed", "misnamed");
    let (code, out) = popkern(&["build", s(&dir)]);
    assert_eq!(code, 5, "{out}");
}

#[test]
fn new_scaffolds_a_project_that_builds() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("project_build")
        .join("scaffolded");
    let _ = std::fs::remove_dir_all(&dir);
    let (code, _) = popkern(&["new", s(&dir), "--name", "Bracket"]);
    assert_eq!(code, 0);
    assert!(dir.join("popkern.toml").exists());
    assert!(dir.join("migrations/001_base.rhai").exists());
    assert!(
        std::fs::read_to_string(dir.join(".gitignore"))
            .expect("gitignore")
            .contains(".popkern/")
    );
    let (code, out) = popkern(&["build", s(&dir)]);
    assert_eq!(code, 0, "{out}");
    // A second `new` never overwrites.
    assert_ne!(popkern(&["new", s(&dir)]).0, 0);
}
