//! PK-412: the examples gallery cannot rot. One table drives three
//! things — the native gate on every script (served, no failures, no
//! downgraded blends, watertight display mesh within the angular facet
//! bound, STEP through the OCC oracle when present), the manifest the popkern live page reads
//! (`web/live/examples.json` + the script copies beside it),
//! and the thumbnails. `POPKERN_UPDATE_GALLERY=1` regenerates the
//! manifest, copies and thumbnails; without it the test fails on any
//! drift, so a script edit that changes a feature count is a diff in
//! review.

use std::path::{Path, PathBuf};

use popkern_cli::{Outcome, execute, parse_args};
use popkern_live::{HostOptions, Session};
use serde_json::{Value, json};

#[path = "gallery/table.rs"]
mod table;
use table::GALLERY;

/// Project examples (live-projects-and-history.md §12): folder name in
/// crates/popkern-project/gallery, title, blurb, expected feature count,
/// and what it demonstrates. The page lists them beside the scripts and
/// builds them with `run_project`.
const GALLERY_PROJECTS: &[(&str, &str, &str, usize, &str)] = &[(
    "boss_plate_in_steps",
    "Boss plate in three steps",
    "The boss plate as a project: plate, boss, rounds — each migration starts from the document the earlier ones built.",
    6,
    "projects, previous(), migrations",
)];

fn project_gallery_dir() -> PathBuf {
    root().join("crates/popkern-project/gallery")
}

/// The files of project folder `dir`, project-relative, sorted.
fn project_files(dir: &Path) -> Vec<String> {
    let mut out = vec!["popkern.toml".to_string()];
    let mut migrations: Vec<String> = std::fs::read_dir(dir.join("migrations"))
        .expect("migrations")
        .filter_map(|e| {
            e.ok()?
                .file_name()
                .to_str()
                .map(|n| format!("migrations/{n}"))
        })
        .collect();
    migrations.sort();
    out.extend(migrations);
    out
}

/// Build every project example; its manifest entry.
fn project_entries() -> Vec<Value> {
    let mut entries = Vec::new();
    for (dir, title, blurb, features, demonstrates) in GALLERY_PROJECTS {
        let folder = project_gallery_dir().join(dir);
        let files: Vec<(String, Vec<u8>)> = project_files(&folder)
            .into_iter()
            .map(|p| {
                let bytes = std::fs::read(folder.join(&p)).expect("project file");
                (p, bytes)
            })
            .collect();
        let mut session = Session::new(HostOptions::default());
        let envelope = session.run_project(files.clone());
        assert!(
            envelope.outcome.is_served(),
            "[{dir}] {:?}",
            envelope.outcome
        );
        assert_eq!(envelope.features, *features, "[{dir}] feature count");
        assert!(
            envelope.bodies.iter().all(|b| b.open_edges == 0),
            "[{dir}] watertight"
        );
        entries.push(json!({
            "project": dir,
            "files": files.iter().map(|(p, _)| p).collect::<Vec<_>>(),
            "title": title,
            "blurb": blurb,
            "features": features,
            "demonstrates": demonstrates,
        }));
    }
    entries
}

/// f32 buffer rounding on unit normals: ~1e-7 rad per component.
const FACET_SPREAD_SLACK: f64 = 1.0e-4;

/// Worst per-triangle vertex-normal spread (radians) over the page's
/// display buffers: the angular facet bound the viewer relies on for
/// smooth shading and straight zebra stripes. Pole vertices count —
/// a flipped pole normal reads as a ~180° spread.
fn worst_facet_spread(normals: &[f32], indices: &[u32]) -> f64 {
    let normal = |i: u32| {
        let k = 3 * i as usize;
        [
            f64::from(normals[k]),
            f64::from(normals[k + 1]),
            f64::from(normals[k + 2]),
        ]
    };
    let angle = |a: [f64; 3], b: [f64; 3]| {
        let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let la = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
        let lb = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
        (dot / (la * lb)).clamp(-1.0, 1.0).acos()
    };
    indices
        .chunks(3)
        .map(|t| {
            let n = [normal(t[0]), normal(t[1]), normal(t[2])];
            angle(n[0], n[1])
                .max(angle(n[1], n[2]))
                .max(angle(n[2], n[0]))
        })
        .fold(0.0, f64::max)
}

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

fn examples_dir() -> PathBuf {
    root().join("crates/popkern-rhai/examples")
}

fn page_dir() -> PathBuf {
    root().join("web/live")
}

fn updating() -> bool {
    std::env::var("POPKERN_UPDATE_GALLERY").is_ok_and(|v| v == "1")
}

fn tmp(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("gallery");
    std::fs::create_dir_all(&dir).expect("dir");
    dir.join(name)
}

fn oracle_passes(step: &Path) -> Option<bool> {
    let oracle = Path::new(&std::env::var("HOME").unwrap_or_default()).join("venvs/ocp/bin/python");
    let check = root().join("tools/occ-oracle/occ_check.py");
    if !(oracle.exists() && check.exists()) {
        return None;
    }
    let output = std::process::Command::new(&oracle)
        .arg(&check)
        .arg(step)
        .output()
        .expect("oracle");
    Some(String::from_utf8_lossy(&output.stdout).contains("\"passed\": true"))
}

#[test]
fn every_gallery_script_serves_clean_watertight_and_step_valid() {
    assert!(GALLERY.len() >= 8, "at least eight gallery scripts");
    let mut entries = Vec::new();
    for (stem, title, blurb, features, demonstrates) in GALLERY {
        let script = examples_dir().join(format!("{stem}.rhai"));
        let source = std::fs::read_to_string(&script).unwrap_or_else(|e| panic!("{stem}: {e}"));
        // Native engine: served, clean, every body watertight.
        let mut session = Session::new(HostOptions::default());
        let envelope = session.run(&source, &format!("{stem}.rhai"));
        assert!(
            envelope.outcome.is_served(),
            "[{stem}] {:?}",
            envelope.outcome
        );
        assert_eq!(envelope.features, *features, "[{stem}] feature count");
        let report = envelope.report.as_ref().expect("report");
        assert_eq!(report["failures"], 0, "[{stem}]");
        assert_eq!(
            report["downgraded_blends"], 0,
            "[{stem}] no script leans on a downgrade"
        );
        assert!(!envelope.bodies.is_empty(), "[{stem}] a body on screen");
        for body in &envelope.bodies {
            assert_eq!(body.open_edges, 0, "[{stem}] watertight mesh");
            let worst = worst_facet_spread(&body.normals, &body.indices);
            assert!(
                worst <= popkern_mesh::DISPLAY_MAX_FACET_ANGLE + FACET_SPREAD_SLACK,
                "[{stem}] display mesh facet spread {:.2}° over the {:.1}° bound",
                worst.to_degrees(),
                popkern_mesh::DISPLAY_MAX_FACET_ANGLE.to_degrees()
            );
        }
        // The command: STEP + thumbnail, oracle on the STEP.
        let step = tmp(&format!("{stem}.step"));
        let png = if updating() {
            page_dir().join("examples").join(format!("{stem}.png"))
        } else {
            tmp(&format!("{stem}.png"))
        };
        let args: Vec<String> = [
            "run",
            script.to_str().expect("utf8"),
            "--step",
            step.to_str().expect("utf8"),
            "--png",
            png.to_str().expect("utf8"),
            "--width",
            "480",
            "--height",
            "360",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let summary = execute(&parse_args(&args).expect("args"));
        assert!(
            matches!(summary.outcome, Outcome::Served),
            "[{stem}] {:?}",
            summary.outcome
        );
        match oracle_passes(&step) {
            Some(passed) => assert!(passed, "[{stem}] OCC oracle rejected the STEP"),
            None => eprintln!("[pk412] {stem}: OCC oracle absent"),
        }
        entries.push(json!({
            "file": format!("{stem}.rhai"),
            "title": title,
            "blurb": blurb,
            "features": features,
            "demonstrates": demonstrates,
            "thumbnail": format!("examples/{stem}.png"),
        }));
    }
    entries.extend(project_entries());
    // The manifest and the script copies the page serves.
    let manifest = json!({ "examples": entries });
    let manifest_text = format!(
        "{}\n",
        serde_json::to_string_pretty(&manifest).expect("json")
    );
    let manifest_path = page_dir().join("examples.json");
    let copies = page_dir().join("examples");
    if updating() {
        std::fs::create_dir_all(&copies).expect("dir");
        std::fs::write(&manifest_path, &manifest_text).expect("write manifest");
        for (stem, ..) in GALLERY {
            std::fs::copy(
                examples_dir().join(format!("{stem}.rhai")),
                copies.join(format!("{stem}.rhai")),
            )
            .expect("copy");
        }
        for (dir, ..) in GALLERY_PROJECTS {
            for path in project_files(&project_gallery_dir().join(dir)) {
                let target = copies.join(dir).join(&path);
                std::fs::create_dir_all(target.parent().expect("parent")).expect("dir");
                std::fs::copy(project_gallery_dir().join(dir).join(&path), target).expect("copy");
            }
        }
        return;
    }
    let committed = std::fs::read_to_string(&manifest_path).unwrap_or_default();
    assert_eq!(
        committed, manifest_text,
        "web/live/examples.json drifted from the gallery table: run with \
         POPKERN_UPDATE_GALLERY=1"
    );
    let committed: Value = serde_json::from_str(&committed).expect("manifest json");
    for (stem, ..) in GALLERY {
        let source = std::fs::read(examples_dir().join(format!("{stem}.rhai"))).expect("read");
        let copy = std::fs::read(copies.join(format!("{stem}.rhai"))).unwrap_or_default();
        assert_eq!(
            source, copy,
            "[{stem}] web/live/examples copy drifted: run with POPKERN_UPDATE_GALLERY=1"
        );
        assert!(
            copies.join(format!("{stem}.png")).exists(),
            "[{stem}] thumbnail missing: run with POPKERN_UPDATE_GALLERY=1"
        );
    }
    for (dir, ..) in GALLERY_PROJECTS {
        for path in project_files(&project_gallery_dir().join(dir)) {
            let source = std::fs::read(project_gallery_dir().join(dir).join(&path)).expect("read");
            let copy = std::fs::read(copies.join(dir).join(&path)).unwrap_or_default();
            assert_eq!(
                source, copy,
                "[{dir}/{path}] web/live/examples copy drifted: run with POPKERN_UPDATE_GALLERY=1"
            );
        }
    }
    assert_eq!(
        committed["examples"].as_array().map(Vec::len),
        Some(GALLERY.len() + GALLERY_PROJECTS.len())
    );
}
