//! PK-1010: `--stl` writes one binary STL of every visible body at its
//! world placement: the layout a slicer reads, watertight, and enclosing
//! the BRep's own volume to within the chord bound.

use std::path::{Path, PathBuf};

use popkern_cli::{Outcome, execute, parse_args};
use popkern_core::tolerance::Tolerance;

fn example(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../popkern-rhai/examples")
        .join(name)
        .canonicalize()
        .expect("example")
}

fn tmp(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("popkern_cli_stl");
    std::fs::create_dir_all(&dir).expect("dir");
    dir.join(name)
}

/// The facets of a binary STL, checking its layout.
fn read_stl(path: &Path) -> Vec<[[f32; 3]; 4]> {
    let bytes = std::fs::read(path).expect("stl bytes");
    assert!(bytes.starts_with(b"popkern binary STL"), "header");
    let count = u32::from_le_bytes(bytes[80..84].try_into().expect("count")) as usize;
    assert_eq!(bytes.len(), 84 + 50 * count, "the count matches the file");
    let f = |o: usize| f32::from_le_bytes(bytes[o..o + 4].try_into().expect("f32"));
    (0..count)
        .map(|k| {
            let base = 84 + 50 * k;
            let v = |i: usize| [f(base + 12 * i), f(base + 12 * i + 4), f(base + 12 * i + 8)];
            [v(0), v(1), v(2), v(3)]
        })
        .collect()
}

/// Enclosed volume, and how many facets' stored normals disagree with
/// their winding.
fn volume_and_normals(facets: &[[[f32; 3]; 4]]) -> (f64, usize) {
    let mut volume = 0.0;
    let mut disagree = 0;
    for [n, a, b, c] in facets {
        let p = |v: &[f32; 3]| [v[0] as f64, v[1] as f64, v[2] as f64];
        let (a, b, c) = (p(a), p(b), p(c));
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let cross = [
            u[1] * w[2] - u[2] * w[1],
            u[2] * w[0] - u[0] * w[2],
            u[0] * w[1] - u[1] * w[0],
        ];
        if cross[0] * n[0] as f64 + cross[1] * n[1] as f64 + cross[2] * n[2] as f64 <= 0.0 {
            disagree += 1;
        }
        volume += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))
            / 6.0;
    }
    (volume, disagree)
}

fn run_stl(script: &Path, out: &Path, extra: &[&str]) -> popkern_cli::RunSummary {
    let mut args: Vec<String> = vec![
        "run".into(),
        script.display().to_string(),
        "--stl".into(),
        out.display().to_string(),
    ];
    args.extend(extra.iter().map(|s| s.to_string()));
    let summary = execute(&parse_args(&args).expect("args"));
    assert!(matches!(summary.outcome, Outcome::Served), "{summary:?}");
    summary
}

/// The boss plate's BRep volume, from the script itself.
fn boss_plate_volume() -> f64 {
    let source = std::fs::read_to_string(example("boss_plate.rhai")).expect("source");
    let mut host = popkern_rhai::ScriptHost::new();
    let mut outcome = host.run(&source, "boss_plate.rhai").expect("runs");
    let mut document = outcome.take_document().expect("document");
    let part = document.to_part_all().expect("part");
    popkern_measure::body_volume(part.brep(), &Tolerance::default()).expect("volume")
}

#[test]
fn a_part_writes_a_watertight_stl_of_its_volume() {
    let exact = boss_plate_volume();
    let mut last = f64::INFINITY;
    for (name, extra) in [
        ("boss_default.stl", vec![]),
        (
            "boss_fine.stl",
            vec!["--stl-chord", "0.0005", "--stl-angle", "4"],
        ),
    ] {
        let out = tmp(name);
        let summary = run_stl(&example("boss_plate.rhai"), &out, &extra);
        let stl = summary.stl.expect("stl summary");
        assert_eq!(stl.open_edges, 0, "{name}: watertight");
        let facets = read_stl(&out);
        assert_eq!(facets.len(), stl.triangles);
        let (volume, disagree) = volume_and_normals(&facets);
        assert_eq!(disagree, 0, "{name}: every normal agrees with its winding");
        assert!(volume > 0.0, "{name}: facets face out");
        assert!(
            (volume - stl.volume).abs() <= 1e-4 * exact,
            "{name}: reported volume"
        );
        let rel = (volume - exact).abs() / exact;
        assert!(
            rel < 2e-3,
            "{name}: volume {volume} vs BRep {exact} (rel {rel:.1e})"
        );
        assert!(rel <= last, "{name}: a finer chord is no further off");
        last = rel;
    }
}

#[test]
fn an_assembly_writes_every_placed_occurrence() {
    let out = tmp("hinge.stl");
    let summary = run_stl(&example("hinge_assembly.rhai"), &out, &[]);
    let stl = summary.stl.expect("stl summary");
    assert!(stl.bodies >= 2, "occurrences placed: {}", stl.bodies);
    assert_eq!(stl.open_edges, 0, "watertight");
    let facets = read_stl(&out);
    let (volume, disagree) = volume_and_normals(&facets);
    assert_eq!(disagree, 0);
    assert!(volume > 0.0 && (volume - stl.volume).abs() <= 1e-6 * volume.abs().max(1.0));
}

#[test]
fn stl_options_must_be_positive() {
    for bad in [["--stl-chord", "0"], ["--stl-angle", "-3"]] {
        let args: Vec<String> = ["run", "x.rhai", bad[0], bad[1]]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(parse_args(&args).is_err(), "{bad:?}");
    }
}

#[test]
fn a_project_build_writes_its_stl() {
    let dir = tmp("stl_project");
    let _ = std::fs::remove_dir_all(&dir);
    let args = |a: &[&str]| -> Vec<String> { a.iter().map(|s| s.to_string()).collect() };
    let scaffold = parse_args(&args(&["new", dir.to_str().expect("utf8")])).expect("args");
    assert_eq!(popkern_cli::project::new_command(&scaffold), 0);
    let out = tmp("stl_project.stl");
    let _ = std::fs::remove_file(&out);
    let build = parse_args(&args(&[
        "build",
        dir.to_str().expect("utf8"),
        "--stl",
        out.to_str().expect("utf8"),
    ]))
    .expect("args");
    assert_eq!(popkern_cli::project::build_command(&build), 0);
    let (volume, disagree) = volume_and_normals(&read_stl(&out));
    assert_eq!(disagree, 0);
    assert!(volume > 0.0);
}
