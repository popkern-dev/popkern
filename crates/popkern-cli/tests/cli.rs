//! PK-406: `popkern run` end to end — exports, JSON report, parameter
//! overrides, exit codes per failure class, operation and wall-clock
//! budgets (a runaway script stops, the watcher survives it), and
//! rebuild-on-save keeping the last good outputs.

use std::path::{Path, PathBuf};
use std::time::Duration;

use popkern_cli::{Cli, Outcome, Watcher, execute, exit_code, parse_args};

fn boss_plate() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../popkern-rhai/examples/boss_plate.rhai")
        .canonicalize()
        .expect("boss plate example")
}

fn tmp(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("popkern_cli");
    std::fs::create_dir_all(&dir).expect("dir");
    dir.join(name)
}

fn parsed(args: &[&str]) -> Cli {
    let owned: Vec<String> = args.iter().map(|s| s.to_string()).collect();
    parse_args(&owned).unwrap_or_else(|e| panic!("{e}"))
}

#[test]
fn boss_plate_runs_and_exports_everything() {
    let (png, step, pkn) = (
        tmp("boss_plate.png"),
        tmp("boss_plate.step"),
        tmp("boss_plate.pkn"),
    );
    let script = boss_plate();
    let cli = parsed(&[
        "run",
        script.to_str().expect("utf8"),
        "--png",
        png.to_str().expect("utf8"),
        "--step",
        step.to_str().expect("utf8"),
        "--pkn",
        pkn.to_str().expect("utf8"),
        "--json",
    ]);
    let summary = execute(&cli);
    assert!(matches!(summary.outcome, Outcome::Served), "{summary:?}");
    assert_eq!(exit_code(&summary), 0);
    assert_eq!(summary.features, 6);
    assert_eq!(
        summary.printed.len(),
        1,
        "the example prints one describe line"
    );
    for path in [&png, &step, &pkn] {
        assert!(path.exists(), "{}", path.display());
        assert!(std::fs::metadata(path).expect("meta").len() > 0);
    }
    let json = serde_json::to_value(&summary).expect("json");
    assert_eq!(json["outcome"], "served");
    assert_eq!(json["report"]["failures"], 0);
    assert!(
        json["report"]["features"]
            .as_array()
            .is_some_and(|f| f.len() == 6)
    );
    assert!(json["elapsed_ms"].as_f64().is_some());
    let text = popkern_cli::report::text(&summary);
    assert!(text.starts_with("served "), "{text}");
    assert!(text.contains("rounds") && text.contains("fillet"), "{text}");
    // The saved document reopens through the recorder.
    let mut model = popkern_script::Model::open(&pkn).expect("open");
    assert_eq!(model.replay().expect("replay").failures().count(), 0);
    // OCC oracle on the STEP when the venv is present.
    let oracle = Path::new(&std::env::var("HOME").unwrap_or_default()).join("venvs/ocp/bin/python");
    let script = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle/occ_check.py");
    if oracle.exists() && script.exists() {
        let output = std::process::Command::new(&oracle)
            .arg(&script)
            .arg(&step)
            .output()
            .expect("oracle");
        let report = String::from_utf8_lossy(&output.stdout);
        assert!(report.contains("\"passed\": true"), "OCC oracle: {report}");
    } else {
        eprintln!("[pk406] OCC oracle absent; STEP round trip checked by popkern only");
    }
}

#[test]
fn parameter_override_replays_before_exporting() {
    let pkn = tmp("boss_plate_tall.pkn");
    let script = boss_plate();
    let cli = parsed(&[
        "run",
        script.to_str().expect("utf8"),
        "--param",
        "extrude-1.distance=0.5",
        "--pkn",
        pkn.to_str().expect("utf8"),
    ]);
    let summary = execute(&cli);
    assert!(matches!(summary.outcome, Outcome::Served), "{summary:?}");
    let model = popkern_script::Model::open(&pkn).expect("open");
    let distance = model
        .params()
        .into_iter()
        .find(|p| p.name == "extrude-1.distance")
        .expect("param");
    assert!((distance.value - 0.5).abs() < 1e-12);
    // An unknown parameter is a script error (exit 2), not a decline.
    let cli = parsed(&["run", script.to_str().expect("utf8"), "--param", "nope=1"]);
    let summary = execute(&cli);
    assert_eq!(exit_code(&summary), 2, "{summary:?}");
}

fn write_script(name: &str, body: &str) -> PathBuf {
    let path = tmp(name);
    std::fs::write(&path, body).expect("write");
    path
}

#[test]
fn exit_codes_distinguish_script_kernel_and_gate_failures() {
    let syntax = write_script("syntax.rhai", "let m = model(\n");
    assert_eq!(
        exit_code(&execute(&parsed(&["run", syntax.to_str().expect("utf8")]))),
        2
    );
    let decline = write_script(
        "decline.rhai",
        r#"let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
cube.fillet(edges(), 0.9);
m"#,
    );
    let summary = execute(&parsed(&["run", decline.to_str().expect("utf8")]));
    assert_eq!(exit_code(&summary), 3, "{summary:?}");
    match &summary.outcome {
        Outcome::KernelDecline { kind, payload, .. } => {
            assert_eq!(kind, "FeatureFailed");
            assert_eq!(payload["error"]["kind"], "BlendDeclined");
        }
        other => panic!("{other:?}"),
    }
    let ok = write_script(
        "ok.rhai",
        r#"let m = model();
m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
m"#,
    );
    let summary = execute(&parsed(&[
        "run",
        ok.to_str().expect("utf8"),
        "--pkn",
        "/nonexistent/dir/for/popkern/out.pkn",
    ]));
    assert_eq!(exit_code(&summary), 4, "{summary:?}");
    let missing = execute(&parsed(&["run", "/nonexistent/script.rhai"]));
    assert_eq!(exit_code(&missing), 4);
    assert!(parse_args(&["run".to_string()]).is_err());
    assert!(parse_args(&["frob".to_string()]).is_err());
}

#[test]
fn version_and_installed_page_layout() {
    assert_eq!(
        parsed(&["--version"]).command,
        popkern_cli::Command::Version
    );
    assert_eq!(parsed(&["version"]).command, popkern_cli::Command::Version);
    assert!(!popkern_cli::VERSION.is_empty());

    // <prefix>/bin/popkern finds <prefix>/share/popkern/live, also
    // through a symlink elsewhere (install.sh links ~/.local/bin).
    let prefix = tmp("installed_prefix");
    let _ = std::fs::remove_dir_all(&prefix);
    let page = prefix.join("share/popkern/live");
    std::fs::create_dir_all(prefix.join("bin")).expect("bin");
    std::fs::create_dir_all(&page).expect("page");
    let exe = prefix.join("bin/popkern");
    std::fs::write(&exe, b"").expect("exe");
    assert_eq!(popkern_cli::installed_page(&exe), None, "no index.html yet");
    std::fs::write(page.join("index.html"), b"<!doctype html>").expect("index");
    let found = popkern_cli::installed_page(&exe).expect("installed page");
    assert_eq!(found, page.canonicalize().expect("canonical"));
    #[cfg(unix)]
    {
        let link = prefix.join("popkern-link");
        let _ = std::fs::remove_file(&link);
        std::os::unix::fs::symlink(&exe, &link).expect("symlink");
        assert_eq!(popkern_cli::installed_page(&link), Some(found));
    }
}

#[test]
fn runaway_scripts_stop_within_their_budgets() {
    let spin = write_script(
        "spin.rhai",
        "let m = model();\nlet i = 0;\nloop { i += 1; }\n",
    );
    let started = std::time::Instant::now();
    let summary = execute(&parsed(&[
        "run",
        spin.to_str().expect("utf8"),
        "--max-ops",
        "200000",
    ]));
    assert!(
        started.elapsed() < Duration::from_secs(10),
        "ops budget must bite"
    );
    assert_eq!(exit_code(&summary), 2, "{summary:?}");
    match &summary.outcome {
        Outcome::ScriptError { kind, .. } => assert_eq!(kind, "BudgetExceeded"),
        other => panic!("{other:?}"),
    }
    let started = std::time::Instant::now();
    let summary = execute(&parsed(&[
        "run",
        spin.to_str().expect("utf8"),
        "--timeout",
        "0.3",
    ]));
    let elapsed = started.elapsed();
    assert!(
        elapsed < Duration::from_secs(5),
        "wall budget must bite, took {elapsed:?}"
    );
    match &summary.outcome {
        Outcome::ScriptError { kind, rendered, .. } => {
            assert_eq!(kind, "BudgetExceeded");
            assert!(rendered.contains("budget"), "{rendered}");
        }
        other => panic!("{other:?}"),
    }
    // An unbounded feature loop is bounded the same way.
    let pattern = write_script(
        "pattern.rhai",
        r#"let m = model();
let i = 0;
loop { m.profile(xy_at(i * 1.0)).rect(1.0, 1.0).sketch("s" + i).extrude(blind(0.5)); i += 1; }
"#,
    );
    let summary = execute(&parsed(&[
        "run",
        pattern.to_str().expect("utf8"),
        "--timeout",
        "1.0",
    ]));
    assert_eq!(exit_code(&summary), 2, "{summary:?}");
}

fn settle(watcher: &mut Watcher) -> usize {
    let before = watcher.history.len();
    for _ in 0..50 {
        watcher.check();
        if watcher.history.len() > before {
            break;
        }
        std::thread::sleep(Duration::from_millis(60));
    }
    watcher.history.len()
}

#[test]
fn watch_rebuilds_on_save_and_keeps_the_last_good_outputs() {
    let script = write_script(
        "watched.rhai",
        "let m = model();\nm.profile(xy()).rect(1.0, 1.0).sketch(\"cube\").extrude(blind(1.0));\nm\n",
    );
    let png = tmp("watched.png");
    let _ = std::fs::remove_file(&png);
    let mut watcher = Watcher::new(parsed(&[
        "run",
        script.to_str().expect("utf8"),
        "--png",
        png.to_str().expect("utf8"),
        "--max-ops",
        "500000",
    ]));
    assert_eq!(settle(&mut watcher), 1, "first run on start");
    assert!(watcher.last_served(), "{:?}", watcher.history.last());
    let first_png = std::fs::read(&png).expect("png written");
    // A broken save: reported, last good PNG kept.
    std::thread::sleep(Duration::from_millis(20));
    std::fs::write(&script, "let m = model(\n").expect("write");
    assert_eq!(settle(&mut watcher), 2, "rebuilt after the save");
    assert!(!watcher.last_served());
    assert_eq!(std::fs::read(&png).expect("png kept"), first_png);
    // A runaway save: the budget stops it and the watcher survives.
    std::thread::sleep(Duration::from_millis(20));
    std::fs::write(&script, "let m = model();\nloop {}\n").expect("write");
    assert_eq!(settle(&mut watcher), 3);
    assert!(matches!(
        watcher.history.last().map(|s| &s.outcome),
        Some(Outcome::ScriptError { kind, .. }) if kind == "BudgetExceeded"
    ));
    // A good save renders again.
    std::thread::sleep(Duration::from_millis(20));
    std::fs::write(
        &script,
        "let m = model();\nm.profile(xy()).rect(1.0, 2.0).sketch(\"slab\").extrude(blind(1.0));\nm\n",
    )
    .expect("write");
    assert_eq!(settle(&mut watcher), 4);
    assert!(watcher.last_served());
    assert_ne!(std::fs::read(&png).expect("png"), first_png, "re-rendered");
    if let Some(ms) = watcher.history.last().and_then(|s| s.elapsed_ms) {
        eprintln!("[pk406] edit-to-run on the slab: {ms:.0} ms (script only)");
    }
}
