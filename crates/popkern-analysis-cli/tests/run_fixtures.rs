//! The `run` command (PH-028): reports on both fixtures validate against
//! the committed schema, the `.vtu` round-trips through a reader, every
//! decline exits non-zero with its typed error, and reports are
//! bit-identical across runs.

use std::path::PathBuf;

use popkern_analysis_cli::run::{REPORT_SCHEMA, RunOptions, run, run_with};
use popkern_analysis_cli::step_model::model_from_step_file;
use serde_json::Value;

#[path = "support/report.rs"]
mod report;
use report::{assert_valid, data_array, repo, schema, validate};

#[test]
fn both_fixtures_report_to_schema_and_their_vtu_round_trips() {
    for (name, field) in [
        ("block_two_bores", "von_mises"),
        ("block_dimple", "temperature"),
    ] {
        let step = repo(&format!("tests/fixtures/analysis/{name}.step"));
        let out = run(
            &step,
            &repo(&format!("tests/fixtures/analysis/{name}.setup.toml")),
        );
        assert!(out.solved(), "{name}: {}", out.report);
        assert_eq!(out.report["schema"], REPORT_SCHEMA);
        assert_valid(&out.report);
        let vtu = out.vtu.unwrap_or_else(|| panic!("{name}: no vtu"));
        let model = model_from_step_file(&step).unwrap_or_else(|e| panic!("{e}"));
        let s = &model.solids[0];
        // Points and connectivity read back bit for bit.
        let points = data_array(&vtu, "points");
        let expect: Vec<f64> = s.vertices.iter().flatten().copied().collect();
        assert_eq!(points.len(), expect.len());
        assert!(
            points
                .iter()
                .zip(&expect)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
        let conn = data_array(&vtu, "connectivity");
        let tris: Vec<f64> = s
            .triangles
            .iter()
            .flatten()
            .map(|&v| f64::from(v))
            .collect();
        assert_eq!(conn, tris);
        let faces = data_array(&vtu, "face");
        assert_eq!(faces.len(), s.triangles.len());
        // The field reads back with one value per vertex, and its
        // extreme matches the report.
        let values = data_array(&vtu, field);
        assert_eq!(values.len(), s.vertices.len());
        if field == "temperature" {
            let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let reported = out.report["result"]["max_temperature_c"].as_f64();
            assert_eq!(Some(max), reported, "{name}");
            // The heat flux field matches the reported maximum, and the
            // coefficient used is reported with its preset.
            let q = data_array(&vtu, "heat_flux");
            assert_eq!(q.len(), 3 * s.vertices.len());
            let qmax = q
                .chunks(3)
                .map(|c| (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt())
                .fold(0.0, f64::max);
            let reported = out.report["result"]["max_heat_flux_w_per_m2"].as_f64();
            assert_eq!(Some(qmax), reported, "{name}");
            assert_eq!(
                out.report["result"]["convection"],
                serde_json::json!([{
                    "faces": "all_other",
                    "h_w_per_m2k": 25.0,
                    "source": "fan",
                    "ambient_c": 22.0,
                }])
            );
        } else {
            let d = data_array(&vtu, "displacement");
            let max = d
                .chunks(3)
                .map(|c| (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt())
                .fold(0.0, f64::max);
            // The report searches the cell-resolution surface samples
            // too, so it is at least the display vertices' largest
            // (block_two_bores: 0.0328 against 0.0314 at the vertices).
            let reported = out.report["result"]["max_displacement_mm"]
                .as_f64()
                .unwrap_or(0.0);
            assert!(
                reported >= max && reported <= 1.1 * max,
                "{name}: {reported} vs {max}"
            );
            let zz = data_array(&vtu, "zz_error");
            assert!(zz.iter().all(|&e| e == -1.0 || e >= 0.0));
        }
    }
}

#[test]
fn the_schema_check_fails_on_a_malformed_report() {
    let out = run(
        &repo("tests/fixtures/analysis/block_two_bores.step"),
        &repo("tests/fixtures/analysis/block_two_bores.setup.toml"),
    );
    let s = schema();
    for (path, bad) in [
        ("/schema", Value::from("popkern-analysis-report/2")),
        ("/result/assessment/safety/estimated", Value::from(false)),
        ("/solve/dof", Value::from("many")),
    ] {
        let mut r = out.report.clone();
        *r.pointer_mut(path).unwrap_or_else(|| panic!("{path}")) = bad;
        let mut errors = Vec::new();
        validate(&s, &s, &r, "$", &mut errors);
        assert!(!errors.is_empty(), "{path} accepted");
    }
}

#[test]
fn reports_are_bit_identical_across_runs() {
    let step = repo("tests/fixtures/analysis/block_two_bores.step");
    let setup = repo("tests/fixtures/analysis/block_two_bores.setup.toml");
    let a = run(&step, &setup);
    let b = run(&step, &setup);
    let text = |v: &Value| serde_json::to_string_pretty(v).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(text(&a.report), text(&b.report));
    assert_eq!(a.vtu, b.vtu);
}

#[test]
fn every_decline_exits_non_zero_with_its_typed_error() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("run_declines");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let header = "units = \"mm\"\nmaterial = \"PLA (FDM)\"\n\n[faces]\n\
                  base = { plane = { normal = [0, 0, -1], through = [0, 0, 0] } }\n\
                  top = { plane = { normal = [0, 0, 1], through = [0, 0, 10] } }\n";
    let load = "\n[[load]]\nfaces = [\"top\"]\nforce_n = 100\ndirection = [1, 0, 0]\n";
    let cases = [
        (
            "unresolved",
            "units = \"mm\"\nmaterial = \"PLA (FDM)\"\n\n[faces]\n\
             base = { plane = { normal = [0, 0, -1], through = [0, 0, 3] } }\n\
             top = { plane = { normal = [0, 0, 1], through = [0, 0, 10] } }\n\
             \n[[support]]\nfaces = [\"base\"]\nkind = \"fixed\"\n"
                .to_string()
                + load,
            "SelectorUnresolved",
        ),
        ("no_support", header.to_string() + load, "NoSupport"),
        (
            "rollers_only",
            header.to_string()
                + "\n[[support]]\nfaces = [\"base\"]\nkind = \"roller\"\nnormal = [0, 0, -1]\n"
                + load,
            "FreeMotion",
        ),
    ];
    let bin = env!("CARGO_BIN_EXE_popkern-analysis");
    for (name, text, kind) in cases {
        let setup = dir.join(format!("{name}.toml"));
        std::fs::write(&setup, &text).unwrap_or_else(|e| panic!("{e}"));
        let status = std::process::Command::new(bin)
            .args(["run"])
            .arg(repo("tests/fixtures/analysis/block_two_bores.step"))
            .arg(&setup)
            .args(["--out"])
            .arg(&dir)
            .status()
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(status.code(), Some(2), "{name}");
        let report: Value = serde_json::from_str(
            &std::fs::read_to_string(dir.join(format!("{name}.report.json")))
                .unwrap_or_else(|e| panic!("{e}")),
        )
        .unwrap_or_else(|e| panic!("{e}"));
        assert_valid(&report);
        assert_eq!(report["status"], "declined", "{name}");
        assert_eq!(report["error"]["kind"], kind, "{name}: {}", report["error"]);
        assert!(
            !dir.join(format!("{name}.vtu")).exists(),
            "{name}: vtu written"
        );
    }
}

/// Runs block_two_bores clamped on its base and sheared on top, with
/// `extra` support text, and returns the worst radial displacement on
/// the bore at (2.5, 5) over the worst displacement anywhere.
fn bore_radial_ratio(name: &str, extra: &str) -> Value {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("bore_roller");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let text = format!(
        "units = \"mm\"\nmaterial = \"PLA (FDM)\"\n\n[faces]\n\
         base = {{ plane = {{ normal = [0, 0, -1], through = [0, 0, 0] }} }}\n\
         top = {{ plane = {{ normal = [0, 0, 1], through = [0, 0, 10] }} }}\n\
         bore = {{ cylinder = {{ axis = [0, 0, 1], radius = 1.5, through = [2.5, 5, 0] }} }}\n\
         \n[[support]]\nfaces = [\"base\"]\nkind = \"fixed\"\n{extra}\
         \n[[load]]\nfaces = [\"top\"]\nforce_n = 100\ndirection = [1, 0, 0]\n"
    );
    let setup = dir.join(format!("{name}.toml"));
    std::fs::write(&setup, text).unwrap_or_else(|e| panic!("{e}"));
    let out = run(
        &repo("tests/fixtures/analysis/block_two_bores.step"),
        &setup,
    );
    if !out.solved() {
        return out.report;
    }
    let vtu = out.vtu.unwrap_or_else(|| panic!("no vtu"));
    let points = data_array(&vtu, "points");
    let d = data_array(&vtu, "displacement");
    let (mut radial, mut all) = (0.0f64, 0.0f64);
    for (p, u) in points.chunks(3).zip(d.chunks(3)) {
        all = all.max((u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt());
        let (dx, dy) = (p[0] - 2.5, p[1] - 5.0);
        let r = dx.hypot(dy);
        if (r - 1.5).abs() < 1e-6 && p[2] > 1e-6 {
            radial = radial.max(((u[0] * dx + u[1] * dy) / r).abs());
        }
    }
    serde_json::json!(radial / all)
}

#[test]
fn a_surface_roller_holds_a_bore_radially_and_a_wrong_normal_declines() {
    // Measured: worst radial displacement on the bore over the worst
    // anywhere, 0.976 free (the bore moves with the shear) and 5.6e-3 on
    // the roller (weakly imposed, read at display vertices).
    let free = bore_radial_ratio("free", "").as_f64().unwrap_or(f64::NAN);
    let held = bore_radial_ratio(
        "held",
        "\n[[support]]\nfaces = [\"bore\"]\nkind = \"roller\"\n",
    )
    .as_f64()
    .unwrap_or(f64::NAN);
    assert!(free > 0.5, "free {free}");
    assert!(held < 0.01 * free, "held {held} vs free {free}");
    // Negative control: `normal` asserts a planar face's normal; on the
    // bore it cannot hold, and the setup declines by name.
    let report = bore_radial_ratio(
        "wrong",
        "\n[[support]]\nfaces = [\"bore\"]\nkind = \"roller\"\nnormal = [1, 0, 0]\n",
    );
    assert_eq!(report["status"], "declined", "{report}");
    let message = report["error"]["message"].as_str().unwrap_or("");
    assert!(message.contains("not the normal"), "{message}");
}

#[test]
fn the_refinement_check_labels_the_clamped_edge_as_a_support_singularity() {
    // block_two_bores clamped on its base and sheared, at 1 mm cells:
    // the clamp's edge is a singularity of the idealised support (the
    // stress there grows without bound), and h/2 shows it. Measured:
    // growth 1.264 in place, peak at (10, 0, 1), one cell above the
    // clamped corner.
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("refine_check");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let base = std::fs::read_to_string(repo("tests/fixtures/analysis/block_two_bores.setup.toml"))
        .unwrap_or_else(|e| panic!("{e}"));
    let setup = dir.join("coarse.toml");
    std::fs::write(&setup, base + "\n[options]\ntarget_cell_size = 1.0\n")
        .unwrap_or_else(|e| panic!("{e}"));
    let step = repo("tests/fixtures/analysis/block_two_bores.step");
    let plain = run(&step, &setup);
    assert_valid(&plain.report);
    assert!(plain.report["result"]["refinement"].is_null());
    let out = run_with(
        &step,
        &setup,
        &RunOptions {
            refinement_check: true,
            ..RunOptions::default()
        },
    );
    assert_valid(&out.report);
    let refinement = &out.report["result"]["refinement"];
    assert_eq!(refinement["status"], "solved", "{refinement}");
    assert_eq!(refinement["cell_size"], 0.5);
    let ratio = refinement["growth"][0]["ratio"].as_f64().unwrap_or(0.0);
    let peak = &out.report["result"]["assessment"]["peaks"][0];
    assert!(ratio > 1.15, "{refinement}");
    // At the clamp's edge: a support label, not "add a fillet".
    assert_eq!(
        peak["singularity"]["SupportEdge"]["ratio"].as_f64(),
        Some(ratio),
        "{peak}"
    );
    // The main result is the coarse solve, unchanged by the check.
    assert_eq!(
        out.report["result"]["max_displacement_mm"],
        plain.report["result"]["max_displacement_mm"]
    );
    assert_eq!(out.vtu, plain.vtu);
}

#[test]
fn the_png_is_deterministic_and_its_legend_spans_the_field() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("png");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let step = repo("tests/fixtures/analysis/block_two_bores.step");
    let setup = repo("tests/fixtures/analysis/block_two_bores.setup.toml");
    let mut bytes = Vec::new();
    for k in 0..2 {
        let path = dir.join(format!("run{k}.png"));
        let options = RunOptions {
            png: Some(path.clone()),
            ..RunOptions::default()
        };
        let out = run_with(&step, &setup, &options);
        assert_valid(&out.report);
        let png = &out.report["result"]["png"];
        assert_eq!(png["field"], "von_mises_mpa");
        // The legend's range is the vertex von Mises range of the .vtu.
        let vm = data_array(&out.vtu.unwrap_or_else(|| panic!("no vtu")), "von_mises");
        let (lo, hi) = vm
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), &v| {
                (l.min(v), h.max(v))
            });
        assert_eq!(
            (png["min"].as_f64(), png["max"].as_f64()),
            (Some(lo), Some(hi))
        );
        bytes.push(std::fs::read(&path).unwrap_or_else(|e| panic!("{e}")));
    }
    assert_eq!(bytes[0], bytes[1], "two runs gave different PNG bytes");
    // A PNG, 800 × 600 (IHDR width and height, big-endian).
    let b = &bytes[0];
    assert_eq!(&b[..8], b"\x89PNG\r\n\x1a\n");
    assert_eq!(&b[16..24], &[0, 0, 3, 32, 0, 0, 2, 88]);
}

#[test]
fn a_build_direction_is_used_for_a_printed_material() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("build_direction");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let base = std::fs::read_to_string(repo("tests/fixtures/analysis/block_two_bores.setup.toml"))
        .unwrap_or_else(|e| panic!("{e}"));
    let setup = dir.join("z.toml");
    let text = base.replace(
        "material = \"PLA (FDM)\"",
        "material = \"PLA (FDM)\"\nbuild_direction = [0, 0, 1]",
    );
    std::fs::write(&setup, text).unwrap_or_else(|e| panic!("{e}"));
    let out = run(
        &repo("tests/fixtures/analysis/block_two_bores.step"),
        &setup,
    );
    assert_valid(&out.report);
    let result = &out.report["result"];
    assert_eq!(result["build_direction_used"], true);
    let safety = &result["assessment"]["safety"];
    assert_eq!(
        safety["build_direction"],
        serde_json::json!([0.0, 0.0, 1.0])
    );
    // PLA's S_layer is sourced (PH-080): layer shear is checked.
    assert_eq!(safety["unchecked"], serde_json::json!([]));
    let checks = safety["checks"].as_array().cloned().unwrap_or_default();
    assert!(
        checks.iter().any(|c| c["check"] == "InterLayerShear"),
        "{checks:?}"
    );
    // Without it the worst-case orientation is reported as before.
    let plain = run(
        &repo("tests/fixtures/analysis/block_two_bores.step"),
        &repo("tests/fixtures/analysis/block_two_bores.setup.toml"),
    );
    assert_eq!(plain.report["result"]["build_direction_used"], false);
}
