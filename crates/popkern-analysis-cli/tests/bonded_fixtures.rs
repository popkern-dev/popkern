//! Bonded bodies (PK-981): bodies that share face area analysed as the
//! one body their fusion makes, under `bodies = "bonded"`.
//!
//! A 40 × 10 × 10 bar cut at x = 20 into two blocks, clamped at x = 0
//! and loaded (or heated) at x = 40. Bonded, the pair answers as the
//! uncut bar does; separate (the default), the loaded block is a free
//! body and the run declines, which is the control that the bond, not
//! the geometry, carries the load. Blocks that overlap decline as an
//! interference; blocks touching only along an edge share no area and
//! stay apart.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use popkern_analysis_cli::run::run;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepProductNode, UnitContext, export_step,
};
use serde_json::Value;

/// The block `[x0, x0 + 20] × [y0, y0 + 10] × [0, 10]`, or `length` long.
fn block(x0: f64, y0: f64, length: f64) -> BRep {
    let square = [
        Point3::new(x0, y0, 0.0),
        Point3::new(x0 + length, y0, 0.0),
        Point3::new(x0 + length, y0 + 10.0, 0.0),
        Point3::new(x0, y0 + 10.0, 0.0),
    ];
    BRep::prism(
        &square,
        Vec3::new(0.0, 0.0, 1.0),
        10.0,
        Tolerance::default(),
    )
    .expect("block")
}

/// `body` written as the one product `pair` of a STEP file in the test
/// temporary directory.
fn step_file(name: &str, body: BRep) -> PathBuf {
    let doc = StepDocument {
        products: vec![StepProductNode {
            name: "pair".to_string(),
            id: "pair".to_string(),
            children: Vec::new(),
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![body],
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    let bytes = export_step(&doc, &StepExportOptions::default()).expect("export");
    let path = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("bonded_{name}.step"));
    std::fs::write(&path, bytes).expect("write");
    path
}

fn pair(name: &str, b_x: f64, b_y: f64) -> PathBuf {
    let body = popkern_bool::concatenate_solids(block(0.0, 0.0, 20.0), block(b_x, b_y, 20.0))
        .expect("two solids");
    step_file(name, body)
}

const STRUCTURAL: &str = r#"
units = "mm"
material = "PLA (FDM)"

[faces]
root = { plane = { normal = [-1, 0, 0], through = [0, 0, 0] } }
tip = { plane = { normal = [1, 0, 0], through = [40, 0, 0] } }

[[support]]
faces = ["root"]
kind = "fixed"

[[load]]
faces = ["tip"]
force_n = 20
direction = [0, 0, -1]
"#;

const THERMAL: &str = r#"
analysis = "thermal"
units = "mm"
material = "Aluminium 6061-T6"

[faces]
root = { plane = { normal = [-1, 0, 0], through = [0, 0, 0] } }
tip = { plane = { normal = [1, 0, 0], through = [40, 0, 0] } }

[[fixed_temperature]]
faces = ["root"]
temperature_c = 20

[[heat_flux]]
faces = ["tip"]
power_w = 1
"#;

/// The report of `step` under `setup`, with `bodies` set when given.
fn solve(step: &Path, setup: &str, bodies: Option<&str>, tag: &str) -> Value {
    let text = match bodies {
        Some(b) => format!("bodies = \"{b}\"\n{setup}"),
        None => setup.to_string(),
    };
    let path = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(format!("bonded_{tag}.toml"));
    std::fs::write(&path, text).expect("write setup");
    run(step, &path).report
}

fn status(report: &Value) -> &str {
    report["status"].as_str().unwrap_or("")
}

#[test]
fn a_bonded_pair_answers_as_the_uncut_bar() {
    let bar = step_file("bar", block(0.0, 0.0, 40.0));
    let pair = pair("pair", 20.0, 0.0);
    let uncut = solve(&bar, STRUCTURAL, None, "bar");
    assert_eq!(status(&uncut), "solved", "{uncut}");
    let bonded = solve(&pair, STRUCTURAL, Some("bonded"), "pair_bonded");
    assert_eq!(status(&bonded), "solved", "{bonded}");
    let interfaces = &bonded["model"]["bonded_interfaces"];
    assert_eq!(interfaces.as_array().map(Vec::len), Some(1), "{interfaces}");
    assert_eq!(
        interfaces[0]["bodies"],
        serde_json::json!(["pair#1", "pair#2"])
    );
    let area = interfaces[0]["area_mm2"].as_f64().unwrap_or(f64::NAN);
    assert!((area - 100.0).abs() < 1e-9, "{area}");
    let tip = |r: &Value| {
        r["result"]["max_displacement_mm"]
            .as_f64()
            .unwrap_or(f64::NAN)
    };
    let (u_bar, u_pair) = (tip(&uncut), tip(&bonded));
    assert!(
        (u_pair - u_bar).abs() <= 1e-6 * u_bar,
        "bonded {u_pair} vs uncut {u_bar}"
    );
    // Control: separate, the loaded block is free.
    let separate = solve(&pair, STRUCTURAL, None, "pair_separate");
    assert_eq!(status(&separate), "declined", "{separate}");
    assert_eq!(separate["error"]["kind"], "FreeMotion", "{separate}");
}

#[test]
fn a_bonded_pair_conducts_as_the_uncut_bar() {
    let bar = step_file("bar_thermal", block(0.0, 0.0, 40.0));
    let pair = pair("pair_thermal", 20.0, 0.0);
    let uncut = solve(&bar, THERMAL, None, "bar_thermal");
    assert_eq!(status(&uncut), "solved", "{uncut}");
    let bonded = solve(&pair, THERMAL, Some("bonded"), "pair_thermal_bonded");
    assert_eq!(status(&bonded), "solved", "{bonded}");
    let hottest = |r: &Value| {
        r["result"]["max_temperature_c"]
            .as_f64()
            .unwrap_or(f64::NAN)
    };
    let (t_bar, t_pair) = (hottest(&uncut), hottest(&bonded));
    assert!(
        (t_pair - t_bar).abs() <= 1e-6 * (t_bar - 20.0),
        "bonded {t_pair} vs uncut {t_bar}"
    );
    // Control: separate, the heated block has nothing to set its level.
    let separate = solve(&pair, THERMAL, None, "pair_thermal_separate");
    assert_eq!(status(&separate), "declined", "{separate}");
    assert_eq!(
        separate["error"]["kind"], "FloatingTemperature",
        "{separate}"
    );
    // The weakest convection alone sets a separate body's level, so the
    // floating-temperature check passes it.
    let cooled = format!(
        "{THERMAL}\n[[convection]]\nfaces = \"all_other\"\npreset = \"still_air_bright_metal\"\nambient_c = 20\n"
    );
    let cooled = solve(&pair, &cooled, None, "pair_thermal_cooled");
    assert_eq!(status(&cooled), "solved", "{cooled}");
}

#[test]
fn overlapping_bodies_decline_as_an_interference() {
    let overlap = pair("overlap", 19.0, 0.0);
    let report = solve(&overlap, STRUCTURAL, Some("bonded"), "overlap");
    assert_eq!(status(&report), "declined", "{report}");
    let message = report["error"]["message"].as_str().unwrap_or("");
    assert!(
        message.contains("pair#1") && message.contains("pair#2") && message.contains("overlap"),
        "{message}"
    );
}

#[test]
fn bodies_touching_along_an_edge_are_not_bonded() {
    let edge = pair("edge", 20.0, 10.0);
    let report = solve(&edge, STRUCTURAL, Some("bonded"), "edge");
    assert_eq!(
        report["model"]["bonded_interfaces"],
        serde_json::json!([]),
        "{report}"
    );
    assert_eq!(status(&report), "declined", "{report}");
    assert_eq!(report["error"]["kind"], "FreeMotion", "{report}");
}
