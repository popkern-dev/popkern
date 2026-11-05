//! Sweeps along paths of several edges (brep-bench §sweep): the runs
//! mitred or rounded at sharp corners, one frame carried round them,
//! each replayed from its emitted script.

use popkern_rhai::{ScriptHost, to_rhai};
use std::f64::consts::PI;

fn run(source: &str) -> Result<popkern_script::Model, String> {
    let mut host = ScriptHost::with_options(Default::default());
    host.run(source, "paths")
        .map_err(|e| e.to_string())?
        .model
        .ok_or_else(|| "no model".into())
}

fn measure(model: &mut popkern_script::Model) -> (f64, f64) {
    let parts = model.to_parts_all().expect("parts");
    let v = parts
        .iter()
        .map(|p| popkern_measure::body_volume(p.brep(), &p.brep().tolerance).expect("volume"))
        .sum();
    let a = parts
        .iter()
        .map(|p| popkern_measure::body_surface_area(p.brep(), &p.brep().tolerance).expect("area"))
        .sum();
    (v, a)
}

fn assert_serves(body: &str, volume: f64, area: Option<f64>) {
    let source = format!("let m = model();\n{body}\nm\n");
    let mut model = run(&source).unwrap_or_else(|e| panic!("{e}\n{source}"));
    let (v, a) = measure(&mut model);
    assert!(
        (v - volume).abs() <= 1e-9 * volume,
        "{body}: volume {v} vs {volume}"
    );
    if let Some(area) = area {
        assert!(
            (a - area).abs() <= 1e-9 * area,
            "{body}: area {a} vs {area}"
        );
    }
    let emitted = to_rhai(&mut model).expect("emit");
    let (v2, _) = measure(&mut run(&emitted).unwrap_or_else(|e| panic!("{e}\n{emitted}")));
    assert!((v - v2).abs() <= 1e-12 * v, "{emitted}");
}

/// brep-bench sweep/mitred-elbow: a 10 × 10 bar round a 60° corner,
/// mitred: w²L = 8000, 4wL + 2w² = 3400.
#[test]
fn a_bar_mitred_round_a_sharp_corner() {
    let s3 = 3f64.sqrt();
    assert_serves(
        &format!(
            "let w = m.wire_polyline([[0.0, 0.0, 0.0], [40.0, 0.0, 0.0], [20.0, {}, 0.0]]);\nlet bar = m.profile(yz()).rect(10.0, 10.0).sketch(\"sq\").sweep(w);",
            20.0 * s3
        ),
        8000.0,
        Some(3400.0),
    );
}

/// A square bar rounded round a right angle: the runs (3000 each)
/// overlap in a 5 × 5 × 10 block inside the corner, and the outer half
/// of the section turns into a quarter cylinder of radius 5, height 10.
/// (A round section's corner is a quarter ball tangent to both runs,
/// which the union does not yet serve: PK-1089.)
#[test]
fn a_bar_rounded_round_a_corner() {
    assert_serves(
        "let w = m.wire_polyline([[-30.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 30.0, 0.0]]);\nlet p = m.profile(yz()).rect(10.0, 10.0).sketch(\"sq\").sweep(w, #{ corner: \"round\" });",
        6000.0 - 250.0 + 62.5 * PI,
        None,
    );
}

/// A closed rectangle path mitred at its four corners, the section
/// centred on it: the section's area times the path's length.
#[test]
fn a_closed_path_is_mitred_at_every_corner() {
    assert_serves(
        "let w = m.wire_polyline([[-30.0, -20.0, 0.0], [30.0, -20.0, 0.0], [30.0, 20.0, 0.0], [-30.0, 20.0, 0.0]], #{ closed: true });\nlet f = m.profile(yz()).rect(8.0, 6.0).sketch(\"bar\").sweep(w);",
        48.0 * 200.0,
        None,
    );
}

#[test]
fn a_u_turn_has_no_mitre() {
    let source = "let m = model();\nlet w = m.wire_polyline([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);\nlet b = m.profile(yz()).rect(1.0, 1.0).sketch(\"sq\").sweep(w);\nm\n";
    match run(source) {
        Ok(_) => panic!("a U-turn served"),
        Err(e) => assert!(e.contains("U-turn"), "{e}"),
    }
}
