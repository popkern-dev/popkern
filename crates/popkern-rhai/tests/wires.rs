//! Wire primitives (PK-1083): each is the exact curve its arguments
//! name (length against the closed form), replays from its emitted
//! script, and is a path: a pipe along it encloses π r² L (the tube
//! formula), a profile swept along it the section's area times L.

use popkern_rhai::{ScriptHost, to_rhai};
use std::f64::consts::PI;

fn run(source: &str) -> Result<popkern_script::Model, String> {
    let mut host = ScriptHost::with_options(Default::default());
    host.run(source, "wires")
        .map_err(|e| e.to_string())?
        .model
        .ok_or_else(|| "no model".into())
}

/// (total wire length, total solid volume) over the model's bodies.
fn measure(model: &mut popkern_script::Model) -> (f64, f64) {
    let mut length = 0.0;
    for wire in model.wire_bodies() {
        let brep = model.brep_of(wire).expect("wire");
        for eh in brep.edges.keys() {
            length += popkern_measure::edge_length(&brep, eh).expect("length");
        }
    }
    let mut volume = 0.0;
    if !model.alive_bodies().is_empty() {
        for p in &model.to_parts_all().expect("parts") {
            let brep = p.brep();
            volume += popkern_measure::body_volume(brep, &brep.tolerance).expect("volume");
        }
    }
    (length, volume)
}

/// `body` (a wire, then optionally a solid built on it) measures
/// `length` and `volume` to `rel`, and so does its emitted script.
fn assert_serves(body: &str, length: f64, volume: f64, rel: f64) {
    let source = format!("let m = model();\n{body}\nm\n");
    let mut model = run(&source).unwrap_or_else(|e| panic!("{e}\n{source}"));
    let (l, v) = measure(&mut model);
    assert!(
        (l - length).abs() <= rel * length,
        "{body}: length {l} vs {length}"
    );
    assert!(
        (v - volume).abs() <= rel * volume.max(1e-300),
        "{body}: volume {v} vs {volume}"
    );
    let emitted = to_rhai(&mut model).expect("emit");
    let (l2, v2) = measure(&mut run(&emitted).unwrap_or_else(|e| panic!("{e}\n{emitted}")));
    assert!(
        (l - l2).abs() <= 1e-12 * l && (v - v2).abs() <= 1e-12 * v.max(1.0),
        "{emitted}"
    );
}

#[test]
fn lines_and_polylines() {
    assert_serves(
        "let w = m.wire_line([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]);",
        5.0,
        0.0,
        1e-12,
    );
    assert_serves(
        "let w = m.wire_polyline([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [3.0, 4.0, 0.0]], #{ closed: true });",
        12.0,
        0.0,
        1e-12,
    );
}

#[test]
fn arcs_and_circles() {
    // A half circle of radius 2 through (0, 2).
    assert_serves(
        "let w = m.wire_arc([2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [-2.0, 0.0, 0.0]);",
        2.0 * PI,
        0.0,
        1e-12,
    );
    // Three quarters: through (0, -2) from (2, 0) to (0, 2) the long way.
    assert_serves(
        "let w = m.wire_arc([2.0, 0.0, 0.0], [0.0, -2.0, 0.0], [0.0, 2.0, 0.0]);",
        3.0 * PI,
        0.0,
        1e-12,
    );
    assert_serves(
        "let w = m.wire_circle(3.0, #{ at: [0.0, 0.0, 5.0], axis: [1.0, 0.0, 0.0], ref: [0.0, 1.0, 0.0] });",
        6.0 * PI,
        0.0,
        1e-12,
    );
}

#[test]
fn helices_and_splines() {
    let l = 3.0 * (4.0 * PI * PI * 4.0 + 1.0f64).sqrt();
    assert_serves("let w = m.wire_helix(2.0, 1.0, 3.0);", l, 0.0, 1e-7);
    assert_serves(
        "let w = m.wire_helix(2.0, 1.0, 3.0, #{ left_hand: true });",
        l,
        0.0,
        1e-7,
    );
    // Collinear points: the interpolant is the straight segment.
    assert_serves(
        "let w = m.wire_spline([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [4.0, 4.0, 4.0]]);",
        4.0 * 3f64.sqrt(),
        0.0,
        1e-9,
    );
    // A quadratic Bézier from its control polygon: (0,0) (1,2) (2,0)
    // is the parabola y = 2x - x², of length sqrt 5 + asinh(2) / 2.
    assert_serves(
        "let w = m.wire_spline([[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]], #{ degree: 2, control: true });",
        5f64.sqrt() + 2f64.asinh() / 2.0,
        0.0,
        1e-12,
    );
}

/// Wires are paths: a pipe along a line, a circle and an arc encloses
/// π r² L; a square swept along a line its area times the length.
#[test]
fn wires_are_paths() {
    assert_serves(
        "let w = m.wire_line([0.0, 0.0, 0.0], [0.0, 0.0, 10.0]);\nlet t = w.pipe(0.5);",
        10.0,
        PI * 0.25 * 10.0,
        1e-9,
    );
    assert_serves(
        "let w = m.wire_circle(4.0);\nlet t = w.pipe(1.0);",
        8.0 * PI,
        2.0 * PI * PI * 4.0,
        1e-9,
    );
    assert_serves(
        "let w = m.wire_arc([4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [-4.0, 0.0, 0.0]);\nlet t = w.pipe(1.0);",
        4.0 * PI,
        PI * 4.0 * PI,
        1e-9,
    );
    assert_serves(
        "let w = m.wire_line([0.0, 0.0, 0.0], [10.0, 0.0, 0.0]);\nlet s = m.profile(xy()).rect(2.0, 3.0).sketch(\"sq\").sweep(w);",
        10.0,
        60.0,
        1e-9,
    );
}

#[test]
fn bad_wires_are_named() {
    for (body, needle) in [
        ("m.wire_line([1.0, 1.0, 1.0], [1.0, 1.0, 1.0])", "coincide"),
        (
            "m.wire_arc([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0])",
            "collinear",
        ),
        (
            "m.wire_spline([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])",
            "at least 4 points",
        ),
        ("m.wire_helix(1.0, 0.0, 2.0)", "positive"),
    ] {
        let source = format!("let m = model();\nlet w = {body};\nm\n");
        match run(&source) {
            Ok(_) => panic!("{body} served"),
            Err(e) => assert!(e.contains(needle), "{body}: {e}"),
        }
    }
}
