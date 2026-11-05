//! A revolve up to a face (brep-bench §4.4): a rib turned until it
//! meets a stop on another body, joined to it. Each profile point's
//! angle is measured from the point itself (the circle's own parameter
//! starts on an arbitrary perpendicular of the axis).

use popkern_rhai::{ScriptHost, to_rhai};
use std::f64::consts::PI;

fn run(source: &str) -> Result<popkern_script::Model, String> {
    let mut host = ScriptHost::with_options(Default::default());
    host.run(source, "revolve")
        .map_err(|e| e.to_string())?
        .model
        .ok_or_else(|| "no model".into())
}

fn volume(model: &mut popkern_script::Model) -> f64 {
    model
        .to_parts_all()
        .expect("parts")
        .iter()
        .map(|p| popkern_measure::body_volume(p.brep(), &p.brep().tolerance).expect("volume"))
        .sum()
}

/// A hub (r 10, z 0..10), a stop (within 25 of the axis between 120°
/// and 135°, built at 0°..15° with its −y face labelled, then turned
/// 120°), and a rib (5 ≤ x ≤ 20, 0 ≤ z ≤ 4) turned from 0° up to the
/// stop's face as a boss into the stop.
fn rib(end: &str, options: &str) -> String {
    let (c, s) = (7.5f64.to_radians(), 15f64.to_radians());
    format!(
        "let m = model();\nlet hub = m.profile(xy()).circle(10.0).sketch(\"hub\").extrude(blind(10.0));\nlet stop = m.profile(xy()).path(0.0, 0.0).line_to(25.0, 0.0).arc_to([{:?}, {:?}], [{:?}, {:?}]).close().sketch(\"stop\").extrude(blind(10.0));\nm.label_face(\"stop_face\", faces().facing_neg(\"y\"));\nstop = stop.rotate(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), 120.0);\nlet rib = m.profile(xz()).rect(15.0, 4.0).at(12.5, 2.0).sketch(\"rib\").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), {end}{options});\nhub = hub.union(stop);\nm\n",
        25.0 * c.cos(),
        25.0 * c.sin(),
        25.0 * s.cos(),
        25.0 * s.sin()
    )
}

#[test]
fn a_rib_turns_to_a_stop_on_another_body() {
    let source = rib(
        "up_to_surface(stop, faces().labeled(\"stop_face\"))",
        ", #{ mode: \"boss\", label: \"rib\" }",
    );
    let mut model = run(&source).unwrap_or_else(|e| panic!("{e}"));
    let v = volume(&mut model);
    // The hub, the stop outside it, the rib's 120° outside it.
    let exact = 1000.0 * PI + 15.0 / 360.0 * PI * 525.0 * 10.0 + PI * 300.0 * 4.0 / 3.0;
    assert!((v - exact).abs() < 1e-6 * exact, "{v} vs {exact}");
    let emitted = to_rhai(&mut model).expect("emit");
    assert!(emitted.contains("up_to_surface("), "{emitted}");
    let v2 = volume(&mut run(&emitted).unwrap_or_else(|e| panic!("{e}\n{emitted}")));
    assert!((v - v2).abs() <= 1e-9 * v, "emitted serves {v2}\n{emitted}");
}

#[test]
fn a_new_body_turned_to_a_slanted_face_declines() {
    // A stop face that the rib's points meet at different angles: the
    // plane x = 2 through y, met by the rib's inner corner long before
    // its outer one. Only a boss into its body could follow it.
    let source = "let m = model();\nlet wall = m.profile(xy()).rect(1.0, 60.0).at(-2.5, 0.0).sketch(\"wall\").extrude(blind(10.0));\nm.label_face(\"wall_face\", faces().facing(\"x\"));\nlet rib = m.profile(xz()).rect(15.0, 4.0).at(12.5, 2.0).sketch(\"rib\").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), up_to_surface(wall, faces().labeled(\"wall_face\")));\nm\n";
    let err = match run(source) {
        Ok(_) => panic!("a new body turned to a slanted face served"),
        Err(e) => e,
    };
    assert!(err.contains("only a boss"), "{err}");
}
