//! Primitives (PK-1083): each serves its closed-form volume (and, for
//! the analytic faces, area), replays from its emitted script, stands
//! where its placement puts it, and names its faces by role.

use popkern_rhai::{ScriptHost, to_rhai};
use std::f64::consts::PI;

fn run(source: &str) -> Result<popkern_script::Model, String> {
    let mut host = ScriptHost::with_options(Default::default());
    host.run(source, "primitives")
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

/// `body` serves `volume` (and `area` when given) to 1e-9, and so does
/// its emitted script.
fn assert_serves(body: &str, volume: f64, area: Option<f64>) {
    let source = format!("let m = model();\nlet b = {body};\nm\n");
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
    let mut again = run(&emitted).unwrap_or_else(|e| panic!("{e}\n{emitted}"));
    let (v2, a2) = measure(&mut again);
    assert!(
        (v - v2).abs() <= 1e-12 * v && (a - a2).abs() <= 1e-12 * a,
        "{emitted}"
    );
}

#[test]
fn boxes_and_cylinders() {
    assert_serves("m.box(2.0, 3.0, 4.0)", 24.0, Some(52.0));
    assert_serves("m.box(2.0, 3.0, 4.0, #{ centred: true })", 24.0, Some(52.0));
    assert_serves(
        "m.cylinder(2.0, 5.0, #{ label: \"post\" })",
        20.0 * PI,
        Some(28.0 * PI),
    );
}

#[test]
fn cones_and_frustums() {
    assert_serves("m.cone(3.0, 1.0, 6.0)", 26.0 * PI, None);
    assert_serves("m.cone(3.0, 0.0, 6.0)", 18.0 * PI, None);
    assert_serves("m.cone(0.0, 3.0, 6.0)", 18.0 * PI, None);
}

#[test]
fn spheres_and_zones() {
    assert_serves("m.sphere(3.0)", 36.0 * PI, Some(36.0 * PI));
    // Between −30° and 60° on r 2: π ∫ (4 − z²) dz from −1 to √3.
    let zone = PI * (3.0 * 3f64.sqrt() + 11.0 / 3.0);
    assert_serves("m.sphere(2.0, #{ lat: [-30.0, 60.0] })", zone, None);
}

#[test]
fn tori_and_sectors() {
    assert_serves("m.torus(5.0, 1.0)", 10.0 * PI * PI, Some(20.0 * PI * PI));
    assert_serves("m.torus(5.0, 1.0, #{ angle: 90.0 })", 2.5 * PI * PI, None);
}

/// A box placed at (10, 0, 0) with its axis along +y and its x along
/// +z: local x (2) runs up z, local y (3) along +x, local z (4) along y.
#[test]
fn a_placed_box_stands_in_its_frame() {
    let source = "let m = model();\nlet b = m.box(2.0, 3.0, 4.0, #{ at: [10.0, 0.0, 0.0], axis: [0.0, 1.0, 0.0], ref: [0.0, 0.0, 1.0] });\nm\n";
    let mut model = run(source).expect("serves");
    let parts = model.to_parts_all().expect("parts");
    let brep = parts[0].brep();
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for v in brep.vertices.values() {
        let p = [v.position.x, v.position.y, v.position.z];
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let expect = ([10.0, 0.0, 0.0], [13.0, 4.0, 2.0]);
    for k in 0..3 {
        assert!(
            (lo[k] - expect.0[k]).abs() < 1e-9 && (hi[k] - expect.1[k]).abs() < 1e-9,
            "{lo:?} {hi:?}"
        );
    }
}

/// The top of a cylinder is its `Cap 1`, by name, whatever it is fed to.
#[test]
fn faces_are_named_by_role() {
    let source = "let m = model();\nlet b = m.cylinder(2.0, 5.0, #{ label: \"post\" });\nb = b.shell(faces().named(\"post/Cap/1\"), 0.5);\nm\n";
    let mut model = run(source).unwrap_or_else(|e| panic!("{e}"));
    let (v, _) = measure(&mut model);
    // Open at the top, 0.5 thick: the cup is the cylinder less the bore
    // of radius 1.5 down to z = 0.5.
    let exact = 20.0 * PI - PI * 2.25 * 4.5;
    assert!((v - exact).abs() <= 1e-9 * exact, "{v} vs {exact}");
}

#[test]
fn bad_dimensions_are_named() {
    for (body, needle) in [
        ("m.torus(1.0, 1.0)", "minor"),
        ("m.cone(0.0, 0.0, 1.0)", "both be 0"),
        ("m.sphere(1.0, #{ lat: [30.0, 10.0] })", "lo < hi"),
        ("m.box(1.0, 0.0, 1.0)", "positive"),
    ] {
        let source = format!("let m = model();\nlet b = {body};\nm\n");
        match run(&source) {
            Ok(_) => panic!("{body} served"),
            Err(e) => assert!(e.contains(needle), "{body}: {e}"),
        }
    }
}

/// A `plane_at` origin is recorded in SI and spelled back in the
/// model's unit (a millimetre script used to come back a thousandth).
#[test]
fn a_custom_plane_round_trips_in_millimetres() {
    let source = "let m = model(#{ unit: \"mm\" });\nlet b = m.profile(plane_at([0.0, 0.0, 15.0], [0.0, 0.0, 1.0])).rect(10.0, 10.0).sketch(\"s\").extrude(blind(5.0));\nlet c = m.box(10.0, 10.0, 15.0, #{ at: [-5.0, -5.0, 0.0] });\nm\n";
    let mut model = run(source).expect("serves");
    let (v, _) = measure(&mut model);
    let emitted = to_rhai(&mut model).expect("emit");
    let (v2, _) = measure(&mut run(&emitted).unwrap_or_else(|e| panic!("{e}\n{emitted}")));
    assert!((v - v2).abs() <= 1e-12 * v, "{v} vs {v2}\n{emitted}");
}

/// A millimetre document emits in millimetres: its header used to say
/// metres over millimetre literals, a script 1000 times too large.
#[test]
fn a_millimetre_model_round_trips_at_its_scale() {
    let source = "let m = model(#{ unit: \"mm\" });\nlet b = m.profile(xy()).rect(10.0, 10.0).sketch(\"s\").extrude(blind(5.0));\nlet c = m.tool(|m| m.profile(xz()).rect(2.0, 2.0).at(20.0, 1.0).sketch(\"r\").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle()));\nm\n";
    let mut model = run(source).expect("serves");
    let (v, _) = measure(&mut model);
    let emitted = to_rhai(&mut model).expect("emit");
    assert!(
        emitted.starts_with("let m = model(#{ unit: \"mm\" });"),
        "{emitted}"
    );
    let (v2, _) = measure(&mut run(&emitted).unwrap_or_else(|e| panic!("{e}\n{emitted}")));
    assert!((v - v2).abs() <= 1e-12 * v, "{v} vs {v2}\n{emitted}");
}
