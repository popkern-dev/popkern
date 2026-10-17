//! PK-460 fastener scripting surface: the Rhai verbs record real
//! fastener features that replay through `popkern-fastener` and emit
//! back to script faithfully.
//!
//! Only the classes the thread-boolean campaign has closed run
//! un-ignored (threaded rod, right-hand single-start); the wider
//! matrix (bolt heads, nuts, holes) unlocks with the walk-class
//! campaign, mirroring `popkern-part`'s ignored suites.

use popkern_rhai::{ScriptHost, to_rhai};

fn run(script: &str) -> popkern_rhai::RunOutcome {
    let mut host = ScriptHost::new();
    host.run(script, "fastener_test.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload))
}

#[test]
fn threaded_rod_records_serves_and_replays() {
    let outcome = run(r#"
let m = model();
let rod = m.threaded_rod("M6", 4.0, #{ label: "rod", tip_chamfer: false });
m.assert_watertight();
m
"#);
    let mut model = outcome.model.expect("model");
    let part = model.to_part().expect("part");
    // Designations bind millimetre-magnitude numbers in model units:
    // an M6 rod's crest edges sit at radius 3.0, its tip at z = 4.0
    // (vertex hulls in x/y under-read the crest — vertices sit at
    // scattered azimuths on the helical edges, so measure radius).
    let (max_r, max_z) = part
        .brep()
        .vertices
        .iter()
        .fold((0.0_f64, f64::MIN), |acc, (_, v)| {
            let p = v.position;
            (acc.0.max((p.x * p.x + p.y * p.y).sqrt()), acc.1.max(p.z))
        });
    assert!((max_r - 3.0).abs() < 1e-3, "major radius 3 units: {max_r}");
    assert!((max_z - 4.0).abs() < 1e-6, "length 4 units: {max_z}");
}

#[test]
fn spec_map_form_binds_model_units() {
    // Dimensionless model: #{ major, pitch } are raw model units.
    let outcome = run(r#"
let m = model();
let rod = m.threaded_rod(#{ major: 8.0, pitch: 1.25 }, 5.0, #{ label: "rod", tip_chamfer: false });
m.assert_watertight();
m
"#);
    let mut model = outcome.model.expect("model");
    let part = model.to_part().expect("part");
    let max_r = part
        .brep()
        .vertices
        .iter()
        .map(|(_, v)| (v.position.x * v.position.x + v.position.y * v.position.y).sqrt())
        .fold(0.0_f64, f64::max);
    assert!((max_r - 4.0).abs() < 1e-3, "major radius 4 units: {max_r}");
}

#[test]
fn threaded_rod_emits_back_to_script() {
    let outcome = run(r#"
let m = model();
let rod = m.threaded_rod("M6", 4.0, #{ label: "rod", tip_chamfer: false });
m
"#);
    let mut model = outcome.model.expect("model");
    let emitted = to_rhai(&mut model).expect("emit");
    assert!(
        emitted.contains("m.threaded_rod(#{ major:"),
        "emits the spec-map form: {emitted}"
    );
}

#[test]
fn unknown_designation_reports_the_grammars() {
    let mut host = ScriptHost::new();
    let err = host
        .run(
            r#"
let m = model();
m.threaded_rod("M7", 4.0);
m
"#,
            "bad_spec.rhai",
        )
        .expect_err("M7 is off-table");
    let text = err.to_string();
    assert!(
        text.contains("M7x"),
        "error suggests the explicit-pitch form: {text}"
    );
}

#[test]
fn bolt_nut_and_tapped_hole_serve() {
    run(r#"
let m = model();
let bolt = m.bolt("M6", 8.0, #{ thread_length: 5.0, label: "bolt" });
let nut = m.nut("M6", #{ clearance: 0.2, label: "nut" });
let plate = m.profile(xy()).rect(14.0, 14.0).sketch("plate").extrude(blind(3.5));
plate = plate.tapped_hole(axis_at([0.0, 0.0, 3.5], [0.0, 0.0, -1.0]), "M6", through_all(), #{ clearance: 0.2 });
m.assert_watertight();
m
"#);
}

/// PK-686: the nut verb exposes the part builder's chamfer toggles,
/// records them, and emits them back to script.
#[test]
fn nut_chamfer_toggles_serve_and_emit() {
    let outcome = run(r#"
let m = model();
let n = m.nut("M6", #{ label: "nut", corner_chamfer: false, mouth_chamfer: false });
m.assert_watertight();
m
"#);
    let mut model = outcome.model.expect("model");
    let emitted = to_rhai(&mut model).expect("emit");
    assert!(
        emitted.contains("corner_chamfer: false") && emitted.contains("mouth_chamfer: false"),
        "emits both toggles: {emitted}"
    );
}

/// PK-847: `cosmetic_thread: true` builds the threaded length as a plain
/// cylinder at the external root diameter d3 = d − 2·(17/24)H (no
/// helix), records the toggle, and emits it back to script.
#[test]
fn cosmetic_threads_serve_at_the_root_diameter_and_emit() {
    use popkern_core::geometry::SurfaceGeometry;
    let h = 0.8 * 3f64.sqrt() / 2.0; // M5: H = (√3/2)·P
    let d3 = 5.0 - 2.0 * (17.0 / 24.0) * h;
    let volume = |model: &mut popkern_script::Model| {
        let part = model.to_part().expect("part");
        let brep = part.brep();
        assert!(brep.validate_all().is_empty(), "{:?}", brep.validate_all());
        let helical = brep.faces.values().any(|f| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Nurbs(_))
            )
        });
        assert!(!helical, "a cosmetic thread cut a helix");
        let root = brep.faces.values().any(|f| {
            matches!(f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Cylinder { radius, .. }) if (radius - d3 / 2.0).abs() < 1e-9)
        });
        assert!(root, "no root-diameter cylinder");
        popkern_measure::body_volume(brep, &brep.tolerance).expect("volume")
    };
    // A rod is the root cylinder exactly.
    let mut rod = run(r#"
let m = model();
m.threaded_rod("M5", 10.0, #{ label: "rod", cosmetic_thread: true });
m.assert_watertight();
m
"#)
    .model
    .expect("model");
    let v = volume(&mut rod);
    let exact = std::f64::consts::PI * (d3 / 2.0).powi(2) * 10.0;
    assert!(((v - exact) / exact).abs() < 1e-9, "rod {v} vs {exact}");
    // A socket screw: head, plain shank, root-diameter thread; it emits
    // the toggle and the emitted script rebuilds the same body.
    let mut screw = run(r#"
let m = model();
m.bolt("M5", 12.0, #{ head: "socket", thread_length: 8.0, cosmetic_thread: true, label: "screw" });
m.assert_watertight();
m
"#)
    .model
    .expect("model");
    let v = volume(&mut screw);
    let emitted = to_rhai(&mut screw).expect("emit");
    assert!(
        emitted.contains("cosmetic_thread: true"),
        "emits the toggle: {emitted}"
    );
    let mut again = run(&emitted).model.expect("model");
    let w = volume(&mut again);
    assert!((v - w).abs() <= 1e-12 * v, "round trip {v} vs {w}");
}
