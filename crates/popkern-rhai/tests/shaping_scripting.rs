//! `offset_solid`, `thicken` and `draft` from Rhai (spec §5 preview):
//! each recorded feature serves the analytic volume, and the document
//! emits back to a script that replays to the same volume.

use popkern_rhai::{ScriptHost, to_rhai};

fn volume(model: &mut popkern_script::Model) -> f64 {
    model
        .to_parts_all()
        .expect("parts")
        .iter()
        .map(|p| popkern_measure::body_volume(p.brep(), &p.brep().tolerance).expect("volume"))
        .sum()
}

/// The volume `source` builds, and the volume its emitted document
/// rebuilds (the round trip must agree).
fn volume_and_round_trip(source: &str, name: &str) -> f64 {
    let mut host = ScriptHost::with_options(Default::default());
    let mut model = host
        .run(source, name)
        .unwrap_or_else(|e| panic!("{name}: {e}"))
        .model
        .expect("model");
    let v = volume(&mut model);
    let emitted = to_rhai(&mut model).unwrap_or_else(|e| panic!("{name}: emit: {e}"));
    let mut again = host
        .run(&emitted, name)
        .unwrap_or_else(|e| panic!("{name}: emitted script: {e}\n{emitted}"))
        .model
        .expect("model");
    let v2 = volume(&mut again);
    assert!(
        (v - v2).abs() <= 1e-9 * v.abs().max(1.0),
        "{name}: emitted document serves {v2}, the script {v}\n{emitted}"
    );
    v
}

fn assert_near(name: &str, v: f64, exact: f64) {
    assert!(
        (v - exact).abs() <= 1e-6 * exact.abs(),
        "{name}: {v} against {exact}"
    );
}

const BOX: &str = "let m = model();\nlet b = m.profile(xy()).rect(10.0, 8.0).sketch(\"box\").extrude(blind(6.0));\n";

#[test]
fn a_box_offset_grows_and_shrinks_by_the_distance() {
    // Sharp corners: the box grown to 12 × 10 × 8 and shrunk to 8 × 6 × 4.
    let grown = format!("{BOX}b = b.offset_solid(1.0, #{{ join: \"intersection\" }});\nm\n");
    assert_near("grown", volume_and_round_trip(&grown, "grown"), 960.0);
    let shrunk = format!("{BOX}b = b.offset_solid(-1.0, #{{ join: \"intersection\" }});\nm\n");
    assert_near("shrunk", volume_and_round_trip(&shrunk, "shrunk"), 192.0);
    // Rounded (the default): Steiner's formula, V + A r + (π/4) L r² +
    // (4π/3) r³ with the box's area 376 and edge length 96.
    let rounded = format!("{BOX}b = b.offset_solid(1.0, #{{ label: \"grow\" }});\nm\n");
    let pi = std::f64::consts::PI;
    let exact = 480.0 + 376.0 + pi / 4.0 * 96.0 + 4.0 / 3.0 * pi;
    assert_near("rounded", volume_and_round_trip(&rounded, "rounded"), exact);
}

#[test]
fn a_box_thickened_is_a_wall_half_each_side() {
    // A closed box thickened by 1 is the 11 × 9 × 7 box less the
    // 9 × 7 × 5 one.
    let source = format!("{BOX}b = b.thicken(1.0, #{{ join: \"intersection\" }});\nm\n");
    assert_near(
        "thickened",
        volume_and_round_trip(&source, "thick"),
        693.0 - 315.0,
    );
}

#[test]
fn a_face_drafted_about_the_base_leans_in() {
    // The +x side drafted 5° about the base (pull +z) leans in by
    // z·tan 5°: a wedge ½ · 6 · 6 tan 5° · 8 comes off.
    let source = format!(
        "{BOX}b = b.draft(faces().facing(\"x\"), 5.0, xy(), #{{ label: \"taper\" }});\nm\n"
    );
    let t = 5.0f64.to_radians().tan();
    assert_near(
        "drafted",
        volume_and_round_trip(&source, "draft"),
        480.0 - 144.0 * t,
    );
}

#[test]
fn a_face_drafted_about_a_face_is_the_same_draft() {
    // The base face as the neutral plane: its outward normal is -z, so
    // the pull is down and the side leans out going up instead.
    let source =
        format!("{BOX}b = b.draft(faces().facing(\"x\"), 5.0, faces().facing_neg(\"z\"));\nm\n");
    let t = 5.0f64.to_radians().tan();
    let v = volume_and_round_trip(&source, "draft_face");
    assert_near("drafted about the base face", v, 480.0 + 144.0 * t);
}

#[test]
fn a_draft_names_a_bad_neutral_plane() {
    let source = format!("{BOX}b = b.draft(faces().facing(\"x\"), 5.0, 3.0);\nm\n");
    let mut host = ScriptHost::with_options(Default::default());
    let err = match host.run(&source, "bad") {
        Ok(_) => panic!("a number is no neutral plane"),
        Err(e) => e.to_string(),
    };
    assert!(err.contains("neutral"), "{err}");
}

/// The box opened at its top is a five-face sheet; thickened by 1 with
/// sharp joins it is a tray whose walls stand outside the old box,
/// inside it, or half each way.
#[test]
fn an_opened_box_thickens_to_the_side_asked() {
    let opened =
        format!("{BOX}b = b.remove_faces(faces().facing(\"z\"), #{{ label: \"open\" }});\n");
    for (side, exact) in [
        ("outward", 12.0 * 10.0 * 7.0 - 480.0),
        ("inward", 480.0 - 240.0),
        ("both", 11.0 * 9.0 * 6.5 - 9.0 * 7.0 * 5.5),
    ] {
        let source = format!(
            "{opened}b = b.thicken(1.0, #{{ side: \"{side}\", join: \"intersection\" }});\nm\n"
        );
        assert_near(side, volume_and_round_trip(&source, side), exact);
    }
}

#[test]
fn a_bad_side_is_named() {
    let source = format!(
        "{BOX}b = b.remove_faces(faces().facing(\"z\"));\nb = b.thicken(1.0, #{{ side: \"up\" }});\nm\n"
    );
    let mut host = ScriptHost::with_options(Default::default());
    let err = match host.run(&source, "bad_side") {
        Ok(_) => panic!("\"up\" is no side"),
        Err(e) => e.to_string(),
    };
    assert!(err.contains("side"), "{err}");
}
