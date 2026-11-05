//! A cut that removes no material declines (brep-bench §4.5): an
//! `on_face` cut placed off its face, or extruded away from the
//! material (a `blind` end runs along the face's outward normal), used
//! to serve the body unchanged. The controls drill the analytic hole
//! and shear a boss off flush, a real cut whose every new face lies on
//! a face of the body.

use popkern_rhai::ScriptHost;

const BOX: &str = "let m = model();\nlet b = m.profile(xy()).rect(10.0, 8.0).sketch(\"box\").extrude(blind(6.0));\n";

fn served_volume(source: &str) -> Result<f64, String> {
    let mut host = ScriptHost::with_options(Default::default());
    let mut model = host
        .run(source, "cut")
        .map_err(|e| e.to_string())?
        .model
        .ok_or("no model")?;
    let parts = model.to_parts_all().map_err(|e| e.to_string())?;
    Ok(parts
        .iter()
        .map(|p| popkern_measure::body_volume(p.brep(), &p.brep().tolerance).expect("volume"))
        .sum())
}

/// A radius-1 cut sketched on the box's top face at `(x, y)`, its
/// blind end `depth` along the face's outward normal.
fn drill(x: f64, y: f64, depth: f64) -> Result<f64, String> {
    served_volume(&format!(
        "{BOX}let h = m.profile(on_face(faces().facing(\"z\"))).circle(1.0).at({x:?}, {y:?}).sketch(\"hole\").extrude(blind({depth:?}), #{{ mode: \"cut\", label: \"drill\" }});\nm\n"
    ))
}

fn assert_declined(result: Result<f64, String>, case: &str) {
    match result {
        Ok(v) => panic!("{case}: served {v}"),
        Err(e) => assert!(e.contains("removes no material"), "{case}: {e}"),
    }
}

#[test]
fn a_cut_into_the_face_drills_the_hole() {
    let v = drill(0.0, 0.0, -2.0).expect("the cut serves");
    let exact = 480.0 - 2.0 * std::f64::consts::PI;
    assert!((v - exact).abs() < 1e-6 * exact, "{v} vs {exact}");
}

#[test]
fn a_cut_off_the_face_says_it_removes_nothing() {
    assert_declined(drill(100.0, 100.0, -2.0), "beside the body");
}

#[test]
fn a_cut_away_from_the_material_says_it_removes_nothing() {
    // The tool stands on the top face: its disc imprints there, flush.
    assert_declined(drill(0.0, 0.0, 2.0), "away from the material");
}

#[test]
fn a_boss_sheared_off_flush_is_a_cut() {
    // A 2 × 2 × 1 boss on the top, then a 4 × 4 cut through the boss's
    // height: the floor lies on the box's top and the cut's top on the
    // boss's, so only the volume tells this from a no-op.
    let source = format!(
        "{BOX}let s = m.profile(on_face(faces().facing(\"z\"))).rect(2.0, 2.0).sketch(\"stud\").extrude(blind(1.0), #{{ mode: \"boss\" }});\nlet c = m.profile(on_face(faces().facing(\"z\").at_height(7.0))).rect(4.0, 4.0).sketch(\"shear\").extrude(blind(-1.0), #{{ mode: \"cut\" }});\nm\n"
    );
    let v = served_volume(&source).expect("the shear serves");
    assert!((v - 480.0).abs() < 1e-6 * 480.0, "{v}");
}

/// A cut sketched on a free plane names no host: it cuts the body it
/// reaches. It used to be kept as a body of its own, a solid cylinder
/// beside the box (480 + 2 pi) instead of the hole.
#[test]
fn a_cut_on_a_free_plane_cuts_the_body_it_reaches() {
    let v = served_volume(&format!(
        "{BOX}let h = m.profile(xy_at(6.0)).circle(1.0).sketch(\"hole\").extrude(blind(-2.0), #{{ mode: \"cut\" }});\nm\n"
    ))
    .expect("the cut serves");
    let exact = 480.0 - 2.0 * std::f64::consts::PI;
    assert!((v - exact).abs() < 1e-6 * exact, "{v} vs {exact}");
}

#[test]
fn a_cut_on_a_free_plane_beside_the_body_removes_nothing() {
    assert_declined(
        served_volume(&format!(
            "{BOX}let h = m.profile(xy_at(6.0)).circle(1.0).at(40.0, 0.0).sketch(\"hole\").extrude(blind(-2.0), #{{ mode: \"cut\" }});\nm\n"
        )),
        "free plane, beside the body",
    );
}

/// Two bodies in reach: the cut does not choose between them.
#[test]
fn a_cut_reaching_two_bodies_names_the_ambiguity() {
    let source = format!(
        "{BOX}let c = m.profile(xy()).rect(10.0, 8.0).at(12.0, 0.0).sketch(\"other\").extrude(blind(6.0));\nlet h = m.profile(xy_at(6.0)).rect(20.0, 1.0).at(6.0, 0.0).sketch(\"slot\").extrude(blind(-2.0), #{{ mode: \"cut\" }});\nm\n"
    );
    match served_volume(&source) {
        Ok(v) => panic!("served {v}"),
        Err(e) => assert!(e.contains("reaches 2 bodies"), "{e}"),
    }
}
