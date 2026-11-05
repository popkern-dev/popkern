//! An extrude up to a face of another body (brep-bench §4.4): the
//! end face is selected on the body named in the end, the boss merges
//! into that body, and the emitted script names it the same way.

use popkern_rhai::{ScriptHost, to_rhai};

fn run(source: &str) -> popkern_script::Model {
    let mut host = ScriptHost::with_options(Default::default());
    host.run(source, "up_to")
        .unwrap_or_else(|e| panic!("{e}\n{source}"))
        .model
        .expect("model")
}

fn volume(model: &mut popkern_script::Model) -> f64 {
    model
        .to_parts_all()
        .expect("parts")
        .iter()
        .map(|p| popkern_measure::body_volume(p.brep(), &p.brep().tolerance).expect("volume"))
        .sum()
}

/// The volume `source` serves, checked to survive emission.
fn served(source: &str) -> f64 {
    let mut model = run(source);
    let v = volume(&mut model);
    let emitted = to_rhai(&mut model).expect("emit");
    assert!(emitted.contains("up_to_surface(") || emitted.contains("offset_from_surface("));
    let v2 = volume(&mut run(&emitted));
    assert!(
        (v - v2).abs() <= 1e-9 * v,
        "emitted serves {v2}, the script {v}\n{emitted}"
    );
    v
}

const BLOCK: &str = "let m = model();\nlet block = m.profile(xy()).rect(20.0, 20.0).sketch(\"block\").extrude(blind(10.0));\n";

#[test]
fn a_pad_runs_down_to_another_bodys_face() {
    // A 6 × 6 pad sketched facing down at z = 15, run down to the
    // block's top at z = 10 and joined to it.
    let source = format!(
        "{BLOCK}let pad = m.profile(plane_at([0.0, 0.0, 15.0], [0.0, 0.0, -1.0])).rect(6.0, 6.0).sketch(\"pad\").extrude(up_to_surface(block, faces().facing(\"z\")), #{{ mode: \"boss\", label: \"pad\" }});\nm\n"
    );
    let v = served(&source);
    assert!((v - (4000.0 + 36.0 * 5.0)).abs() < 1e-6 * v, "{v}");
}

#[test]
fn a_pad_stops_short_of_another_bodys_face() {
    // The same pad stopped 2 short of the block: a separate body 3 tall.
    let source = format!(
        "{BLOCK}let pad = m.profile(plane_at([0.0, 0.0, 15.0], [0.0, 0.0, -1.0])).rect(6.0, 6.0).sketch(\"pad\").extrude(offset_from_surface(block, faces().facing(\"z\"), 2.0));\nm\n"
    );
    let v = served(&source);
    assert!((v - (4000.0 + 36.0 * 3.0)).abs() < 1e-6 * v, "{v}");
}

#[test]
fn the_end_face_is_selected_on_the_named_body() {
    // On the pad's own body there is no face facing +z to stop at.
    let source = format!(
        "{BLOCK}let pad = m.profile(plane_at([0.0, 0.0, 15.0], [0.0, 0.0, -1.0])).rect(6.0, 6.0).sketch(\"pad\").extrude(up_to_surface(faces().facing(\"z\").at_height(10.0)));\nm\n"
    );
    let mut host = ScriptHost::with_options(Default::default());
    let err = match host.run(&source, "own") {
        Ok(_) => panic!("the pad's own body has no face at z = 10"),
        Err(e) => e.to_string(),
    };
    assert!(
        err.contains("Selection") || err.contains("selection"),
        "{err}"
    );
}

const DRUM: &str = "let m = model();\nlet drum = m.profile(yz_at(-20.0)).circle(10.0).sketch(\"drum\").extrude(blind(40.0));\n";
const PAD: &str =
    "m.profile(plane_at([0.0, 0.0, 15.0], [0.0, 0.0, -1.0])).rect(6.0, 6.0).sketch(\"pad\")";

/// A pad run down to a drum's side ends ON it: the pad's volume above
/// the cylinder, 6 ∫ (15 − √(100 − y²)) dy over |y| ≤ 3, on the drum.
#[test]
fn a_pad_follows_a_drums_side() {
    let source = format!(
        "{DRUM}let pad = {PAD}.extrude(up_to_surface(drum, faces().kind(\"cylinder\")), #{{ mode: \"boss\", label: \"pad\" }});\nm\n"
    );
    let pi = std::f64::consts::PI;
    let under = 3.0 * 91f64.sqrt() + 100.0 * 0.3f64.asin();
    let exact = 4000.0 * pi + 6.0 * (90.0 - under);
    let v = served(&source);
    assert!((v - exact).abs() < 1e-6 * exact, "{v} vs {exact}");
}

/// A post up to a hollow dome's inner face (the sphere of radius 18)
/// ends on it: the post's volume is 2π ∫ ρ √(324 − ρ²) dρ to ρ = 3.
/// Under a concave face the centroid's ray is already the farthest, so
/// the old flat end served this too; the drum is the case that tells
/// the two apart.
#[test]
fn a_post_follows_a_domes_inner_face() {
    let source = "let m = model();\nfn half_ball(m, r) {\n  m.profile(xz()).path(0.0, 0.0).line_to(r, 0.0).arc_to([r * 0.6, r * 0.8], [0.0, r]).close().sketch(\"half_disc\").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())\n}\nlet dome = half_ball(m, 20.0);\nlet cavity = m.tool(|m| half_ball(m, 18.0));\nm.set_active(cavity);\nm.label_face(\"inner\", faces().kind(\"sphere\"));\ndome = dome.subtract(cavity);\nlet plate = m.tool(|m| m.profile(xy_at(-2.0)).circle(20.0).sketch(\"plate\").extrude(blind(2.0)));\ndome = dome.union(plate);\nlet post = m.profile(xy()).circle(3.0).sketch(\"post\").extrude(up_to_surface(dome, faces().labeled(\"inner\")), #{ mode: \"boss\", label: \"post\" });\nm\n";
    let pi = std::f64::consts::PI;
    let post = 2.0 * pi / 3.0 * (5832.0 - 315.0 * 315f64.sqrt());
    let exact = 800.0 * pi + 2.0 / 3.0 * pi * (8000.0 - 5832.0) + post;
    let v = served(source);
    assert!((v - exact).abs() < 1e-6 * exact, "{v} vs {exact}");
}

fn declined(source: &str) -> String {
    let mut host = ScriptHost::with_options(Default::default());
    match host.run(source, "declined") {
        Ok(_) => panic!("served\n{source}"),
        Err(e) => e.to_string(),
    }
}

#[test]
fn a_new_body_up_to_a_curved_face_declines() {
    // Only a boss into the drum can end on its side; a separate body
    // would need trimming by the cylinder (it used to end flat).
    let err = declined(&format!(
        "{DRUM}let pad = {PAD}.extrude(up_to_surface(drum, faces().kind(\"cylinder\")));\nm\n"
    ));
    assert!(err.contains("only a boss"), "{err}");
}

#[test]
fn a_boss_through_a_thin_wall_declines() {
    // A tube 0.1 thick: the pad's end spans 10 − √91 ≈ 0.46 of depth,
    // more than the wall behind the face.
    let err = declined(&format!(
        "{DRUM}m.label_face(\"skin\", faces().kind(\"cylinder\"));\nlet bore = m.tool(|m| m.profile(yz_at(-21.0)).circle(9.9).sketch(\"bore\").extrude(blind(42.0)));\ndrum = drum.subtract(bore);\nlet pad = {PAD}.extrude(up_to_surface(drum, faces().labeled(\"skin\")), #{{ mode: \"boss\" }});\nm\n"
    ));
    assert!(err.contains("thinner behind it"), "{err}");
}
