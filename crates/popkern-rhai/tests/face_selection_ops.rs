//! Face queries joined and negated (brep-bench §4.1, §4.2): `or`,
//! `except` and `not` on `faces()`, and `facing` on flat faces bounded
//! by circles. The shells are the bench's cases, against their closed
//! forms (box 10 × 8 × 6, cylinder r 5 × 10, walls 1 thick).

use popkern_rhai::{ScriptHost, to_rhai};

fn run(source: &str) -> popkern_script::Model {
    let mut host = ScriptHost::with_options(Default::default());
    host.run(source, "faces")
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

/// The volume a box shelled 1 thick open at `opening` serves, checked
/// to survive emission.
fn shelled_box(opening: &str) -> f64 {
    let source = format!(
        "let m = model();\nlet b = m.profile(xy()).rect(10.0, 8.0).sketch(\"box\").extrude(blind(6.0));\nb = b.shell({opening}, 1.0);\nm\n"
    );
    let mut model = run(&source);
    let v = volume(&mut model);
    let emitted = to_rhai(&mut model).expect("emit");
    let mut again = run(&emitted);
    assert!((volume(&mut again) - v).abs() < 1e-9 * v, "{emitted}");
    v
}

#[test]
fn two_openings_join_with_or() {
    let v = shelled_box(r#"faces().facing("z").or(faces().facing_neg("z"))"#);
    assert!((v - 192.0).abs() < 1e-6 * 192.0, "{v}");
    let v = shelled_box(r#"faces().facing("z").or(faces().facing_neg("y"))"#);
    assert!((v - 200.0).abs() < 1e-6 * 200.0, "{v}");
}

#[test]
fn every_face_but_one_is_except_or_not() {
    let v = shelled_box(r#"faces().except(faces().facing_neg("z"))"#);
    assert!((v - 80.0).abs() < 1e-6 * 80.0, "{v}");
    let v = shelled_box(r#"faces().facing_neg("z").not()"#);
    assert!((v - 80.0).abs() < 1e-6 * 80.0, "{v}");
}

#[test]
fn a_query_spells_and_parses_back() {
    let text = r#"faces().facing(z).or(faces().except(faces().facing_neg(z)))"#;
    let sel = popkern_script::selection_parse::parse_faces(text).expect("parses");
    assert_eq!(sel.describe(), text);
}

/// A cylinder's end faces have no vertex polygon (one closed circle
/// each): `facing` reads their normal from the circle.
#[test]
fn a_circular_end_face_faces_its_axis() {
    let source = "let m = model();\nlet c = m.profile(xy()).circle(5.0).sketch(\"cylinder\").extrude(blind(10.0));\nc = c.shell(faces().facing(\"z\"), 1.0);\nm\n";
    let mut model = run(source);
    let v = volume(&mut model);
    let exact = 106.0 * std::f64::consts::PI;
    assert!((v - exact).abs() < 1e-6 * exact, "{v} vs {exact}");
}
