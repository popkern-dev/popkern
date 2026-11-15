//! PK-1197: a loft between two convex sections of different side counts
//! in parallel planes is their convex hull: flat triangles and
//! trapezoids where a side faces a corner or a parallel side, oblique
//! cone pieces where a corner faces an arc. The skinned surface paired
//! the sections by chord length and read under half the hull.

use popkern_core::tolerance::Tolerance;

fn volume(source: &str) -> f64 {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host.run(source, "loft").unwrap_or_else(|e| panic!("{e}"));
    let body = run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone();
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    popkern_measure::body_volume(&body, &Tolerance::default()).unwrap()
}

#[test]
fn a_square_to_round_transition_is_its_hull() {
    let v = volume(
        "let m = model();\n\
         let square = m.profile(xy()).rect(40.0, 40.0).sketch(\"square\");\n\
         let round = m.profile(xy_at(30.0)).circle(15.0).sketch(\"round\");\n\
         m.loft([square, round]);\nm",
    );
    let want = 28000.0 + 2250.0 * std::f64::consts::PI;
    assert!((v - want).abs() < 1e-6 * want, "{v} want {want}");
}

#[test]
fn a_square_to_octagon_hood_is_its_hull() {
    let v = volume(
        "let m = model();\n\
         let square = m.profile(xy()).rect(40.0, 40.0).sketch(\"square\");\n\
         let octagon = m.profile(xy_at(30.0)).polygon([[6.0, -14.0], [14.0, -6.0], [14.0, 6.0], \
         [6.0, 14.0], [-6.0, 14.0], [-14.0, 6.0], [-14.0, -6.0], [-6.0, -14.0]]).sketch(\"octagon\");\n\
         m.loft([square, octagon]);\nm",
    );
    assert!((v - 33760.0).abs() < 1e-6 * 33760.0, "{v}");
}
