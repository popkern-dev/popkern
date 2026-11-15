//! PK-1217: a corner where three bevels meet closes with a flat face by
//! default, through the points where the setback lines meet on each
//! original face (a cube of side 10 chamfered 2 on every edge: the
//! triangles x + y + z <= 26, volume a^3 - 6ad^2 + 16d^3/3 = 2408/3);
//! `corner: "apex"` runs the bevel planes on until they meet (808).

use popkern_core::tolerance::Tolerance;

fn cube(options: &str) -> f64 {
    let source = format!(
        "let m = model();\n\
         let c = m.profile(xy()).rect(10.0, 10.0).at(5.0, 5.0).sketch(\"cube\").extrude(blind(10.0));\n\
         c = c.chamfer(edges(), 2.0{options});\nm"
    );
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host.run(&source, "cube").unwrap_or_else(|e| panic!("{e}"));
    let body = run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone();
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    assert_eq!(
        body.faces.len(),
        if options.contains("apex") { 18 } else { 26 },
        "{options}"
    );
    popkern_measure::body_volume(&body, &Tolerance::default()).unwrap()
}

#[test]
fn a_cube_chamfered_on_every_edge_closes_its_corners_flat() {
    let v = cube("");
    assert!((v - 2408.0 / 3.0).abs() < 1e-9 * v, "{v}");
}

#[test]
fn apex_corners_run_the_bevel_planes_on_until_they_meet() {
    let v = cube(", #{ corner: \"apex\" }");
    assert!((v - 808.0).abs() < 1e-9 * v, "{v}");
}

#[test]
fn an_unknown_corner_is_refused() {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let Err(e) = host.run(
        "let m = model();\nlet c = m.box(10.0, 10.0, 10.0);\nc = c.chamfer(edges(), 1.0, #{ corner: \"round\" });\nm",
        "cube",
    ) else {
        panic!("an unknown corner must be refused");
    };
    assert!(e.to_string().contains("flat | apex"), "{e}");
}
