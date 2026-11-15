//! A horn torus (a circle revolved about a line it touches) is a valid
//! body: its chart normal vanishes only at the point on the axis, which
//! the orientation check sampled alone at the chart's middle.

use popkern_core::tolerance::Tolerance;

#[test]
fn a_horn_torus_validates_and_reads_its_volume() {
    let source = "let m = model();\n\
                  m.profile(xz()).circle(5.0).at(5.0, 0.0).sketch(\"section\")\n\
                      .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());\nm";
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host.run(source, "horn").unwrap_or_else(|e| panic!("{e}"));
    let body = run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone();
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    let v = popkern_measure::body_volume(&body, &Tolerance::default()).unwrap();
    let want = 250.0 * std::f64::consts::PI.powi(2);
    assert!((v - want).abs() < 1e-6 * want, "{v} want {want}");
}
