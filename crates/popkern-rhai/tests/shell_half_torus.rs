//! PK-1218: a half torus shelled open at its cut ends. The plane of each
//! end holds the torus's axis and meets the offset torus in a meridian
//! circle of the offset minor radius; the coaxial arm took planes square
//! to the axis only, and the inner rims kept the outer radius. Volume
//! π²R(r² − (r − t)²) = 21π²/2.

use popkern_core::tolerance::Tolerance;

#[test]
fn a_half_torus_shelled_open_at_its_ends_reads_its_closed_form() {
    let source = "let m = model();\n\
        let t = m.profile(xz()).circle(2.0).at(6.0, 0.0).sketch(\"tube\")\n\
            .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());\n\
        let half = m.tool(|m| m.profile(xy_at(-5.0)).rect(30.0, 15.0).at(0.0, -7.5)\n\
            .sketch(\"half\").extrude(blind(10.0)));\n\
        t = t.subtract(half);\n\
        t = t.shell(faces().kind(\"plane\").at_height(0.0), 0.5);\nm";
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(source, "half_torus")
        .unwrap_or_else(|e| panic!("{e}"));
    let body = run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone();
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    let v = popkern_measure::body_volume(&body, &Tolerance::default()).unwrap();
    let want = 10.5 * std::f64::consts::PI.powi(2);
    assert!((v - want).abs() < 1e-6 * want, "{v} want {want}");
}
