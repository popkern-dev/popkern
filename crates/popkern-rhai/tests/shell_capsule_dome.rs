//! PK-1216: a capsule shelled open at its dome. The rim vertex, seeded at
//! the old rim where the offset wall and the dome's sphere are both
//! radial, is reseeded inside the opening; the rim circle's Newton is
//! seeded at its recomputed end, off the dome's flat station; and the
//! wall's vertical on the sphere is its great-circle arc. Volume
//! 1000π/3 − (16π·12 + π·2²(3·5 − 2)/3) = 124π.

use popkern_core::tolerance::Tolerance;

#[test]
fn a_capsule_shelled_open_at_its_dome_reads_its_closed_form() {
    let source = "let m = model();\n\
        let c = m.profile(xz()).path(0.0, 0.0).line_to(5.0, 0.0).line_to(5.0, 10.0)\n\
            .arc_to([5.0 * cos(PI() / 4.0), 10.0 + 5.0 * sin(PI() / 4.0)], [0.0, 15.0]).close()\n\
            .sketch(\"capsule\").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());\n\
        c = c.shell(faces().kind(\"sphere\"), 1.0);\nm";
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(source, "capsule")
        .unwrap_or_else(|e| panic!("{e}"));
    let body = run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone();
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    let v = popkern_measure::body_volume(&body, &Tolerance::default()).unwrap();
    let want = 124.0 * std::f64::consts::PI;
    assert!((v - want).abs() < 1e-6 * want, "{v} want {want}");
}
