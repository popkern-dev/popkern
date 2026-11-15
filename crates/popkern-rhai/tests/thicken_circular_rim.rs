//! A sheet whose free rims trace circles as NURBS (a half-pipe made by
//! removing an extrude's flat faces) thickens both ways: each end arc's
//! wall is the annulus between the arc and its offset. Read as no circle,
//! the arc's quad collapsed (its ends and their offsets are collinear)
//! and the thicken declined.

use popkern_core::tolerance::Tolerance;

fn volume(side: &str) -> f64 {
    let source = format!(
        "let m = model();\n\
         let p = m.profile(xy()).path(5.0, 0.0).arc_to([0.0, 5.0], [-5.0, 0.0]).line_to(5.0, 0.0).close()\n\
             .sketch(\"half_disc\").extrude(blind(10.0));\n\
         p = p.remove_faces(faces().kind(\"plane\"));\n\
         p = p.thicken(1.0, #{{ side: \"{side}\" }});\nm"
    );
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(&source, "half_pipe")
        .unwrap_or_else(|e| panic!("{side}: {e}"));
    let b = run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone();
    let errors = b.validate_all();
    assert!(errors.is_empty(), "{side}: {errors:?}");
    popkern_measure::body_volume(&b, &Tolerance::default()).unwrap()
}

#[test]
fn a_half_pipe_thickens_outward_and_inward() {
    let (r, t, h) = (5.0_f64, 1.0_f64, 10.0_f64);
    let pi = std::f64::consts::PI;
    for (side, want) in [
        ("outward", pi * h * ((r + t).powi(2) - r * r) / 2.0),
        ("inward", pi * h * (r * r - (r - t).powi(2)) / 2.0),
    ] {
        let v = volume(side);
        assert!((v - want).abs() < 1e-6 * want, "{side}: {v} want {want}");
    }
}
