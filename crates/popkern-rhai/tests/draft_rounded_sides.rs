//! A box's sides with rounded vertical edges drafted together: each
//! quarter cylinder drafts into a quarter cone narrowing along the pull
//! whichever way its own axis runs, and the vertical edge where it meets
//! a drafted plane recomputes as the cone's ruling in that plane. Left
//! stale, those edges opened the body by 10·tan 5° = 0.8749 at its top.

use popkern_core::tolerance::Tolerance;

#[test]
fn a_rounded_box_drafted_reads_its_closed_form() {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(
            r#"let m = model();
let b = m.profile(xy()).rect(30.0, 20.0).at(15.0, 10.0).sketch("box").extrude(blind(10.0));
b = b.fillet(edges().direction("z"), 4.0, #{ label: "rounds" });
b = b.draft(faces().except(faces().facing("z").or(faces().facing_neg("z"))), 5.0, xy(), #{ label: "sides" });
m"#,
            "draft_rounded_sides",
        )
        .unwrap_or_else(|e| panic!("{e}"));
    let b = run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone();
    let errors = b.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
    // A₀H − P₀tH²/2 + πt²H³/3: A₀ = ab − (4 − π)r², P₀ = 2(a + b) − (8 − 2π)r.
    let (a, bb, r, h) = (30.0_f64, 20.0_f64, 4.0_f64, 10.0_f64);
    let t = (5.0_f64).to_radians().tan();
    let pi = std::f64::consts::PI;
    let a0 = a * bb - (4.0 - pi) * r * r;
    let p0 = 2.0 * (a + bb) - (8.0 - 2.0 * pi) * r;
    let want = a0 * h - p0 * t * h * h / 2.0 + pi * t * t * h.powi(3) / 3.0;
    let v = popkern_measure::body_volume(&b, &Tolerance::default()).unwrap();
    assert!((v - want).abs() < 1e-6 * want, "{v} want {want}");
}

/// A box drilled through its ends, the ends drafted 8 degrees about the
/// floor: each hole rim recomputes as the exact ellipse the drafted
/// plane cuts from the hole's wall, its vertex slides along the hole's
/// seam to the plane, and the seam re-reads its range.
#[test]
fn a_cross_drilled_box_with_drafted_ends_reads_its_closed_form() {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(
            r#"let m = model();
let b = m.profile(xy()).rect(30.0, 20.0).at(15.0, 10.0).sketch("box").extrude(blind(20.0));
let drill = m.tool(|m| m.profile(yz_at(-1.0)).circle(4.0).at(10.0, 10.0).sketch("drill").extrude(blind(32.0)));
b = b.subtract(drill, #{ label: "hole" });
b = b.draft(faces().facing("x").or(faces().facing_neg("x")), 8.0, xy(), #{ label: "ends" });
m"#,
            "draft_cross_drilled",
        )
        .unwrap_or_else(|e| panic!("{e}"));
    let b = run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone();
    let errors = b.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
    // b(aH − tH²) − πr²(a − 2tz₀)
    let (a, bb, h, z0, r) = (30.0_f64, 20.0_f64, 20.0_f64, 10.0_f64, 4.0_f64);
    let t = (8.0_f64).to_radians().tan();
    let want = bb * (a * h - t * h * h) - std::f64::consts::PI * r * r * (a - 2.0 * t * z0);
    let v = popkern_measure::body_volume(&b, &Tolerance::default()).unwrap();
    assert!((v - want).abs() < 1e-6 * want, "{v} want {want}");
}
