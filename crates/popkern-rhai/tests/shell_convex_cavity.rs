//! PK-1212: a convex body shelled where its cavity never reaches the
//! opening, is empty, or runs past a consumed round is the body less its
//! cavity, the intersection of its kept faces' offset regions; and a
//! sphere cap's rim solves where its plane meets the offset sphere at a
//! shallow angle (Newton after projection). Each reads its closed form.

use popkern_core::tolerance::Tolerance;
use std::f64::consts::PI;

fn volume(name: &str, body: &str) -> f64 {
    let source = format!("let m = model();\n{body}\nm");
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(&source, name)
        .unwrap_or_else(|e| panic!("{name}: {e}"));
    let parts = run.model.expect("model").to_parts_all().expect("parts");
    parts
        .iter()
        .map(|p| {
            let b = p.brep();
            let errors = b.validate_all();
            assert!(errors.is_empty(), "{name}: {errors:?}");
            popkern_measure::body_volume(b, &Tolerance::default()).unwrap()
        })
        .sum()
}

fn cap(h: f64) -> String {
    format!(
        "let h = {h};\n\
         let b = m.profile(xz()).path(0.0, -5.0).arc_to([5.0, 0.0], [sqrt(25.0 - h * h), h]).line_to(0.0, h)\n\
             .close().sketch(\"cap\").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());\n\
         b = b.shell(faces().kind(\"plane\").at_height(h), 1.0);"
    )
}

fn rounded_plate(t: f64) -> String {
    format!(
        "let p = m.profile(xy()).rect(20.0, 14.0).at(10.0, 7.0).sketch(\"plate\").extrude(blind(10.0));\n\
         p = p.fillet(edges().straight().direction(\"z\"), 2.0);\n\
         p = p.shell(faces().facing(\"z\").at_height(10.0), {t});"
    )
}

fn check(name: &str, body: &str, want: f64) {
    let v = volume(name, body);
    assert!((v - want).abs() < 1e-6 * want, "{name}: {v} want {want}");
}

#[test]
fn a_ball_cut_below_its_offset_sphere_keeps_a_spherical_void() {
    check("sealed", &cap(4.5), 641.0 * PI / 8.0);
}

#[test]
fn a_ball_cut_just_under_its_offset_sphere_opens_a_pinhole() {
    check("pinhole", &cap(3.9999), 240.852_609_341_829_2);
}

#[test]
fn frustums_shelled_past_their_offset_apex() {
    let shell = |r0: f64, r1: f64, h: f64, t: f64, open: &str| {
        format!("let b = m.cone({r0}, {r1}, {h});\nb = b.shell({open}, {t});")
    };
    let top = |h: f64| format!("faces().facing(\"z\").at_height({h})");
    check(
        "sealed frustum",
        &shell(6.0, 3.0, 4.0, 2.7, &top(4.0)),
        83.904 * PI,
    );
    check(
        "tall frustum",
        &shell(7.0, 1.0, 8.0, 2.0, &top(8.0)),
        140.0 * PI,
    );
    check(
        "too thick",
        &shell(6.0, 3.0, 4.0, 3.5, &top(4.0)),
        84.0 * PI,
    );
    check(
        "open both ends",
        &shell(6.0, 3.0, 4.0, 2.7, "faces().kind(\"plane\")"),
        238.638_323_209_793_4,
    );
}

#[test]
fn a_rounded_plate_shelled_past_its_corner_radius_meets_in_sharp_corners() {
    check("consumed", &rounded_plate(3.0), 1856.0 + 40.0 * PI);
    check("to an edge", &rounded_plate(2.0), 1360.0 + 40.0 * PI);
}

#[test]
fn walls_that_meet_or_a_side_opened_whole() {
    check(
        "walls meeting",
        "let b = m.profile(xy()).rect(10.0, 8.0).at(5.0, 4.0).sketch(\"box\").extrude(blind(6.0));\n\
         b = b.shell(faces().facing(\"z\").at_height(6.0), 4.0);",
        480.0,
    );
    check(
        "side open",
        "let c = m.profile(xy()).circle(5.0).sketch(\"cylinder\").extrude(blind(10.0));\n\
         c = c.shell(faces().kind(\"cylinder\"), 1.0);",
        50.0 * PI,
    );
}
