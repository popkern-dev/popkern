//! A rounded cube is the same shape at every scale: filleting all
//! twelve edges of a cube with r = 0.04375 s serves, validates, and
//! encloses the exact rounded-cube volume from s = 8 mm to 1.6 m, and in
//! a millimetre model.
//!
//! The corner surface repair skipped any face whose current surface fit
//! its vertices within an absolute 0.01, while a corner's stale
//! cube-corner vertex sits r(√2−1) = 0.414 r off the true sphere: below
//! r ≈ 0.024 the outlier hid inside the window, the corner kept its
//! wrong surface and the fillet declined (cubes 0.008 to 0.25 failed,
//! 0.3 and up served; a 0.7 mm round on a 16 mm die failed outright).
//! It also achieves the G1 it asks for at every scale.

use popkern_core::topology::BRep;
use popkern_rhai::ScriptHost;

/// The rounded cube and how many blends delivered less continuity than
/// they asked for.
fn rounded_cube(unit: &str, s: f64, r: f64) -> (BRep, usize) {
    let source = format!(
        "let m = model(#{{ unit: \"{unit}\" }});\n\
         let h = {s} / 2.0;\n\
         let die = m.profile(xy_at(-h)).rect({s}, {s}).sketch(\"block\").extrude(blind({s}));\n\
         die = die.fillet(edges(), {r}, #{{ label: \"rounds\" }});\nm"
    );
    let mut host = ScriptHost::new();
    let mut outcome = host
        .run(&source, "rounded_cube.rhai")
        .unwrap_or_else(|e| panic!("s {s} {unit}: {}", e.render()));
    let model = outcome.model.as_mut().expect("model");
    let downgraded = model.report().downgraded_blends().count();
    let brep = model
        .to_part()
        .unwrap_or_else(|e| panic!("s {s} {unit}: {e}"))
        .brep()
        .clone();
    (brep, downgraded)
}

/// The exact volume of a cube of side `s` with every edge rounded `r`:
/// the inner box, six face slabs, twelve quarter-cylinders and eight
/// sphere octants.
fn exact_volume(s: f64, r: f64) -> f64 {
    let pi = std::f64::consts::PI;
    let core = s - 2.0 * r;
    core.powi(3) + 6.0 * r * core * core + 3.0 * pi * r * r * core + 4.0 / 3.0 * pi * r.powi(3)
}

fn check(unit: &str, metres_per_unit: f64, s: f64) {
    let r = 0.04375 * s;
    let (brep, downgraded) = rounded_cube(unit, s, r);
    assert!(brep.validate().is_ok(), "s {s} {unit}: invalid");
    // The seam normals are read scale-free: a 0.7 mm corner sphere's
    // finite-difference normal fell under f64::EPSILON and the G1
    // fillet reported itself downgraded to G0.
    assert_eq!(
        downgraded, 0,
        "s {s} {unit}: the G1 fillet reads as downgraded"
    );
    let volume = popkern_measure::body_volume(&brep, &brep.tolerance).expect("volume");
    let expected = exact_volume(s * metres_per_unit, r * metres_per_unit);
    assert!(
        ((volume - expected) / expected).abs() < 1e-4,
        "s {s} {unit}: volume {volume} vs {expected}"
    );
}

#[test]
fn a_rounded_cube_serves_at_every_scale() {
    for s in [0.008, 0.016, 0.16, 0.25, 0.3, 1.6] {
        check("m", 1.0, s);
    }
}

#[test]
fn a_millimetre_die_takes_its_rounds() {
    check("mm", 1e-3, 16.0);
}

/// A 1.6 mm cube in a metre model (r = 70 um, 70x the 1e-6 linear
/// tolerance) declines fast with the corner-resolution error instead of
/// serving a corner the kernel cannot resolve. Measured before the
/// gate: 46 s (corner-cap self-intersection checks marched isolated
/// contacts through the whole 200k-step budget), then a body check 7
/// falsely rejected (it read the raw chart normal against
/// f64::EPSILON); with check 7 scale-free, a self-intersecting body 2%
/// short of the exact volume. Lifting the gate is PK-860.
#[test]
fn a_corner_radius_near_the_tolerance_declines_fast() {
    let (s, r) = (0.0016, 0.04375 * 0.0016);
    let source = format!(
        "let m = model();\n\
         let h = {s} / 2.0;\n\
         let die = m.profile(xy_at(-h)).rect({s}, {s}).sketch(\"block\").extrude(blind({s}));\n\
         die = die.fillet(edges(), {r}, #{{ label: \"rounds\" }});\nm"
    );
    let watch = std::time::Instant::now();
    let error = match ScriptHost::new().run(&source, "tiny_cube.rhai") {
        Ok(_) => panic!("a 70-tolerance corner radius serves"),
        Err(e) => e.render(),
    };
    assert!(error.contains("linear tolerance"), "{error}");
    assert!(
        watch.elapsed().as_secs_f64() < 5.0,
        "took {:?}",
        watch.elapsed()
    );
}
