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
/// tolerance) serves exactly. Measured on the way (PK-838, PK-860):
/// 46 s (corner-cap self-intersection checks marched isolated
/// contacts through the whole 200k-step budget); a body check 7
/// falsely rejected (the raw chart normal against f64::EPSILON); then
/// a self-intersecting body 2% short, because the corner carve matched
/// anchors within 100 tolerances, which below r/tol = 241 swallowed the
/// rail crossings at r(sqrt2 - 1) from them; then exact but reported
/// G0, because pcurve recompute read every cylinder under r = 1e-4 as
/// a chart pole and fitted the strips' corner arcs 3.4e-6 off.
#[test]
fn a_rounded_cube_seventy_tolerances_across_serves_exact_and_fast() {
    let (s, r) = (0.0016, 0.04375 * 0.0016);
    let watch = std::time::Instant::now();
    let (brep, downgraded) = rounded_cube("m", s, r);
    assert!(
        watch.elapsed().as_secs_f64() < 10.0,
        "took {:?}",
        watch.elapsed()
    );
    assert_eq!(downgraded, 0, "the G1 fillet reads as downgraded");
    let issues = brep.validate_all();
    assert!(issues.is_empty(), "{issues:?}");
    let volume = popkern_measure::body_volume(&brep, &brep.tolerance).expect("volume");
    let expected = exact_volume(s, r);
    assert!(
        ((volume - expected) / expected).abs() < 1e-6,
        "volume {volume} vs {expected}"
    );
}

/// A corner radius within the tolerance band declines fast with the
/// corner-resolution error.
#[test]
fn a_corner_radius_inside_the_tolerance_declines_fast() {
    let (s, r) = (0.0016, 5e-7);
    let source = format!(
        "let m = model();\n\
         let h = {s} / 2.0;\n\
         let die = m.profile(xy_at(-h)).rect({s}, {s}).sketch(\"block\").extrude(blind({s}));\n\
         die = die.fillet(edges(), {r}, #{{ label: \"rounds\" }});\nm"
    );
    let watch = std::time::Instant::now();
    let error = match ScriptHost::new().run(&source, "sub_tolerance.rhai") {
        Ok(_) => panic!("a sub-tolerance corner radius serves"),
        Err(e) => e.render(),
    };
    assert!(error.contains("linear tolerance"), "{error}");
    assert!(
        watch.elapsed().as_secs_f64() < 5.0,
        "took {:?}",
        watch.elapsed()
    );
}

/// The same 1.6 mm cube at a 1e-7 tolerance (700 tolerances across its
/// 70 um rounds, inside the corner resolution) serves exact and clean.
/// The corner-cap contact collapse read "on the ball" against an
/// absolute 1e-3, which every vertex of this cube satisfied: it
/// contracted the wrong chords and the rounds' end arcs carried
/// reversed pcurves 3.8e-5 off (PK-860).
#[test]
fn a_tiny_rounded_cube_at_a_fine_tolerance_serves_exact() {
    let (s, r) = (0.0016, 0.04375 * 0.0016);
    let source = format!(
        "let m = model(#{{ tolerance: 1e-7 }});\n\
         let h = {s} / 2.0;\n\
         let die = m.profile(xy_at(-h)).rect({s}, {s}).sketch(\"block\").extrude(blind({s}));\n\
         die = die.fillet(edges(), {r}, #{{ label: \"rounds\" }});\nm"
    );
    let mut run = ScriptHost::new()
        .run(&source, "tiny_fine.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    let brep = run.model.take().expect("model").to_brep().expect("brep");
    let issues = brep.validate_all();
    assert!(issues.is_empty(), "{issues:?}");
    let volume = popkern_measure::body_volume(&brep, &brep.tolerance).expect("volume");
    let expected = exact_volume(s, r);
    assert!(
        ((volume - expected) / expected).abs() < 1e-6,
        "volume {volume} vs {expected}"
    );
}

/// The same 70 um rounds at coarse tolerances, down to two tolerances
/// across, serve exact and clean. At ten tolerances validate_all
/// reported the two strips meeting the corner ball as crossing: the
/// cylinder pair check counted grid points within ten tolerances of
/// the other strip (a full radius here) as contacts (PK-866).
#[test]
fn a_tiny_rounded_cube_serves_at_coarse_tolerances() {
    let (s, r) = (0.0016, 0.04375 * 0.0016);
    for tolerance in [7e-6, 2e-5, 3.5e-5] {
        let source = format!(
            "let m = model(#{{ tolerance: {tolerance} }});\n\
             let h = {s} / 2.0;\n\
             let die = m.profile(xy_at(-h)).rect({s}, {s}).sketch(\"block\").extrude(blind({s}));\n\
             die = die.fillet(edges(), {r}, #{{ label: \"rounds\" }});\nm"
        );
        let mut run = ScriptHost::new()
            .run(&source, "tiny_coarse.rhai")
            .unwrap_or_else(|e| panic!("tol {tolerance}: {}", e.render()));
        let brep = run.model.take().expect("model").to_brep().expect("brep");
        let issues = brep.validate_all();
        assert!(issues.is_empty(), "tol {tolerance}: {issues:?}");
        let volume = popkern_measure::body_volume(&brep, &brep.tolerance).expect("volume");
        let expected = exact_volume(s, r);
        assert!(
            ((volume - expected) / expected).abs() < 1e-6,
            "tol {tolerance}: volume {volume} vs {expected}"
        );
    }
}
