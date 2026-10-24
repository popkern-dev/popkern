//! A sphere drilled across and along its pole axis, served by
//! graph-merge (PK-940, PK-944).
//!
//! The cellular path serves this fixture (PK-937); graph-merge is its
//! fallback and must serve it too. A hole circle is one closed record
//! whose start node joins its last piece to its first. The dead-end
//! weld's far tier read that pass-through as a record ending in both
//! directions, and where the start sat near a sphere-seam junction on
//! the same circle (5.4e-2 away on the mixed rotation) it hauled the
//! start onto it: both hole halves then ended where they began, the
//! face's outer loop no longer chained, and the face served wound
//! inside out. Before the fix 26 of these 315 cases declined and one
//! served a valid body whose sphere area was 5.4e-7 off.
//!
//!
//! A bore along the pole axis cuts two latitude rims, which the
//! seam meridian and the chart's own seam split into pieces. The exact
//! rim pass re-spells each piece as an arc anchored at the curve's
//! ends, and those sat 4e-10 (2e-9 in angle) from the vertices they
//! meet: on a piece whose forward half-edge runs its curve backward,
//! the vertex anchoring missed and the rim's pcurve stopped short of
//! the seam meridian's, so the face's parameter-space boundary did not
//! close and the sphere area came out 3.4e-10 off under rotation_z(0.7)
//! and 1e-10 under rotation_y(-0.4) (PK-944).
//!
//! This binary holds one test: it forces graph-merge for its whole
//! process (nextest runs each test in its own process).

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_measure::face_area;

const PI: f64 = std::f64::consts::PI;

fn rotations() -> Vec<(&'static str, Transform3)> {
    vec![
        ("unrotated", Transform3::identity()),
        ("z 0.7", Transform3::rotation_z(0.7)),
        ("x 1.1", Transform3::rotation_x(1.1)),
        ("y -0.4", Transform3::rotation_y(-0.4)),
        (
            "mixed",
            Transform3::rotation_x(0.3)
                .compose(&Transform3::rotation_y(1.9))
                .compose(&Transform3::rotation_z(-2.2)),
        ),
    ]
}

#[test]
fn graph_merge_serves_a_sphere_drilled_across_and_along_its_pole_axis() {
    // SAFETY-adjacent note: set_var is fine here — this binary's only
    // test, and nextest runs it in its own process.
    unsafe { std::env::set_var("POPKERN_CELLULAR", "0") };
    let tol = Tolerance::default();
    // Two caps of height 0.5 - sqrt(0.21) leave the radius-0.5 sphere.
    let want = PI - 2.0 * 2.0 * PI * 0.5 * (0.5 - 0.21_f64.sqrt());
    // The bore turns about the sphere's pole axis: across it, at turn 0
    // on the sphere's seam meridian (63 turns); along it, only the
    // bore's own seam moves (16 turns).
    let across = (0..63).map(|k| {
        let turn = 0.1 * f64::from(k);
        let tilt = Transform3::rotation_z(turn).compose(&Transform3::rotation_y(PI / 2.0));
        (format!("across, turn {turn:.1}"), tilt)
    });
    let along = (0..16).map(|k| {
        let turn = 0.1 * f64::from(k);
        (
            format!("along, turn {turn:.1}"),
            Transform3::rotation_z(turn),
        )
    });
    for (turn, tilt) in across.chain(along) {
        for (name, r) in rotations() {
            let ball = BRep::sphere(0.5, tol)
                .and_then(|s| s.transform(&r))
                .expect("sphere");
            let bore = BRep::cylinder(0.2, 2.0, tol)
                .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -1.0)))
                .and_then(|c| c.transform(&tilt))
                .and_then(|c| c.transform(&r))
                .expect("bore");
            let drilled = boolean(BoolOperation::Subtraction, ball, bore, tol, None)
                .unwrap_or_else(|e| panic!("{turn}, {name}: {e}"));
            drilled
                .validate()
                .unwrap_or_else(|e| panic!("{turn}, {name}: invalid: {e:?}"));
            let got: f64 = drilled
                .faces
                .iter()
                .filter(|(_, f)| {
                    f.surface
                        .as_ref()
                        .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Sphere { .. }))
                })
                .map(|(h, _)| face_area(&drilled, h, &tol).expect("area"))
                .sum();
            assert!(
                (got - want).abs() < 1e-10 * want,
                "{turn}, {name}: drilled sphere {got} vs {want}"
            );
        }
    }
    let fallbacks = popkern_bool::merge_fallback::take_merge_fallback_events();
    assert!(
        !fallbacks.is_empty(),
        "the kill switch must route every case through graph-merge"
    );
}
