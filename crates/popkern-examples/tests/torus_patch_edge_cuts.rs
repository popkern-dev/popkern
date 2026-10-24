//! Booleans that cut the primitive torus exactly along its patch edges
//! (PK-942).
//!
//! `BRep::torus` is a 4 x 4 grid of patches whose edges are meridians
//! at u = k pi/2 and parallels at v = k pi/2. A box face through the
//! axis cuts along two meridian edges, a plane through the centre
//! perpendicular to the axis along both equators. Both bodies are
//! rotated together, so every orientation is the same scene, but a
//! rotation moves the chart seam off the patch edges (`BRep::transform`
//! shifts the pcurves). Before the fix, 11 of 16 turns about the axis
//! declined the quadrant cut and one served a valid body whose torus
//! area was 5.3x off; the equatorial half-space declined under every
//! turn and the mixed rotation.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_measure::face_area;

const PI: f64 = std::f64::consts::PI;
/// Whole torus area 4 pi^2 R r for R = 1, r = 0.3.
const WHOLE: f64 = 4.0 * PI * PI * 0.3;

fn rotations() -> Vec<(String, Transform3)> {
    vec![
        ("unrotated".into(), Transform3::identity()),
        ("z 0.7".into(), Transform3::rotation_z(0.7)),
        ("x 1.1".into(), Transform3::rotation_x(1.1)),
        ("y -0.4".into(), Transform3::rotation_y(-0.4)),
        (
            "mixed".into(),
            Transform3::rotation_x(0.3)
                .compose(&Transform3::rotation_y(1.9))
                .compose(&Transform3::rotation_z(-2.2)),
        ),
    ]
}

/// Turns about the torus axis, 0.1 to 1.5.
fn axis_turns() -> Vec<(String, Transform3)> {
    (1..=15)
        .map(|k| {
            let t = 0.1 * k as f64;
            (format!("z {t:.1}"), Transform3::rotation_z(t))
        })
        .collect()
}

/// A box of side 4 at `offset`.
fn block(offset: (f64, f64, f64), tol: Tolerance) -> BRep {
    BRep::unit_cube(tol)
        .and_then(|c| c.transform(&Transform3::scale_uniform(4.0)))
        .and_then(|c| c.transform(&Transform3::translation(offset.0, offset.1, offset.2)))
        .expect("block")
}

/// The x > 0, y > 0 quadrant: faces through the axis on two meridian
/// edges.
const QUADRANT: (f64, f64, f64) = (0.0, 0.0, -2.0);
/// The z > 0 half-space: a face on both equators.
const UPPER: (f64, f64, f64) = (-2.0, -2.0, 0.0);

fn assert_serves(
    what: &str,
    op: BoolOperation,
    offset: (f64, f64, f64),
    turns: &[(String, Transform3)],
    want: f64,
) {
    let tol = Tolerance::default();
    for (name, r) in turns {
        let torus = BRep::torus(1.0, 0.3, tol)
            .and_then(|t| t.transform(r))
            .expect("torus");
        let tool = block(offset, tol).transform(r).expect("rotate");
        let body =
            boolean(op, torus, tool, tol, None).unwrap_or_else(|e| panic!("{what} {name}: {e}"));
        body.validate()
            .unwrap_or_else(|e| panic!("{what} {name}: invalid: {e:?}"));
        let area: f64 = body
            .faces
            .iter()
            .filter(|(_, f)| {
                f.surface
                    .as_ref()
                    .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Torus { .. }))
            })
            .map(|(h, _)| face_area(&body, h, &tol).expect("area"))
            .sum();
        assert!(
            (area - want).abs() < 1e-10 * want,
            "{what} {name}: torus area {area} vs {want}"
        );
    }
}

#[test]
fn a_quadrant_cut_on_the_meridian_edges_serves_in_every_orientation() {
    let all: Vec<_> = rotations().into_iter().chain(axis_turns()).collect();
    assert_serves(
        "quadrant cut",
        BoolOperation::Subtraction,
        QUADRANT,
        &all,
        0.75 * WHOLE,
    );
    let some = rotations();
    assert_serves(
        "quadrant kept",
        BoolOperation::Intersection,
        QUADRANT,
        &some,
        0.25 * WHOLE,
    );
    assert_serves(
        "quadrant union",
        BoolOperation::Union,
        QUADRANT,
        &some,
        0.75 * WHOLE,
    );
}

#[test]
fn an_equatorial_half_space_serves_in_every_orientation() {
    let all: Vec<_> = rotations().into_iter().chain(axis_turns()).collect();
    assert_serves(
        "lower half",
        BoolOperation::Subtraction,
        UPPER,
        &all,
        0.5 * WHOLE,
    );
    assert_serves(
        "upper half",
        BoolOperation::Intersection,
        UPPER,
        &rotations(),
        0.5 * WHOLE,
    );
}
