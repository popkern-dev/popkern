//! Holes and threads refuse a target that is not a solid by its kind
//! (PK-1101 S5), before any tool body is swept.

use popkern_bool::BoolError;
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};
use popkern_fastener::holes::{clearance_hole, external_thread, tapped_hole};
use popkern_fastener::{ClearanceHoleOptions, FastenerError, HoleDepth, TapOptions, ThreadSpec};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// A cube beside a sheet: lumps of two kinds, a general body.
fn general_body() -> BRep {
    let mut body = BRep::unit_cube(tol()).expect("cube");
    let mut sheet = BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(3.0, 0.0, 0.0)))
        .expect("cube");
    let face = sheet.faces.keys().next().expect("a face");
    sheet.retire_faces_to_boundary(&[face]).expect("retires");
    sheet.settle_open_shells();
    body.append_body(&sheet);
    assert_eq!(body.body_kind(), BodyKind::General);
    body
}

fn assert_refused(result: Result<BRep, FastenerError>, operation: &str) {
    match result {
        Err(FastenerError::Boolean { stage, source }) => {
            assert_eq!(stage, operation);
            assert!(matches!(
                source,
                BoolError::UnsupportedBodyKind {
                    kind: BodyKind::General,
                    ..
                }
            ));
        }
        other => panic!("{operation}: expected the kind gate, got {other:?}"),
    }
}

#[test]
fn holes_and_threads_refuse_a_general_body() {
    let body = general_body();
    let spec = ThreadSpec::parse("M6").expect("M6");
    let origin = Point3::new(0.5, 0.5, 1.0);
    let down = Vec3::new(0.0, 0.0, -1.0);
    assert_refused(
        tapped_hole(
            &body,
            &spec,
            origin,
            down,
            HoleDepth::Blind(0.5),
            &TapOptions::default(),
            tol(),
        ),
        "tapped hole",
    );
    assert_refused(
        clearance_hole(
            &body,
            &spec,
            origin,
            down,
            HoleDepth::Blind(0.5),
            &ClearanceHoleOptions::default(),
            tol(),
        ),
        "clearance hole",
    );
    assert_refused(
        external_thread(&body, &spec, origin, down, 0.0, 0.5, tol()),
        "external thread",
    );
}
