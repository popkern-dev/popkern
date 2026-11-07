//! Every boolean operation refuses a general body by its kind
//! (PK-1101 S5), before any geometry is read.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

use crate::BoolOperation;
use crate::error::BoolError;

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

fn cube() -> BRep {
    BRep::unit_cube(tol()).expect("cube")
}

/// `result` is the kind gate's refusal of operand `operand` of
/// `operation`, the body found general.
fn assert_refused<T: std::fmt::Debug>(
    result: Result<T, BoolError>,
    operation: &str,
    operand: usize,
) {
    match result {
        Err(BoolError::UnsupportedBodyKind {
            operation: op,
            operand: i,
            kind: BodyKind::General,
            ..
        }) => {
            assert_eq!((op, i), (operation, operand));
        }
        other => panic!("{operation}: expected the kind gate, got {other:?}"),
    }
}

#[test]
fn booleans_refuse_a_general_body() {
    for op in [
        BoolOperation::Union,
        BoolOperation::Intersection,
        BoolOperation::Subtraction,
    ] {
        assert_refused(
            crate::boolean(op, cube(), general_body(), tol(), None),
            "boolean",
            1,
        );
        assert_refused(
            crate::boolean(op, general_body(), cube(), tol(), None),
            "boolean",
            0,
        );
        assert_refused(
            crate::boolean_with_sheets(op, &general_body(), &cube(), tol()),
            "sheet boolean",
            0,
        );
    }
    let message = crate::union(cube(), general_body(), tol(), None)
        .expect_err("refused")
        .to_string();
    assert_eq!(message, "boolean needs solids; operand 1 is a general body");
}

#[test]
fn partition_operations_refuse_a_general_body() {
    assert_refused(
        crate::partition(&[&cube(), &general_body()], tol()).map(|_| ()),
        "partition",
        1,
    );
    assert_refused(
        crate::ops::imprint::imprint(&general_body(), &cube(), &tol()).map(|_| ()),
        "imprint",
        0,
    );
}

#[test]
fn plane_operations_refuse_a_general_body() {
    let body = general_body();
    let origin = Point3::new(0.5, 0.5, 0.5);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    assert_refused(
        crate::ops::split_slice::split_by_plane(
            &body,
            &crate::ops::split_slice::SplitPlane { origin, normal },
            &tol(),
        )
        .map(|_| ()),
        "split by plane",
        0,
    );
    assert_refused(
        crate::sheet::split::split_solid_by_plane(&body, origin, normal, tol()),
        "split by plane",
        0,
    );
    assert_refused(
        crate::ops::plane_section::plane_section(
            &body,
            &crate::ops::plane_section::SectionPlane {
                origin,
                normal,
                u_axis: Vec3::new(1.0, 0.0, 0.0),
            },
            0.01,
            &tol(),
        )
        .map(|_| ()),
        "plane section",
        0,
    );
    let face = cube().faces.keys().next().expect("a face");
    assert_refused(
        crate::ops::face_level::face_level_boolean(
            &cube(),
            &body,
            BoolOperation::Union,
            &[crate::ops::face_level::FacePairSpec {
                face_a: face,
                face_b: face,
            }],
            &tol(),
        )
        .map(|_| ()),
        "face-level boolean",
        1,
    );
}
