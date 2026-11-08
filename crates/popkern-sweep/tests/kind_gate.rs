//! Sweeps into a body refuse a target that is not a solid by its kind
//! (PK-1101 S5), before the tool is swept.

use popkern_core::error::CoreError;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

fn tol() -> Tolerance {
    Tolerance::default()
}

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::Point3;
use popkern_sweep::{SweepError, SweepOptions, SweepProfile, sweep_add, sweep_cut};

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

fn line(a: Point3, b: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    }
}

#[test]
fn sweeps_into_a_general_body_are_refused() {
    let body = general_body();
    let profile = SweepProfile::Wire(line(Point3::new(0.0, 0.0, 0.0), Point3::new(0.1, 0.0, 0.0)));
    let path = line(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 1.0));
    for (operation, result) in [
        (
            "sweep add",
            sweep_add(
                &body,
                &profile,
                &path,
                &SweepOptions::default(),
                tol(),
                None,
            ),
        ),
        (
            "sweep cut",
            sweep_cut(
                &body,
                &profile,
                &path,
                &SweepOptions::default(),
                tol(),
                None,
            ),
        ),
    ] {
        match result {
            Err(SweepError::Core(CoreError::UnsupportedBodyKind {
                operation: op,
                found: BodyKind::General,
                ..
            })) => assert_eq!(op, operation),
            other => panic!("{operation}: expected the kind gate, got {other:?}"),
        }
    }
}
