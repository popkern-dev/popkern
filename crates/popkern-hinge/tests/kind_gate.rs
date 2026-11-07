//! A seam hinge refuses a target that is not a solid by its kind
//! (PK-1101 S5), before any knuckle is built.

use popkern_bool::BoolError;
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};
use popkern_hinge::{
    HingeComponent, HingeError, HingeKind, HingeParams, Seam, seam_hinge_component,
};

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

#[test]
fn a_seam_hinge_refuses_a_general_body() {
    let seam = Seam {
        origin: Point3::new(0.0, 1.0, 1.0),
        axis: Vec3::new(1.0, 0.0, 0.0),
        outward: Vec3::new(0.0, 1.0, 0.0),
        up: Vec3::new(0.0, 0.0, 1.0),
    };
    let params = HingeParams {
        pin_diameter: 0.1,
        thickness: 0.05,
        ..HingeParams::print_in_place(0.8)
    };
    let result = seam_hinge_component(
        HingeKind::PrintInPlace,
        HingeComponent::LeafA,
        &seam,
        &params,
        Some(&general_body()),
        tol(),
    );
    match result {
        Err(HingeError::Boolean { stage, source }) => {
            assert_eq!(stage, "seam hinge");
            assert!(matches!(
                *source,
                BoolError::UnsupportedBodyKind {
                    kind: BodyKind::General,
                    ..
                }
            ));
        }
        other => panic!("expected the kind gate, got {other:?}"),
    }
}
