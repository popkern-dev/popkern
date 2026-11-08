//! Drafts refuse a body that is not a solid by its kind (PK-1101 S5).

use popkern_core::error::CoreError;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

fn tol() -> Tolerance {
    Tolerance::default()
}

use popkern_core::geometry::{Point3, Vec3};
use popkern_draft::{
    AngleVariation, DraftError, DraftOptions, DraftSpec, NeutralReference, draft_faces,
};

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
fn a_draft_refuses_a_general_body() {
    let body = general_body();
    let face = body.faces.keys().next().expect("a face");
    let spec = DraftSpec {
        face,
        angle: 0.05,
        pull_direction: Vec3::new(0.0, 0.0, 1.0),
        neutral: NeutralReference::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        },
        angle_variation: AngleVariation::Constant,
    };
    let result = draft_faces(&body, &[spec], &DraftOptions::default(), tol(), None);
    assert!(
        matches!(
            result,
            Err(DraftError::Core(CoreError::UnsupportedBodyKind {
                operation: "draft",
                found: BodyKind::General,
                ..
            }))
        ),
        "{result:?}"
    );
}
