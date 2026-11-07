//! Centroid and inertia refuse a general body by its kind (PK-1101 S5);
//! volume measures it as the sum of its solid regions (S8), area any.

use popkern_core::error::CoreError;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

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
fn body_measures_refuse_a_general_body() {
    let body = general_body();
    // The sheet bounds no volume: the cube's region is the body's.
    let volume = popkern_measure::body_volume(&body, &tol()).expect("volume");
    assert!((volume - 1.0).abs() < 1e-9, "{volume}");
    for result in [
        popkern_measure::body_centroid(&body, &tol()).map(|_| ()),
        popkern_measure::body_moment_of_inertia(&body, popkern_measure::Axis::Z, &tol())
            .map(|_| ()),
    ] {
        assert!(
            matches!(
                result,
                Err(popkern_measure::MeasureError::Core(
                    CoreError::UnsupportedBodyKind {
                        found: BodyKind::General,
                        ..
                    }
                ))
            ),
            "{result:?}"
        );
    }
    assert!(popkern_measure::body_surface_area(&body, &tol()).is_ok());
}
