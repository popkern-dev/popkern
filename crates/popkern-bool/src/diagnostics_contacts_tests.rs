use popkern_core::geometry::intersection::contact_classify::SurfaceContactKind;
use popkern_core::geometry::{Point3, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::{BoolOperation, boolean_with_diagnostics};

fn cube(size: f64, at: (f64, f64, f64)) -> BRep {
    BRep::unit_cube(Tolerance::default())
        .and_then(|c| c.transform(&Transform3::scale_uniform(size)))
        .and_then(|c| c.transform(&Transform3::translation(at.0, at.1, at.2)))
        .expect("cube")
}

#[test]
fn an_inscribed_ball_touches_each_face_at_its_middle() {
    let tol = Tolerance::default();
    let ball = BRep::sphere(5.0, tol)
        .and_then(|b| b.transform(&Transform3::translation(5.0, 5.0, 5.0)))
        .expect("ball");
    let (_, diag) = boolean_with_diagnostics(
        BoolOperation::Intersection,
        &cube(10.0, (0.0, 0.0, 0.0)),
        &ball,
        tol,
    )
    .expect("serves the ball");
    let mut at: Vec<Point3> = diag.contacts.iter().map(|c| c.location).collect();
    at.dedup_by(|p, q| p.distance_to(q) <= 1e-6);
    assert_eq!(at.len(), 6, "{:?}", diag.contacts);
    for c in &diag.contacts {
        assert!(
            matches!(
                c.classification.kind,
                SurfaceContactKind::TangentPoint { .. }
            ),
            "{c:?}"
        );
        let off = c.location - Point3::new(5.0, 5.0, 5.0);
        let axes = [off.x.abs(), off.y.abs(), off.z.abs()];
        assert!(
            axes.iter().filter(|d| (**d - 5.0).abs() <= 1e-6).count() == 1
                && axes.iter().filter(|d| **d <= 1e-6).count() == 2,
            "a face middle: {c:?}"
        );
    }
}

#[test]
fn crossing_cubes_touch_nowhere() {
    let (_, diag) = boolean_with_diagnostics(
        BoolOperation::Union,
        &cube(1.0, (0.0, 0.0, 0.0)),
        &cube(1.0, (0.5, 0.5, 0.5)),
        Tolerance::default(),
    )
    .expect("unites");
    assert!(diag.contacts.is_empty(), "{:?}", diag.contacts);
}
