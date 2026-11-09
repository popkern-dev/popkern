//! Records of different body pairs touching on a face they share, in
//! the N-body partition (PK-1130's acceptance fixture for moving the
//! partition's own crossing split onto the face arrangement): a box
//! [-10, 10]² × [-5, 0] and two balls of radius 3 on its top face at
//! (±3, 0, 0). Each ball cuts the top face in a circle; the two circles,
//! of the pairs (box, ball A) and (box, ball B), touch at the origin,
//! where the balls touch too. The union is the box and two domes, its
//! top face touching itself at the origin (a pinch use): one valid body
//! of volume 2000 + 36π, at identity and turned, under both policies.

use std::f64::consts::PI;

use popkern_bool::{BoolOptions, NonManifoldPolicy, SetExpr, nary_boolean};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn two_circles_of_different_pairs_touching_on_a_face_serve() {
    let tol = Tolerance::default();
    let rect = |z: f64| {
        vec![
            Point3::new(-10.0, -10.0, z),
            Point3::new(10.0, -10.0, z),
            Point3::new(10.0, 10.0, z),
            Point3::new(-10.0, 10.0, z),
        ]
    };
    let turned = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.7, 0.3),
        0.645_771_823_237_9,
    )
    .expect("axis");
    let want = 2000.0 + 36.0 * PI;
    for place in [Transform3::identity(), turned] {
        let bx = BRep::polygon_frustum(&rect(-5.0), &rect(0.0), tol)
            .and_then(|b| b.transform(&place))
            .expect("box");
        let ball = |x: f64| {
            BRep::sphere(3.0, tol)
                .and_then(|b| b.transform(&Transform3::translation(x, 0.0, 0.0)))
                .and_then(|b| b.transform(&place))
                .expect("ball")
        };
        let (a, b) = (ball(-3.0), ball(3.0));
        for policy in [NonManifoldPolicy::Decline, NonManifoldPolicy::General] {
            let body = nary_boolean(
                &[&bx, &a, &b],
                &SetExpr::Any,
                tol,
                &BoolOptions {
                    non_manifold: policy,
                },
            )
            .expect("serves");
            assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
            let v = popkern_measure::body_volume(&body, &tol).expect("volume");
            assert!(
                (v - want).abs() < 1e-6 * want,
                "{policy:?}: volume {v} for {want}"
            );
        }
    }
}
