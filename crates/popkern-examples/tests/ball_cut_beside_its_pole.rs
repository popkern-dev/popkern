//! A ball cut by a plane whose circle passes beside its pole (PK-1184).
//! `BRep::sphere` (radius 3) has its poles on z; a half-space whose
//! plane is tilted 0.3 or 0.6 rad from the xy plane removes a cap whose
//! circle passes `eps` from the south pole (the pole removed for
//! `eps > 0`, kept for `eps < 0`). With the ball turned about z the
//! seam no longer meets the circle's nearest point to the pole, and
//! the circle encircling the pole sweeps a whole turn of longitude
//! within `eps` of it: the chart walk's projection continued it
//! station to station, lost the turn (17 samples spanning 2.99 rad of
//! the 2 pi), found no window for the seam's borders, and kept the cap
//! inside out (a valid-looking boolean, two faces, the volume 98% to
//! 99.9% short).

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

const R: f64 = 3.0;

/// The ball, turned `turn` about z, less the half-space `n . x > p`, `n`
/// tilted `alpha` from -z, `p` placing the circle `eps` from the pole;
/// and `p`.
fn cut(alpha: f64, eps: f64, turn: f64) -> (BRep, f64) {
    let tol = Tolerance::default();
    let n = Vec3::new(alpha.sin(), 0.0, -alpha.cos());
    let p = R * alpha.cos() - eps;
    let ball = BRep::sphere(R, tol)
        .and_then(|b| b.transform(&Transform3::rotation_z(turn)))
        .expect("ball");
    let square = [
        Point3::new(-10.0, -10.0, 0.0),
        Point3::new(10.0, -10.0, 0.0),
        Point3::new(10.0, 10.0, 0.0),
        Point3::new(-10.0, 10.0, 0.0),
    ];
    let z = Vec3::new(0.0, 0.0, 1.0);
    let tilt =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), z.cross(&n), z.dot(&n).acos())
            .expect("axis");
    let slab = BRep::prism(&square, z, 20.0, tol)
        .and_then(|b| b.transform(&tilt))
        .and_then(|b| b.transform(&Transform3::translation(n.x * p, n.y * p, n.z * p)))
        .expect("slab");
    let body = boolean(BoolOperation::Subtraction, ball, slab, tol, None)
        .unwrap_or_else(|e| panic!("alpha {alpha} eps {eps} turn {turn}: {e}"));
    (body, p)
}

/// Every placement serves the ball less its cap: valid, its sphere face
/// the sphere's area less the cap's (`4 pi R^2 - 2 pi R h`) and its flat
/// face the cut disc (`pi (R^2 - p^2)`), each to 1e-6. (The face areas
/// judge the body; its volume as `body_volume` measures it near the pole
/// is PK-1178's.)
#[test]
fn a_ball_cut_beside_its_pole_keeps_the_ball() {
    let tol = Tolerance::default();
    let pi = std::f64::consts::PI;
    for alpha in [0.3, 0.6] {
        for eps in [1e-3, 1e-4, 1e-5, -1e-5, -1e-4, -1e-3] {
            for turn in [0.0, 0.4, 2.0] {
                let what = format!("alpha {alpha} eps {eps:e} turn {turn}");
                let (body, p) = cut(alpha, eps, turn);
                let issues = body.validate_all();
                assert!(issues.is_empty(), "{what}: {issues:?}");
                let h = R - p;
                assert_eq!(body.faces.len(), 2, "{what}");
                for (face, f) in body.faces.iter() {
                    let plane = matches!(
                        f.surface.as_ref().expect("surface").geometry,
                        SurfaceGeometry::Plane { .. }
                    );
                    let want = if plane {
                        pi * (R * R - p * p)
                    } else {
                        4.0 * pi * R * R - 2.0 * pi * R * h
                    };
                    let area = popkern_measure::face_area(&body, face, &tol).expect("area");
                    assert!(
                        ((area - want) / want).abs() < 1e-6,
                        "{what}: {} face area {area} for {want}",
                        if plane { "flat" } else { "sphere" }
                    );
                }
            }
        }
    }
}
