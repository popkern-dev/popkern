use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::*;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    BRep::unit_cube(tol())
        .expect("cube builds")
        .transform(&Transform3::translation(x, y, z))
        .expect("translates")
}

#[test]
fn contact_report_names_point_curve_and_area_touches() {
    let a = cube_at(0.0, 0.0, 0.0);
    let ball = BRep::sphere(0.4, tol())
        .expect("sphere")
        .transform(&Transform3::translation(0.5, 0.5, 1.4))
        .expect("translates");
    let touches = contact_report(&a, &ball, &tol()).expect("reports");
    assert_eq!(touches.len(), 1);
    assert_eq!(touches[0].extent, TouchExtent::Point);
    assert!(a.faces.contains_key(touches[0].face_a));
    assert!(ball.faces.contains_key(touches[0].face_b));

    let flush = contact_report(&a, &cube_at(1.0, 0.0, 0.0), &tol()).expect("reports");
    assert!(flush.iter().any(|c| c.extent == TouchExtent::Area));

    let edge = contact_report(&a, &cube_at(1.0, 1.0, 0.0), &tol()).expect("reports");
    assert!(!edge.is_empty());
    assert!(
        edge.iter().all(|c| c.extent != TouchExtent::Area),
        "{edge:?}"
    );
}

/// Negative control: bodies one hundred tolerances apart do not touch.
#[test]
fn bodies_apart_report_no_contact() {
    let apart = contact_report(&cube_at(0.0, 0.0, 0.0), &cube_at(1.0001, 0.0, 0.0), &tol())
        .expect("reports");
    assert!(apart.is_empty(), "{apart:?}");
}

fn volume(b: &BRep) -> f64 {
    // Planar bodies: the divergence theorem over the faces' fan
    // triangles is exact.
    use popkern_core::geometry::Vec3;
    let mut v = 0.0;
    for (fh, _) in &b.faces {
        let ring = b.face_outer_loop_positions(fh);
        if ring.len() < 3 {
            continue;
        }
        let p0 = ring[0];
        let o = Vec3::new(p0.x, p0.y, p0.z);
        for w in ring[1..].windows(2) {
            let (p1, p2) = (
                Vec3::new(w[0].x, w[0].y, w[0].z),
                Vec3::new(w[1].x, w[1].y, w[1].z),
            );
            v += o.dot(&p1.cross(&p2)) / 6.0;
        }
    }
    v.abs()
}

#[test]
fn symmetric_difference_splits_into_the_two_differences() {
    let split = BoolOptions {
        non_manifold: NonManifoldPolicy::SplitToManifold,
    };
    let a = cube_at(0.0, 0.0, 0.0);
    let b = cube_at(0.5, 0.5, 0.5);
    assert!(matches!(
        symmetric_difference(a.clone(), b.clone(), tol(), &split),
        Err(crate::BoolError::SplitIntoBodies { bodies: 2, .. })
    ));
    let bodies = crate::symmetric_difference_bodies(a, b, tol(), &split).expect("serves split");
    assert_eq!(bodies.len(), 2);
    let mut total = 0.0;
    for body in &bodies {
        assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
        total += volume(body);
    }
    // vol(A ⊕ B) = vol(A ∪ B) − vol(A ∩ B) = 1.875 − 0.125.
    assert!((total - 1.75).abs() < 1e-9, "{total}");
}

#[test]
fn symmetric_difference_declines_by_default_and_disjoint_is_both() {
    let r = symmetric_difference(
        cube_at(0.0, 0.0, 0.0),
        cube_at(0.5, 0.5, 0.5),
        tol(),
        &BoolOptions::default(),
    );
    assert!(
        matches!(r, Err(BoolError::NonManifoldResult { .. })),
        "{r:?}"
    );
    let apart = symmetric_difference(
        cube_at(0.0, 0.0, 0.0),
        cube_at(3.0, 0.0, 0.0),
        tol(),
        &BoolOptions::default(),
    )
    .expect("disjoint serves");
    assert_eq!(apart.back_shell_count(), 2);
    assert!((volume(&apart) - 2.0).abs() < 1e-9);
}
