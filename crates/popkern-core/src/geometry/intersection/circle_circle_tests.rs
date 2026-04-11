use super::*;
use std::f64::consts::TAU;

fn cc(
    c1: Point3,
    n1: Vec3,
    r1: f64,
    a1: (f64, f64),
    c2: Point3,
    n2: Vec3,
    r2: f64,
    a2: (f64, f64),
) -> Vec<CurveIntersection> {
    intersect_circle_circle(&c1, &n1, r1, a1, &c2, &n2, r2, a2, &Tolerance::default())
}
fn zn() -> Vec3 {
    Vec3::new(0.0, 0.0, 1.0)
}
fn o() -> Point3 {
    Point3::new(0.0, 0.0, 0.0)
}
fn full() -> (f64, f64) {
    (0.0, TAU)
}

#[test]
fn coplanar_cases() {
    // Two points
    let r = cc(
        o(),
        zn(),
        1.0,
        full(),
        Point3::new(1.0, 0.0, 0.0),
        zn(),
        1.0,
        full(),
    );
    assert_eq!(r.len(), 2);
    // Tangent
    assert_eq!(
        cc(
            o(),
            zn(),
            1.0,
            full(),
            Point3::new(2.0, 0.0, 0.0),
            zn(),
            1.0,
            full()
        )
        .len(),
        1
    );
    // Disjoint
    assert!(
        cc(
            o(),
            zn(),
            1.0,
            full(),
            Point3::new(5.0, 0.0, 0.0),
            zn(),
            1.0,
            full()
        )
        .is_empty()
    );
    // Identical => overlap
    let r = cc(o(), zn(), 1.0, full(), o(), zn(), 1.0, full());
    assert!(matches!(&r[0], CurveIntersection::Overlap { .. }));
    // Concentric different radii
    assert!(cc(o(), zn(), 1.0, full(), o(), zn(), 2.0, full()).is_empty());
}

#[test]
fn non_coplanar_cases() {
    let xn = Vec3::new(1.0, 0.0, 0.0);
    assert_eq!(cc(o(), zn(), 1.0, full(), o(), xn, 1.0, full()).len(), 2);
    // Disjoint non-coplanar
    assert!(
        cc(
            o(),
            zn(),
            1.0,
            full(),
            Point3::new(10.0, 10.0, 10.0),
            xn,
            1.0,
            full()
        )
        .is_empty()
    );
    // Restricted arc range
    assert!(cc(o(), zn(), 1.0, (0.0, 1.0), o(), xn, 1.0, full()).len() <= 2);
}

#[test]
fn degenerate_and_offset_cases() {
    let z = Vec3::new(0.0, 0.0, 0.0);
    let pz = Point3::new(0.0, 0.0, 5.0);
    assert!(cc(o(), z, 1.0, full(), o(), zn(), 1.0, full()).is_empty());
    assert!(cc(o(), zn(), 1.0, full(), o(), z, 1.0, full()).is_empty());
    assert!(cc(o(), zn(), 1.0, full(), pz, zn(), 1.0, full()).is_empty());
}

#[test]
fn angle_in_range_cases() {
    assert!(angle_in_range(-0.1, (0.0, TAU), 0.2));
    assert!(angle_in_range(0.1, (-0.1, 1.0), 0.01));
    // Wrapping range
    assert!(angle_in_range(5.5, (5.0, 1.0), 0.01));
    assert!(angle_in_range(0.5, (5.0, 1.0), 0.01));
    assert!(!angle_in_range(3.0, (5.0, 1.0), 0.01));
}
