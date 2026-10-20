//! Plane sections against analytic loops (PK-856), and the plane split
//! built on the same half-space cutter.

use super::*;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn section(body: &BRep, origin: Point3, normal: Vec3) -> PlaneSection {
    let plane = SectionPlane::new(origin, normal).expect("plane");
    plane_section(body, &plane, 1e-4, &tol()).expect("section")
}

/// The centroid of the section's largest region, in 3-D.
fn centroid_3d(s: &PlaneSection) -> Point3 {
    let c = s.regions[0].centroid;
    s.plane.origin + s.plane.u_axis * c[0] + s.plane.v_axis() * c[1]
}

#[test]
fn a_cube_cut_square_on_and_across_its_diagonal() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    let s = section(&cube, Point3::new(0.0, 0.0, 0.5), Vec3::new(0.0, 0.0, 1.0));
    assert_eq!(s.regions.len(), 1);
    assert!((s.area() - 1.0).abs() < 1e-12, "{}", s.area());
    assert!(centroid_3d(&s).distance_to(&Point3::new(0.5, 0.5, 0.5)) < 1e-12);
    assert_eq!(s.regions[0].outer.points.len(), 4);
    // Through the centre normal to the diagonal: the regular hexagon of
    // side s = sqrt(2)/2, area (3 sqrt(3) / 2) s^2 = 3 sqrt(3) / 4.
    let s = section(&cube, Point3::new(0.5, 0.5, 0.5), Vec3::new(1.0, 1.0, 1.0));
    assert_eq!(s.regions[0].outer.points.len(), 6);
    let hexagon = 3.0 * 3.0_f64.sqrt() / 4.0;
    assert!(
        (s.area() - hexagon).abs() < 1e-12,
        "{} vs {hexagon}",
        s.area()
    );
    assert!(centroid_3d(&s).distance_to(&Point3::new(0.5, 0.5, 0.5)) < 1e-12);
}

#[test]
fn a_cylinder_cut_across_and_along_its_axis() {
    let (r, h) = (10.0, 4.0);
    let cyl = BRep::cylinder(r, h, tol()).expect("cylinder");
    let s = section(&cyl, Point3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 0.0, 1.0));
    let disc = std::f64::consts::PI * r * r;
    assert!(
        ((s.area() - disc) / disc).abs() < 1e-12,
        "{} vs {disc}",
        s.area()
    );
    assert!(centroid_3d(&s).distance_to(&Point3::new(0.0, 0.0, 1.0)) < 1e-9);
    // The loop is sampled within the chord: every point on the circle.
    for p in &s.regions[0].outer.points {
        let q = s.plane.origin + s.plane.u_axis * p[0] + s.plane.v_axis() * p[1];
        assert!(((q.x * q.x + q.y * q.y).sqrt() - r).abs() < 1e-9);
    }
    let s = section(&cyl, Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
    assert!((s.area() - 2.0 * r * h).abs() < 1e-9, "{}", s.area());
}

#[test]
fn a_bored_plate_section_keeps_its_hole_and_a_boss_its_disc() {
    let plate = BRep::prism(
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
            Point3::new(0.0, 10.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        2.0,
        tol(),
    )
    .expect("plate");
    let post = |z0: f64, height: f64| {
        let c = BRep::cylinder(2.0, height, tol()).expect("post");
        c.transform(&popkern_core::geometry::Transform3::translation(
            5.0, 5.0, z0,
        ))
        .expect("placed")
    };
    let bored = crate::subtraction_ref(&plate, &post(-1.0, 4.0), tol(), None).expect("bore");
    let s = section(&bored, Point3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 0.0, 1.0));
    assert_eq!(s.regions.len(), 1);
    assert_eq!(s.regions[0].holes.len(), 1);
    let ring = 100.0 - 4.0 * std::f64::consts::PI;
    assert!((s.area() - ring).abs() < 1e-9, "{} vs {ring}", s.area());
    let boss = crate::union_ref(&plate, &post(2.0, 4.0), tol(), None).expect("boss");
    let s = section(&boss, Point3::new(0.0, 0.0, 4.0), Vec3::new(0.0, 0.0, 1.0));
    assert_eq!(s.regions.len(), 1);
    let disc = 4.0 * std::f64::consts::PI;
    assert!((s.area() - disc).abs() < 1e-9, "{} vs {disc}", s.area());
    assert!(centroid_3d(&s).distance_to(&Point3::new(5.0, 5.0, 4.0)) < 1e-9);
    // A plane that misses the body has no section.
    assert!(
        section(&boss, Point3::new(0.0, 0.0, 9.0), Vec3::new(0.0, 0.0, 1.0))
            .regions
            .is_empty()
    );
}

#[test]
fn a_plane_split_halves_a_tilted_cylinder_larger_than_a_unit_cube() {
    use super::super::split_slice::{SplitPlane, split_by_plane};
    let cyl = BRep::cylinder(10.0, 4.0, tol()).expect("cylinder");
    let plane = SplitPlane {
        origin: Point3::new(1.0, 0.0, 2.0),
        normal: Vec3::new(1.0, 0.0, 1.0),
    };
    let halves = split_by_plane(&cyl, &plane, &tol()).expect("split");
    let n = Vec3::new(1.0, 0.0, 1.0).normalized().expect("unit");
    for (half, sign) in [(&halves.positive, 1.0), (&halves.negative, -1.0)] {
        let half = half.as_ref().expect("both halves exist");
        assert!(half.validate().is_ok(), "{:?}", half.validate().err());
        // Every vertex on its own side of the plane.
        for v in half.vertices.values() {
            assert!(sign * (v.position - plane.origin).dot(&n) > -1e-9);
        }
    }
}
