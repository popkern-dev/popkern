use super::*;
use crate::geometry::Point3;
use crate::geometry::surface::SurfaceGeometry;

const SCALE: ContactScale = ContactScale {
    linear: 1e-6,
    extent: 10.0,
};

fn sphere(c: (f64, f64, f64), r: f64) -> SurfaceGeometry {
    SurfaceGeometry::Sphere {
        center: Point3::new(c.0, c.1, c.2),
        radius: r,
        axis: Vec3::new(0.0, 0.0, 1.0),
    }
}

fn cylinder_z(origin: (f64, f64, f64), r: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(origin.0, origin.1, origin.2),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: r,
    }
}

fn plane_z(z: f64, up: bool) -> SurfaceGeometry {
    let s = if up { 1.0 } else { -1.0 };
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, z),
        normal: Vec3::new(0.0, 0.0, s),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, s, 0.0),
    }
}

fn uv_of(s: &SurfaceGeometry, p: Point3) -> (f64, f64) {
    s.closest_point(&p).expect("closest point")
}

fn classify(s1: &SurfaceGeometry, s2: &SurfaceGeometry, p: Point3) -> ContactClassification {
    classify_contact(s1, uv_of(s1, p), s2, uv_of(s2, p), SCALE).expect("classifies")
}

/// Whether unit `d` is parallel to the direction `(x, y, z)`.
fn along(d: Vec3, x: f64, y: f64, z: f64) -> bool {
    let w = Vec3::new(x, y, z).normalized().expect("direction");
    d.cross(&w).length() < 1e-6
}

/// Viviani's singular point: a ball of radius 10 and a cylinder of
/// radius 5 on x = 5 are tangent at (10, 0, 0); their relative
/// curvature is +0.1 round the cylinder and −0.1 along its axis, so two
/// branches cross there at ±45° (the figure eight's tangents
/// (0, 10, ±10) at t = 0 and π).
#[test]
fn viviani_point_is_a_crossing_at_45_degrees() {
    let c = classify(
        &sphere((0.0, 0.0, 0.0), 10.0),
        &cylinder_z((5.0, 0.0, 0.0), 5.0),
        Point3::new(10.0, 0.0, 0.0),
    );
    let [lo, hi] = c.eigenvalues.expect("regular point");
    assert!(
        (lo + 0.1).abs() < 1e-9 && (hi - 0.1).abs() < 1e-9,
        "{lo} {hi}"
    );
    let SurfaceContactKind::CrossingPoint { branches } = &c.kind else {
        panic!("{:?}", c.kind);
    };
    assert_eq!(branches.len(), 2);
    assert!(branches.iter().any(|&b| along(b, 0.0, 1.0, 1.0)));
    assert!(branches.iter().any(|&b| along(b, 0.0, 1.0, -1.0)));
}

/// Negative control: the first classifier subtracted each surface's
/// second fundamental form in its OWN parameter basis. At Viviani's
/// point those coefficients are (L, N) = (−10, −10) on the sphere and
/// (−5, 0) on the cylinder: the difference (−5, −10) is definite, which
/// reads as a touch where the surfaces cross.
#[test]
fn the_parameter_basis_difference_misreads_viviani_as_a_touch() {
    let s1 = sphere((0.0, 0.0, 0.0), 10.0);
    let s2 = cylinder_z((5.0, 0.0, 0.0), 5.0);
    let p = Point3::new(10.0, 0.0, 0.0);
    let n = Vec3::new(1.0, 0.0, 0.0);
    let raw = |s: &SurfaceGeometry| {
        let (u, v) = uv_of(s, p);
        (
            s.second_partial_uu(u, v).dot(&n),
            s.second_partial_uv(u, v).dot(&n),
            s.second_partial_vv(u, v).dot(&n),
        )
    };
    let (a, b) = (raw(&s1), raw(&s2));
    let (dl, dm, dn) = (a.0 - b.0, a.1 - b.1, a.2 - b.2);
    assert!(dl * dn - dm * dm > 0.0, "the raw difference reads definite");
    assert!(matches!(
        classify(&s1, &s2, p).kind,
        SurfaceContactKind::CrossingPoint { .. }
    ));
}

/// The side a surface lying toward `toward` from the contact is on,
/// against the normal the classification measured along (a sphere's
/// evaluated normal points inward; the face's orientation flag, applied
/// by callers, makes it outward).
fn side_toward(c: &ContactClassification, toward: Vec3) -> ContactSide {
    if c.normal.dot(&toward) > 0.0 {
        ContactSide::AlongNormal
    } else {
        ContactSide::AgainstNormal
    }
}

/// A cube's inscribed ball touches its bottom face (z = 0) at its own
/// pole, a chart singularity of the sphere: an isolated touch, the ball
/// on the plane's +z side, whichever surface comes first.
#[test]
fn inscribed_ball_touches_the_face_at_its_pole() {
    let ball = sphere((5.0, 5.0, 5.0), 5.0);
    let p = Point3::new(5.0, 5.0, 0.0);
    let up = Vec3::new(0.0, 0.0, 1.0);
    for plane in [plane_z(0.0, true), plane_z(0.0, false)] {
        let c = classify(&plane, &ball, p);
        assert_eq!(
            c.kind,
            SurfaceContactKind::TangentPoint {
                side: side_toward(&c, up)
            }
        );
    }
    // From the ball, the plane lies below the contact.
    let c = classify(&ball, &plane_z(0.0, true), p);
    assert_eq!(
        c.kind,
        SurfaceContactKind::TangentPoint {
            side: side_toward(&c, up * -1.0)
        }
    );
}

/// A ball of radius 5 in a cylinder of radius 5 touches it along the
/// equator: rank-one relative curvature, a tangent curve along the
/// equator with the ball inside the cylinder.
#[test]
fn ball_in_band_is_tangent_along_the_equator() {
    let band = cylinder_z((0.0, 0.0, 0.0), 5.0);
    let ball = sphere((0.0, 0.0, 0.0), 5.0);
    let c = classify(&band, &ball, Point3::new(5.0, 0.0, 0.0));
    let SurfaceContactKind::TangentCurve { direction, side } = c.kind else {
        panic!("{:?}", c.kind);
    };
    assert!(along(direction, 0.0, 1.0, 0.0), "{direction:?}");
    // The ball is inside the band: toward the axis.
    assert_eq!(side, Some(side_toward(&c, Vec3::new(-1.0, 0.0, 0.0))));
    assert!(c.ring.is_some());
}

/// The ice-cream cone: a cone of half-angle atan(3/4) with its apex at
/// (0, 0, −25) is tangent to the ball of radius 15 along the circle of
/// radius 12 at z = −9: a tangent curve along the circle, the ball
/// outside the cone below the circle's plane and inside it above.
#[test]
fn ice_cream_cone_is_tangent_along_its_rim() {
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, -25.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: 0.75f64.patan(),
    };
    let ball = sphere((0.0, 0.0, 0.0), 15.0);
    let c = classify(&cone, &ball, Point3::new(12.0, 0.0, -9.0));
    let SurfaceContactKind::TangentCurve { direction, side } = c.kind else {
        panic!("{:?}", c.kind);
    };
    assert!(along(direction, 0.0, 1.0, 0.0), "{direction:?}");
    // The ball lies inside the cone, toward its axis, on both sides
    // of the circle.
    assert_eq!(side, Some(side_toward(&c, Vec3::new(-1.0, 0.0, 0.0))));
}

/// Two copies of one plane coincide.
#[test]
fn coincident_planes() {
    let c = classify(
        &plane_z(0.0, true),
        &plane_z(0.0, true),
        Point3::new(1.0, 2.0, 0.0),
    );
    assert_eq!(c.kind, SurfaceContactKind::Coincident);
}

/// Balls of equal radius touching from outside: an isolated touch,
/// the second beyond the contact in +x.
#[test]
fn equal_balls_touching_from_outside() {
    let a = sphere((0.0, 0.0, 0.0), 5.0);
    let b = sphere((10.0, 0.0, 0.0), 5.0);
    let c = classify(&a, &b, Point3::new(5.0, 0.0, 0.0));
    assert_eq!(
        c.kind,
        SurfaceContactKind::TangentPoint {
            side: side_toward(&c, Vec3::new(1.0, 0.0, 0.0))
        }
    );
}

/// Surfaces crossing transversally are no contact.
#[test]
fn transversal_surfaces_are_refused() {
    let s = sphere((0.0, 0.0, 0.0), 5.0);
    let p = plane_z(0.0, true);
    let q = Point3::new(5.0, 0.0, 0.0);
    assert!(classify_contact(&s, uv_of(&s, q), &p, uv_of(&p, q), SCALE).is_err());
}

/// The floor is 2ε / L².
#[test]
fn the_curvature_floor_is_the_tolerance_over_the_extent() {
    assert!((SCALE.curvature_floor() - 2e-8).abs() < 1e-20);
}

/// Balls of radius 5 and 5 + 1e-9 tangent from inside differ in
/// curvature by 4e-11, under the floor 2ε / L²: second order cannot
/// tell, and the ring sees them within the tolerance all round:
/// coincident at this scale. Radii 5 and 6 differ by 1/30: a touch.
#[test]
fn nearly_equal_inner_balls_coincide_and_unequal_ones_touch() {
    let p = Point3::new(5.0, 0.0, 0.0);
    let a = sphere((0.0, 0.0, 0.0), 5.0);
    let near = sphere((-1e-9, 0.0, 0.0), 5.0 + 1e-9);
    let c = classify(&a, &near, p);
    let [lo, hi] = c.eigenvalues.expect("regular");
    assert!(lo.abs().max(hi.abs()) < c.floor, "{lo} {hi}");
    assert_eq!(c.kind, SurfaceContactKind::Coincident);
    assert!(c.ring.is_some());
    let far = sphere((-1.0, 0.0, 0.0), 6.0);
    let c = classify(&a, &far, p);
    assert!(matches!(c.kind, SurfaceContactKind::TangentPoint { .. }));
    assert!(c.ring.is_none());
}
