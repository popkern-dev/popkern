//! NURBS twins (PK-1126, the plan's §8.3): each contact of the
//! tangency corpus found on the analytic pair is found again with the
//! first surface converted exactly to NURBS, where the search runs on a
//! spline rather than any closed form; the two agree in count, place
//! and kind. `to_nurbs_in_range` converts each exactly (a plane's
//! bilinear patch; cylinders, cones and spheres as rational surfaces of
//! revolution in their chart frames since PK-1142, whose sampled
//! control net left the surface and these twins found nothing).

use super::*;
use crate::geometry::Vec3;
use crate::geometry::intersection::contact_classify::SurfaceContactKind;
use std::f64::consts::TAU;

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

/// The kind's name, for comparing two classifications.
fn kind_of(c: &SurfaceContact) -> &'static str {
    match c.classification.kind {
        SurfaceContactKind::TangentPoint { .. } => "touch",
        SurfaceContactKind::CrossingPoint { .. } => "crossing",
        SurfaceContactKind::TangentCurve { .. } => "tangent curve",
        _ => "other",
    }
}

/// Contacts of `(s1 over w1, s2)` and of its NURBS twin agree: every
/// twin contact lies within 1e-5 of an analytic one of the same kind
/// and every analytic one of a twin one. Tangent curves are sampled
/// at different points, so for them only the kinds and that each point
/// lies on the other's curve (`on_curve`) are compared.
fn assert_twins_agree(
    s1: &SurfaceGeometry,
    w1: ParamWindow,
    s2: &SurfaceGeometry,
    on_curve: impl Fn(&Point3) -> bool,
) {
    let w2 = s2.parameter_range();
    let analytic = find_surface_contacts(s1, w1, s2, w2, SCALE);
    let twin = SurfaceGeometry::Nurbs(s1.to_nurbs_in_range(w1.0, w1.1));
    let tw = twin.parameter_range();
    let nurbs = find_surface_contacts(&twin, tw, s2, w2, SCALE);
    assert!(
        !analytic.contacts.is_empty(),
        "the analytic pair finds its contact"
    );
    assert!(!nurbs.contacts.is_empty(), "the twin finds it too");
    let matched = |from: &[SurfaceContact], to: &[SurfaceContact]| {
        from.iter().all(|a| {
            if kind_of(a) == "tangent curve" {
                return on_curve(&a.point) && to.iter().any(|b| kind_of(b) == "tangent curve");
            }
            to.iter()
                .any(|b| kind_of(b) == kind_of(a) && a.point.distance_to(&b.point) < 1e-5)
        })
    };
    assert!(
        matched(&analytic.contacts, &nurbs.contacts),
        "analytic {:?} vs twin {:?}",
        analytic.contacts,
        nurbs.contacts
    );
    assert!(
        matched(&nurbs.contacts, &analytic.contacts),
        "twin {:?} vs analytic {:?}",
        nurbs.contacts,
        analytic.contacts
    );
}

/// The inscribed ball's touch at a cube face's centre, from the face
/// (a bilinear twin) and from the ball (a rational twin, the touch at
/// its pole).
#[test]
fn the_inscribed_ball_touch_is_found_on_both_twins() {
    let ball = sphere((5.0, 5.0, 5.0), 5.0);
    let face = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert_twins_agree(&face, ((0.0, 10.0), (0.0, 10.0)), &ball, |_| false);
    assert_twins_agree(&ball, ball.parameter_range(), &face, |_| false);
}

/// Viviani's crossing at (10, 0, 0), from the cylinder's twin.
#[test]
fn vivianis_crossing_is_found_on_the_cylinders_twin() {
    let ball = sphere((0.0, 0.0, 0.0), 10.0);
    let rod = SurfaceGeometry::Cylinder {
        origin: Point3::new(5.0, 0.0, -12.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 5.0,
    };
    assert_twins_agree(&rod, ((0.0, TAU), (0.0, 24.0)), &ball, |_| false);
}

/// The ball in its band, tangent along the equator, from the band's
/// twin: tangent-curve points on z = 0 either way.
#[test]
fn the_ball_in_its_band_is_tangent_along_the_equator_on_the_twin() {
    let ball = sphere((0.0, 0.0, 0.0), 5.0);
    let band = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, -3.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 5.0,
    };
    assert_twins_agree(&band, ((0.0, TAU), (0.0, 6.0)), &ball, |p| p.z.abs() < 1e-5);
}

/// The ice-cream cone tangent to its ball along the rim circle (z = -9,
/// radius 12), from the cone's twin.
#[test]
fn the_ice_cream_cone_is_tangent_along_its_rim_on_the_twin() {
    let ball = sphere((0.0, 0.0, 0.0), 15.0);
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, -25.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: 0.75f64.atan(),
    };
    assert_twins_agree(&cone, ((0.0, TAU), (0.0, 18.0)), &ball, |p| {
        (p.z + 9.0).abs() < 1e-4 && (p.x.hypot(p.y) - 12.0).abs() < 1e-4
    });
}

/// The cylinder inscribed through a cube, tangent to a face along a
/// ruling (z = -5 under the axis), from the face's twin and the
/// cylinder's: tangent-curve points on that line either way.
#[test]
fn the_inscribed_cylinder_is_tangent_along_a_ruling_on_both_twins() {
    let rod = SurfaceGeometry::Cylinder {
        origin: Point3::new(-10.0, 0.0, 0.0),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: 5.0,
    };
    let floor = SurfaceGeometry::Plane {
        origin: Point3::new(-5.0, -5.0, -5.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let on_ruling = |p: &Point3| p.y.abs() < 1e-4 && (p.z + 5.0).abs() < 1e-4;
    assert_twins_agree(&floor, ((0.0, 10.0), (0.0, 10.0)), &rod, on_ruling);
    assert_twins_agree(&rod, ((0.0, TAU), (5.0, 15.0)), &floor, on_ruling);
}

/// Villarceau's bitangent plane 3x - 4z = 0 through a torus (5, 3):
/// its two circles cross at (±16/5, 0, ±12/5), from the plane's twin.
#[test]
fn villarceaus_crossings_are_found_on_the_planes_twin() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 3.0,
    };
    let n = Vec3::new(3.0, 0.0, -4.0) * 0.2;
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0)
            - (Vec3::new(4.0, 0.0, 3.0) * 0.2) * 9.0
            - Vec3::new(0.0, 1.0, 0.0) * 9.0,
        normal: n,
        u_axis: Vec3::new(4.0, 0.0, 3.0) * 0.2,
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert_twins_agree(&plane, ((0.0, 18.0), (0.0, 18.0)), &torus, |_| false);
    let analytic = find_surface_contacts(
        &plane,
        ((0.0, 18.0), (0.0, 18.0)),
        &torus,
        torus.parameter_range(),
        SCALE,
    );
    for want in [Point3::new(3.2, 0.0, 2.4), Point3::new(-3.2, 0.0, -2.4)] {
        assert!(
            analytic
                .contacts
                .iter()
                .any(|c| c.point.distance_to(&want) < 1e-5),
            "{want:?} in {:?}",
            analytic.contacts
        );
    }
}

/// A torus (5, 1) inside a drum of radius 6, tangent along its outer
/// equator, from the drum's twin.
#[test]
fn the_torus_in_its_drum_is_tangent_along_its_equator_on_the_twin() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 1.0,
    };
    let drum = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, -3.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 6.0,
    };
    assert_twins_agree(&drum, ((0.0, TAU), (0.0, 6.0)), &torus, |p| {
        p.z.abs() < 1e-4 && (p.x.hypot(p.y) - 6.0).abs() < 1e-4
    });
}

/// Two balls of radius 3 touching at the origin, from one's twin.
#[test]
fn touching_balls_are_found_on_a_twin() {
    let a = sphere((-3.0, 0.0, 0.0), 3.0);
    let b = sphere((3.0, 0.0, 0.0), 3.0);
    assert_twins_agree(&a, a.parameter_range(), &b, |_| false);
}

/// The paraboloid z = x²/2 (a quadratic Bézier profile revolved, its
/// tip's curvature radius 1) touching the plane under its tip: the
/// normal-cone subdivision finds the touch from the paraboloid's side
/// and classifies it, no box left undecided (the plan's E4). A ball of
/// radius 1 in the cup (a fourth-order touch) and the plane searched
/// from the plane's side are filed with their evidence (PK-1179).
#[test]
fn the_paraboloid_tip_on_a_plane_is_decided() {
    use crate::geometry::NurbsCurve;
    let profile = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 2.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let cup = SurfaceGeometry::Nurbs(
        crate::geometry::revolution::revolution_surface(
            &profile,
            &Point3::new(0.0, 0.0, 0.0),
            &Vec3::new(0.0, 0.0, 1.0),
            TAU,
        )
        .expect("paraboloid"),
    );
    let floor = SurfaceGeometry::Plane {
        origin: Point3::new(-3.0, -3.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let found = find_surface_contacts(
        &cup,
        cup.parameter_range(),
        &floor,
        ((0.0, 6.0), (0.0, 6.0)),
        SCALE,
    );
    assert_eq!(found.contacts.len(), 1, "{:?}", found.contacts);
    let c = &found.contacts[0];
    assert!(
        c.point.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1e-6,
        "{:?}",
        c.point
    );
    assert_eq!(kind_of(c), "touch", "{:?}", c.classification);
}
