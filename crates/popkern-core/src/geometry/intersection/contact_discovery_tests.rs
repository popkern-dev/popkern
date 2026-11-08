use super::*;
use crate::geometry::Vec3;
use crate::geometry::intersection::contact_classify::SurfaceContactKind;

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

fn plane_z(z: f64) -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, z),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

fn whole(s: &SurfaceGeometry) -> ParamWindow {
    s.parameter_range()
}

/// The cube's bottom face against its inscribed ball: one isolated
/// touch at the face's centre, found from the plane's side and from
/// the ball's (where it is the ball's pole).
#[test]
fn inscribed_ball_touches_the_face_once() {
    let ball = sphere((5.0, 5.0, 5.0), 5.0);
    let face = plane_z(0.0);
    let face_window = ((0.0, 10.0), (0.0, 10.0));
    let found = find_surface_contacts(&face, face_window, &ball, whole(&ball), SCALE);
    assert_eq!(found.contacts.len(), 1, "{found:?}");
    let c = &found.contacts[0];
    assert!(
        c.point.distance_to(&Point3::new(5.0, 5.0, 0.0)) < 1e-6,
        "{:?}",
        c.point
    );
    assert!(
        matches!(
            c.classification.kind,
            SurfaceContactKind::TangentPoint { .. }
        ),
        "{c:?}"
    );
    let back = find_surface_contacts(&ball, whole(&ball), &face, face_window, SCALE);
    assert_eq!(back.contacts.len(), 1, "{back:?}");
    assert!(
        back.contacts[0]
            .point
            .distance_to(&Point3::new(5.0, 5.0, 0.0))
            < 1e-6
    );
}

/// Viviani: the ball and the cylinder meet in a figure eight; the one
/// contact is its crossing at (10, 0, 0), the transversal curve is not
/// rediscovered as contacts.
#[test]
fn viviani_has_one_contact_at_the_crossing() {
    let ball = sphere((0.0, 0.0, 0.0), 10.0);
    let rod = SurfaceGeometry::Cylinder {
        origin: Point3::new(5.0, 0.0, -12.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 5.0,
    };
    let rod_window = ((0.0, std::f64::consts::TAU), (0.0, 24.0));
    let found = find_surface_contacts(&rod, rod_window, &ball, whole(&ball), SCALE);
    assert_eq!(found.contacts.len(), 1, "{found:?}");
    let c = &found.contacts[0];
    assert!(
        c.point.distance_to(&Point3::new(10.0, 0.0, 0.0)) < 1e-5,
        "{:?}",
        c.point
    );
    assert!(matches!(
        c.classification.kind,
        SurfaceContactKind::CrossingPoint { .. }
    ));
}

/// A plane through a ball's middle crosses it: no contact.
#[test]
fn a_crossing_plane_is_no_contact() {
    let ball = sphere((0.0, 0.0, 0.0), 5.0);
    let face = plane_z(0.0);
    let found = find_surface_contacts(
        &face,
        ((-10.0, 10.0), (-10.0, 10.0)),
        &ball,
        whole(&ball),
        SCALE,
    );
    assert!(found.contacts.is_empty(), "{found:?}");
}

/// A ball 2e-5 above a plane is a near miss, not a contact.
#[test]
fn a_ball_just_above_the_plane_is_a_near_miss() {
    let ball = sphere((5.0, 5.0, 5.0 + 2e-5), 5.0);
    let face = plane_z(0.0);
    let found = find_surface_contacts(
        &face,
        ((0.0, 10.0), (0.0, 10.0)),
        &ball,
        whole(&ball),
        SCALE,
    );
    assert!(found.contacts.is_empty(), "{found:?}");
    assert_eq!(found.near_misses.len(), 1, "{found:?}");
    assert!((found.near_misses[0].gap.abs() - 2e-5).abs() < 1e-7);
}

/// A ball inside a cylinder of its own radius touches it along the
/// equator: every contact found is a point of that tangent curve.
#[test]
fn ball_in_band_contacts_lie_on_the_equator() {
    let ball = sphere((0.0, 0.0, 0.0), 5.0);
    let band = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, -3.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 5.0,
    };
    let found = find_surface_contacts(
        &band,
        ((0.0, std::f64::consts::TAU), (0.0, 6.0)),
        &ball,
        whole(&ball),
        SCALE,
    );
    assert!(!found.contacts.is_empty());
    for c in &found.contacts {
        assert!(c.point.z.abs() < 1e-5, "{:?}", c.point);
        assert!(
            matches!(
                c.classification.kind,
                SurfaceContactKind::TangentCurve { .. }
            ),
            "{:?}",
            c.classification.kind
        );
    }
}
