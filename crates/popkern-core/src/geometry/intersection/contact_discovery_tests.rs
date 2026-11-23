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

/// The cone touching a cylinder at scale 1000 (below).
fn cone_on_cylinder_at_scale_1000() -> (SurfaceGeometry, SurfaceGeometry, Point3) {
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, -2000.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1000.0,
    };
    // The cone's outward normal at the touch is -x, its generator g
    // there square to x: the axis is cos α g + sin α x.
    let (alpha, beta) = (0.3_f64, 0.7_f64);
    let g = Vec3::new(0.0, beta.sin(), beta.cos());
    let x = Vec3::new(1.0, 0.0, 0.0);
    let touch = Point3::new(1000.0, 0.0, 0.0);
    let cone = SurfaceGeometry::Cone {
        apex: touch - g * 1500.0,
        axis: g * alpha.cos() + x * alpha.sin(),
        half_angle: alpha,
    };
    (cylinder, cone, touch)
}

/// A cylinder of radius 1000 and a cone touching it outside at
/// (1000, 0, 0), the cone's generator there skew to the cylinder's axis
/// (scale 1000). Each chart pairs an angle with a length of thousands:
/// one difference step for both (a sixth of a radian in the angle)
/// differenced an unsymmetric gap over a wide arc, and Newton stalled
/// short of the touch, which read as a near miss (corpus YC-E 943).
#[test]
fn a_cone_touching_a_cylinder_at_scale_1000_is_found_at_the_exact_point() {
    let (cylinder, cone, touch) = cone_on_cylinder_at_scale_1000();
    let tau = std::f64::consts::TAU;
    let scale = ContactScale {
        linear: 1e-6,
        extent: 5000.0,
    };
    let found = find_surface_contacts(
        &cylinder,
        ((0.0, tau), (0.0, 4000.0)),
        &cone,
        ((0.0, tau), (0.0, 4000.0)),
        scale,
    );
    assert_eq!(found.contacts.len(), 1, "{found:?}");
    let off = found.contacts[0].point.distance_to(&touch);
    assert!(off < 1e-3, "the touch found {off:.3e} from the exact point");
}

/// From the cylinder's side alone, seeded a grid cell off the touch
/// (an angle of 0.13 and 250 along the axis): Newton's steps along the
/// gap's shallow valley were measured against a gradient in parameter
/// units, where the angle's component (at radius 1000) outweighs the
/// length's a thousandfold; each read as worse, the damping crawled,
/// and the walk stopped with the gap 2.5.
#[test]
fn newton_from_the_cylinder_side_reaches_the_cone_touch() {
    let (cylinder, cone, touch) = cone_on_cylinder_at_scale_1000();
    let tau = std::f64::consts::TAU;
    let window = ((0.0, tau), (0.0, 4000.0));
    let scale = ContactScale {
        linear: 1e-6,
        extent: 5000.0,
    };
    let found = contact_from(
        (&cylinder, window),
        (&cone, window),
        (0.13, 2250.0),
        0.17,
        scale,
    )
    .expect("the touch from the cylinder's side");
    let off = found.point.distance_to(&touch);
    assert!(off < 1e-3, "the touch found {off:.3e} from the exact point");
}

/// A cylinder of radius 0.001 and a cone crossing it tangentially at
/// scale 0.001 (corpus YC-V 966, its faces' surfaces and windows), Newton
/// seeded on a branch 0.013 rad round the cylinder from the crossing with
/// a march's difference step (a millionth of the window). Damped in
/// parameter units, where the cylinder's curvature per radian² and per
/// length² stand a millionfold apart, a millionth of the larger
/// cancelled the saddle's curvature along the angle: the step overshot
/// threefold and then turned the wrong way, Newton stopped on the branch,
/// and the march crawled on past the crossing for minutes.
#[test]
fn newton_reaches_a_tangent_crossing_at_scale_a_thousandth() {
    let cylinder = SurfaceGeometry::Cylinder {
        origin: Point3::new(
            -0.00031317792761719317,
            0.0018470743392259641,
            0.0018623162365610452,
        ),
        axis: Vec3::new(
            0.019997835712565692,
            -0.8020045909915312,
            0.5969830170074524,
        ),
        radius: 0.001,
    };
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(
            0.0019066377683225639,
            0.0009372593284713874,
            0.000538019358331209,
        ),
        axis: Vec3::new(-0.7188112200556341, 0.18635206449919964, 0.6697636433690745),
        half_angle: 0.4636476090008061,
    };
    let tau = std::f64::consts::TAU;
    let w1 = ((0.0, tau), (-0.00991489975031797, 0.012414899750317971));
    let w2 = ((0.0, tau), (-0.009428289554577885, 0.012028289554577885));
    let scale = ContactScale {
        linear: 1e-6,
        extent: 0.00499873694343911,
    };
    let crossing = Point3::new(
        -0.0005492838599147503,
        0.0007184050674709597,
        0.0014953458170282213,
    );
    for seed in [
        (4.473026747641222, 0.0006618649879867174),
        (4.4678942441946425, 0.000669348027939345),
        (4.427747644569971, 0.0007291146903645988),
    ] {
        let found = contact_from((&cylinder, w1), (&cone, w2), seed, 1e-6 * tau, scale)
            .expect("the crossing");
        assert!(
            matches!(
                found.classification.kind,
                SurfaceContactKind::CrossingPoint { .. }
            ),
            "{:?}",
            found.classification.kind
        );
        let off = found.point.distance_to(&crossing);
        assert!(off < 1e-9, "the crossing found {off:.3e} from it");
    }
}
