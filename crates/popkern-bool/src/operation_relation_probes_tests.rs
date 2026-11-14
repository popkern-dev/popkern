use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::*;
use crate::error::BoolOperation;

/// The cube [0, 10]³.
fn cube10() -> BRep {
    BRep::unit_cube(Tolerance::default())
        .expect("cube builds")
        .transform(&Transform3::scale_uniform(10.0))
        .expect("scales")
}

/// The ball of radius 5 inscribed in [`cube10`], turned by `turn`
/// about its centre first: unturned, its poles lie on the cube's top
/// and bottom faces and its seam runs through the side contacts.
fn inscribed_ball(turn: Option<([f64; 3], f64)>) -> BRep {
    let mut ball = BRep::sphere(5.0, Tolerance::default()).expect("sphere builds");
    if let Some((axis, angle)) = turn {
        let axis = Vec3::new(axis[0], axis[1], axis[2]);
        let turn = Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), axis, angle)
            .expect("axis");
        ball = ball.transform(&turn).expect("turns");
    }
    ball.transform(&Transform3::translation(5.0, 5.0, 5.0))
        .expect("translates")
}

fn surfaces(body: &BRep) -> Vec<&SurfaceGeometry> {
    body.faces
        .values()
        .filter_map(|f| f.surface.as_ref().map(|s| &s.geometry))
        .collect()
}

const TURNS: [Option<([f64; 3], f64)>; 3] = [
    None,
    Some(([1.0, 0.7, 0.3], 0.645_771_823_237_9)),
    Some(([0.0, 1.0, 0.0], std::f64::consts::FRAC_PI_2)),
];

/// The common part of a cube and its inscribed ball is the ball, at
/// every placement of the ball's poles and seam. Unturned it served an
/// EMPTY body: every probe of the ball (its pole vertices, seam
/// samples) lay on a cube face, and the OnBoundary reverse read as
/// "apart" (PK-1124).
#[test]
fn inscribed_ball_common_is_the_ball_at_every_placement() {
    for turn in TURNS {
        let r = crate::boolean(
            BoolOperation::Intersection,
            cube10(),
            inscribed_ball(turn),
            Tolerance::default(),
            None,
        )
        .unwrap_or_else(|e| panic!("{turn:?}: {e}"));
        let s = surfaces(&r);
        assert!(!s.is_empty(), "{turn:?}: the common part came back empty");
        assert!(
            s.iter()
                .all(|g| matches!(g, SurfaceGeometry::Sphere { .. })),
            "{turn:?}: the common part is not the ball"
        );
    }
}

/// The union of a cube and its inscribed ball is the cube, six planar
/// faces, at every placement. Unturned it declined with six
/// zero-thickness contacts, the operands combined side by side.
#[test]
fn inscribed_ball_union_is_the_cube_at_every_placement() {
    for turn in TURNS {
        let r = crate::boolean(
            BoolOperation::Union,
            cube10(),
            inscribed_ball(turn),
            Tolerance::default(),
            None,
        )
        .unwrap_or_else(|e| panic!("{turn:?}: {e}"));
        let s = surfaces(&r);
        assert_eq!(s.len(), 6, "{turn:?}");
        assert!(s.iter().all(|g| matches!(g, SurfaceGeometry::Plane { .. })));
    }
}

/// Cube minus its inscribed ball leaves a wall of no thickness at the
/// six face centres, where the cavity touches the outside: a
/// non-manifold result, declined under the default policy naming
/// exactly those six points.
#[test]
fn inscribed_ball_cut_declines_at_the_six_face_centres() {
    let e = crate::boolean(
        BoolOperation::Subtraction,
        cube10(),
        inscribed_ball(None),
        Tolerance::default(),
        None,
    )
    .expect_err("non-manifold");
    let BoolError::NonManifoldResult { contacts } = e else {
        panic!("{e}");
    };
    let centres = [
        (5.0, 5.0, 0.0),
        (5.0, 5.0, 10.0),
        (5.0, 0.0, 5.0),
        (5.0, 10.0, 5.0),
        (0.0, 5.0, 5.0),
        (10.0, 5.0, 5.0),
    ];
    assert_eq!(contacts.len(), 6);
    for (x, y, z) in centres {
        let c = Point3::new(x, y, z);
        assert!(
            contacts.iter().any(|k| k.location.distance_to(&c) < 1e-6),
            "no contact at {c:?}"
        );
    }
}

/// A ball resting on the cube's top face from outside touches it at
/// one point: apart. The common part is empty (regularised), which is
/// right here and must stay so.
#[test]
fn a_ball_resting_on_the_cube_is_apart() {
    let ball = BRep::sphere(5.0, Tolerance::default())
        .expect("sphere builds")
        .transform(&Transform3::translation(5.0, 5.0, 15.0))
        .expect("translates");
    let r = crate::boolean(
        BoolOperation::Intersection,
        cube10(),
        ball,
        Tolerance::default(),
        None,
    )
    .expect("common serves");
    assert!(r.faces.is_empty());
}

#[test]
fn the_verdict_needs_agreement() {
    let t = |inside, outside, on_boundary| ProbeTally {
        probes: inside + outside + on_boundary,
        inside,
        outside,
        on_boundary,
    };
    assert_eq!(t(3, 0, 2).verdict(), Some(PointClassification::Inside));
    assert_eq!(t(0, 1, 0).verdict(), Some(PointClassification::Outside));
    assert_eq!(t(0, 0, 4).verdict(), Some(PointClassification::OnBoundary));
    assert_eq!(t(2, 1, 0).verdict(), None);
    assert!(matches!(
        t(0, 0, 4).undecided(1),
        BoolError::RelationUndecided {
            operand: 1,
            on_boundary: 4,
            ..
        }
    ));
}

/// Face-interior probes lie off the face's boundary: on the unturned
/// ball none is a pole or on the seam.
#[test]
fn interior_probes_avoid_the_seam_and_poles() {
    let ball = BRep::sphere(5.0, Tolerance::default()).expect("sphere builds");
    let face = ball.faces.keys().next().expect("a face");
    let pts = face_interior_points(&ball, face, INTERIOR_PER_FACE);
    assert_eq!(pts.len(), INTERIOR_PER_FACE);
    for p in pts {
        assert!((p.x * p.x + p.y * p.y).sqrt() > 0.1, "{p:?} at a pole");
    }
}
