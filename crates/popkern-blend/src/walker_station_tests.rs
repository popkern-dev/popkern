use std::f64::consts::{PI, TAU};

use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;

use crate::contact::{BallSide, ContactSolution};
use crate::walker::{WalkOptions, WalkRequest, WalkerError, walk_spine};

const R_BOSS: f64 = 2.0;
const R_BALL: f64 = 0.3;

/// A boss on a plate: cylinder radius 2 along z, plane z = 0, the
/// intersection circle as a closed spine (the arc's own angle basis
/// decides where angle 0 sits; every expectation is evaluated).
fn boss_rim() -> (SurfaceGeometry, SurfaceGeometry, CurveGeometry) {
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: R_BOSS,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let spine = CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: R_BOSS,
        start_angle: 0.0,
        end_angle: TAU,
    };
    (cyl, plane, spine)
}

fn request<'a>(
    cyl: &'a SurfaceGeometry,
    plane: &'a SurfaceGeometry,
    spine: &'a CurveGeometry,
    radius_at: &'a dyn Fn(f64) -> f64,
    seed: Option<((f64, f64), (f64, f64))>,
) -> WalkRequest<'a> {
    WalkRequest {
        surface_0: cyl,
        surface_1: plane,
        side_0: BallSide::Outside,
        side_1: BallSide::Outside,
        radius_at,
        spine,
        t_range: (0.0, TAU),
        uv0_range: ((0.0, TAU), (-10.0, 10.0)),
        uv1_range: ((-10.0, 10.0), (-10.0, 10.0)),
        initial_seed: seed,
        surface_0_periodicity: (Some(TAU), None),
        surface_1_periodicity: (None, None),
    }
}

/// The convex rolling ball at spine point `p`: centre radially
/// outside by the ball radius, lifted by the ball radius.
fn ball_at(p: Point3) -> ContactSolution {
    let radial = Vec3::new(p.x, p.y, 0.0) * (1.0 / R_BOSS);
    let center = Point3::new(0.0, 0.0, R_BALL) + radial * (R_BOSS + R_BALL);
    ContactSolution {
        uv0: (0.0, 0.0),
        uv1: (0.0, 0.0),
        center,
        contact0: Point3::new(0.0, 0.0, R_BALL) + radial * R_BOSS,
        contact1: Point3::new(0.0, 0.0, 0.0) + radial * (R_BOSS + R_BALL),
        iterations: 0,
        step_rejections: 0,
        residual: 0.0,
        condition: 1.0,
        jacobian: [[0.0; 4]; 4],
    }
}

/// The antipode of the spine start lies in the same spine plane, so
/// a seed there converges to the antipodal ball; the walker must
/// re-seed from the start and deliver the ball AT the start.
#[test]
fn a_seed_at_the_antipode_of_a_closed_rim_is_walked_back_to_its_station() {
    let (cyl, plane, spine) = boss_rim();
    let radius_at = |_: f64| R_BALL;
    let antipode = spine.point_at(PI);
    let seed = (
        cyl.closest_point(&antipode).expect("on the cylinder"),
        plane.closest_point(&antipode).expect("on the plane"),
    );
    let seeded = request(&cyl, &plane, &spine, &radius_at, Some(seed));
    let walk = walk_spine(&seeded, &WalkOptions::default(), &Tolerance::default())
        .expect("re-seeded walk serves");
    let first = walk.samples.first().expect("samples").solution.center;
    let expected = ball_at(spine.point_at(0.0)).center;
    assert!(
        first.distance_to(&expected) < 1e-6,
        "first ball centre {first:?} is at the spine start's station {expected:?}"
    );
    // The unseeded walk agrees to roundoff.
    let unseeded = request(&cyl, &plane, &spine, &radius_at, None);
    let reference =
        walk_spine(&unseeded, &WalkOptions::default(), &Tolerance::default()).expect("walk serves");
    let reference = reference.samples.first().expect("samples").solution.center;
    assert!(first.distance_to(&reference) < 1e-9);
}

/// The certificate itself: the antipodal ball is nearer the station
/// half a turn away than to its own, by far more than the margin.
#[test]
fn the_certificate_names_the_station_the_ball_belongs_to() {
    let (cyl, plane, spine) = boss_rim();
    let radius_at = |_: f64| R_BALL;
    let req = request(&cyl, &plane, &spine, &radius_at, None);
    let antipodal = ball_at(spine.point_at(PI));
    match super::certify_station(&req, 0.0, &antipodal) {
        Err(WalkerError::StationMismatch {
            t_spine,
            nearest_t,
            center,
        }) => {
            assert_eq!(t_spine, 0.0);
            assert!((nearest_t - PI).abs() < 1e-9, "nearest station {nearest_t}");
            assert_eq!(center, antipodal.center);
        }
        other => panic!("expected a station mismatch, got {other:?}"),
    }
    let own = ball_at(spine.point_at(0.0));
    assert!(super::certify_station(&req, 0.0, &own).is_ok());
    // The same ball asked for at the end of the closed spine (2π ≡ 0)
    // is at its station too: the neighbourhood wraps.
    assert!(super::certify_station(&req, TAU, &own).is_ok());
    // A quarter turn is unambiguous either way.
    let quarter = ball_at(spine.point_at(PI / 2.0));
    assert!(super::certify_station(&req, PI / 2.0, &quarter).is_ok());
    assert!(matches!(
        super::certify_station(&req, 0.0, &quarter),
        Err(WalkerError::StationMismatch { .. })
    ));
}
