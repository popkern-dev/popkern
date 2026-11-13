//! Two balls cutting a box's top in circles a few tolerances apart
//! (PK-1180). Box [0, 10]^3, balls of radius 2 centred 1 over the top at
//! x = 5 -+ d/2 (y = 5), d = 2 sqrt 3 + gap: each cuts the top in a circle
//! of radius sqrt 3, the two circles `gap` apart at P = (5, 5, 10), and
//! the balls overlap in a circle (plane x = 5) whose lowest point lies
//! about 0.87 gap over the top. Unions declined from one tolerance to
//! twice it, as built (the first ball's seam meridian through P) and
//! with the seams turned away: the clip read the overlap circle's near
//! miss of the first ball's rim latitude (7.7e-7 rad, 1.5e-6 in space)
//! as two crossings and kept the 3e-4 between them; a result's
//! entities carry the model's tolerance, read as their own, so in the
//! second union the seam's crossing of the second ball joined the
//! seam's vertex 1.5e-6 down it; a vertex 1.5e-6 off a record split it,
//! and two records 1.5e-6 apart met.

use popkern_bool::{
    BoolError, BoolOperation, BoolOptions, NonManifoldPolicy, SetExpr, boolean_with_options,
    nary_boolean,
};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// The box and the two balls for `gap`, each ball turned `own` about its
/// vertical axis (moving its seam off P), the scene turned `turn` about
/// an axis through P.
fn scene(gap: f64, own: f64, turn: f64) -> [BRep; 3] {
    let tol = Tolerance::default();
    let spin = Transform3::rotation_around_axis(
        Point3::new(5.0, 5.0, 10.0),
        Vec3::new(0.3, 1.0, 0.2),
        turn,
    )
    .expect("axis");
    let cube = BRep::prism(
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 10.0, 0.0),
            Point3::new(0.0, 10.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        10.0,
        tol,
    )
    .and_then(|b| b.transform(&spin))
    .expect("box");
    let d = 2.0 * 3f64.sqrt() + gap;
    let ball = |x: f64| {
        BRep::sphere(2.0, tol)
            .and_then(|b| b.transform(&Transform3::rotation_z(own)))
            .and_then(|b| b.transform(&Transform3::translation(x, 5.0, 11.0)))
            .and_then(|b| b.transform(&spin))
            .expect("ball")
    };
    [cube, ball(5.0 - d / 2.0), ball(5.0 + d / 2.0)]
}

fn volume(body: &BRep) -> f64 {
    popkern_measure::body_volume(body, &Tolerance::default()).expect("volume")
}

/// The box united with one ball and then the other, and the three
/// united at once, under `policy`.
fn both_paths([cube, a, b]: &[BRep; 3], policy: NonManifoldPolicy) -> [Result<BRep, BoolError>; 2] {
    let tol = Tolerance::default();
    let options = BoolOptions {
        non_manifold: policy,
    };
    let union = |x: BRep, y: &BRep| {
        boolean_with_options(BoolOperation::Union, x, y.clone(), tol, None, &options)
    };
    [
        union(cube.clone(), a).and_then(|x| union(x, b)),
        nary_boolean(&[cube, a, b], &SetExpr::Any, tol, &options),
    ]
}

/// Within the tolerance the rims are one place at P, where the space
/// outside the three pinches (its channels either side of P meet only
/// there): a point contact, as a ball resting on a face within the
/// tolerance is. Both paths decline by default with the contact at P,
/// and under the general policy serve one valid general body holding
/// it, of one volume, continuous with the gap just past the tolerance.
/// The N-ary path placed each of the three records' pairwise touches at
/// the middle of its stretch within the tolerance (4e-13, 4.2e-5 and
/// 1.9e-4 from P) and served a solid whose rims shared a span; the
/// two-body path, with the seams turned, kept a point of the overlap
/// circle against the first ball's cap (the stretch's ends read as an
/// exit and a re-entry) and served both caps whole, the lens counted
/// twice (PK-1181). Exactly at the tolerance the gap measures within
/// the reach by less than the records' fit, and is not tested here
/// (PK-1183).
#[test]
fn within_the_tolerance_the_unions_touch_at_one_point() {
    let p = Point3::new(5.0, 5.0, 10.0);
    for (own, turn) in [(0.0, 0.0), (0.9, 0.0), (0.0, 0.4)] {
        let past = both_paths(&scene(1.2e-6, own, turn), NonManifoldPolicy::Decline);
        let past = volume(past[0].as_ref().expect("1.2e-6 serves"));
        for gap in [0.0, 5e-7, 9e-7] {
            let what = format!("gap {gap:e} own {own} turn {turn}");
            let scene = scene(gap, own, turn);
            for (path, r) in ["two-body", "N-ary"]
                .iter()
                .zip(both_paths(&scene, NonManifoldPolicy::Decline))
            {
                match r {
                    Err(BoolError::NonManifoldResult { contacts }) => {
                        assert!(
                            contacts.len() == 1 && contacts[0].location.distance_to(&p) < 1e-5,
                            "{what}: {path} contacts {contacts:?}"
                        );
                    }
                    other => panic!("{what}: {path} by default: {:?}", other.map(|_| "served")),
                }
            }
            let general = both_paths(&scene, NonManifoldPolicy::General);
            let mut volumes = Vec::new();
            for (path, r) in ["two-body", "N-ary"].iter().zip(general) {
                let body = r.unwrap_or_else(|e| panic!("{what}: {path} general declined: {e}"));
                let issues = body.validate_all();
                assert!(issues.is_empty(), "{what}: {path} invalid: {issues:?}");
                let contacts = body.contacts();
                assert!(
                    contacts.len() == 1 && contacts[0].location.distance_to(&p) < 1e-5,
                    "{what}: {path} contacts {contacts:?}"
                );
                volumes.push(volume(&body));
            }
            assert!(
                (volumes[0] - volumes[1]).abs() < 1e-6,
                "{what}: volumes {volumes:?}"
            );
            assert!(
                (volumes[0] - past).abs() < 1e-5,
                "{what}: {} against {past}",
                volumes[0]
            );
        }
    }
}

/// From 1.2 to 100 tolerances apart, at each placement, the box united
/// with one ball and then the other, and the three united at once, serve
/// one valid solid; the two agree in volume and the volume runs linear
/// in the gap (each ball's cap over the top grows as the overlap lens
/// shrinks).
#[test]
fn the_unions_serve_from_just_past_the_tolerance() {
    let tol = Tolerance::default();
    let gaps = [1.2e-6, 1.5e-6, 2e-6, 3e-6, 5e-6, 1e-5, 1e-4];
    for (own, turn) in [(0.0, 0.0), (0.9, 0.0), (0.0, 0.4)] {
        let mut volumes = Vec::new();
        for gap in gaps {
            let [cube, a, b] = scene(gap, own, turn);
            let what = format!("gap {gap:e} own {own} turn {turn}");
            let two =
                popkern_bool::boolean(BoolOperation::Union, cube.clone(), a.clone(), tol, None)
                    .and_then(|x| {
                        popkern_bool::boolean(BoolOperation::Union, x, b.clone(), tol, None)
                    })
                    .unwrap_or_else(|e| panic!("{what}: two-body union declined: {e}"));
            let all = nary_boolean(
                &[&cube, &a, &b],
                &SetExpr::Any,
                tol,
                &BoolOptions::default(),
            )
            .unwrap_or_else(|e| panic!("{what}: N-ary union declined: {e}"));
            for (path, body) in [("two-body", &two), ("N-ary", &all)] {
                let issues = body.validate_all();
                assert!(issues.is_empty(), "{what}: {path} invalid: {issues:?}");
            }
            let (v2, vn) = (volume(&two), volume(&all));
            assert!((v2 - vn).abs() < 1e-6, "{what}: volumes {v2} and {vn}");
            volumes.push(v2);
        }
        // Linear through the widest two gaps.
        let n = gaps.len();
        let slope = (volumes[n - 1] - volumes[n - 2]) / (gaps[n - 1] - gaps[n - 2]);
        for (gap, v) in gaps.iter().zip(&volumes) {
            let line = volumes[n - 1] + slope * (gap - gaps[n - 1]);
            assert!(
                (v - line).abs() < 2e-6,
                "own {own} turn {turn} gap {gap:e}: volume {v} off the line {line}"
            );
        }
    }
}
