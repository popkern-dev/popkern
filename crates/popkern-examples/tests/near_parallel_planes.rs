//! Planes nearly parallel but apart past the tolerance over their faces
//! (PK-1137).
//!
//! Plane x plane decided parallelism by the angle between the normals
//! against the angular tolerance (1e-4 rad) and then coincidence by one
//! origin's distance: planes a ten-thousandth of a radian apart through
//! a common line, 5e-4 apart across a face 10 wide, read coincident. A
//! wedge cut a ten-thousandth steep lost its cut along the top face and
//! declined; a millionth steep served an EMPTY body; a cube united with
//! itself turned 0.001 degrees declined. Coincidence is now distance
//! over the faces' boxes. The wedge's end faces, 5e-5 thick at a
//! hundred-thousandth, then went to sliver removal, which absorbed them
//! into a neighbour 50 tolerances away: absorption moves no geometry,
//! so it now takes only slivers within the tolerance of their absorber.

use std::f64::consts::PI;

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn slab(lo: (f64, f64, f64), hi: (f64, f64, f64), tol: Tolerance) -> BRep {
    let rect = |z: f64| {
        vec![
            Point3::new(lo.0, lo.1, z),
            Point3::new(hi.0, lo.1, z),
            Point3::new(hi.0, hi.1, z),
            Point3::new(lo.0, hi.1, z),
        ]
    };
    BRep::polygon_frustum(&rect(lo.2), &rect(hi.2), tol).expect("slab")
}

/// The cube `[-5, 5]² × [-10, 0]` cut by the half-space above
/// `z = slope·y`.
fn wedge(slope: f64) -> Result<BRep, popkern_bool::BoolError> {
    let tol = Tolerance::default();
    let tilt = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        slope.atan(),
    )
    .expect("tilt");
    boolean(
        BoolOperation::Intersection,
        slab((-5.0, -5.0, -10.0), (5.0, 5.0, 0.0), tol),
        slab((-20.0, -20.0, 0.0), (20.0, 20.0, 20.0), tol)
            .transform(&tilt)
            .expect("tilted"),
        tol,
        None,
    )
}

#[test]
fn a_shallow_wedge_serves_its_volume() {
    let tol = Tolerance::default();
    for slope in [1e-3, 1e-4, 1e-5] {
        let body = wedge(slope).unwrap_or_else(|e| panic!("{slope}: {e:?}"));
        assert!(
            body.validate_all().is_empty(),
            "{slope}: {:?}",
            body.validate_all()
        );
        assert_eq!(body.faces.len(), 5, "{slope}");
        let v = popkern_measure::body_volume(&body, &tol).expect("volume");
        let want = 125.0 * slope;
        assert!(
            (v - want).abs() < 1e-9 * want,
            "{slope}: volume {v} for {want}"
        );
    }
}

/// A millionth steep the wedge is at most five tolerances thick (PK-1131
/// owns it); it must decline, never serve an empty body.
#[test]
fn a_wedge_a_millionth_steep_never_serves_empty() {
    if let Ok(body) = wedge(1e-6) {
        let v = popkern_measure::body_volume(&body, &Tolerance::default()).expect("volume");
        assert!(
            body.faces.len() >= 5 && (v - 1.25e-4).abs() < 1e-9,
            "served {} faces, volume {v}",
            body.faces.len()
        );
    }
}

#[test]
fn a_cube_united_with_itself_turned_a_hair() {
    let tol = Tolerance::default();
    let cube = || slab((-5.0, -5.0, 0.0), (5.0, 5.0, 10.0), tol);
    let turn = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        0.001 * PI / 180.0,
    )
    .expect("turn");
    let body = boolean(
        BoolOperation::Union,
        cube(),
        cube().transform(&turn).expect("turned"),
        tol,
        None,
    )
    .expect("serves");
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    let v = popkern_measure::body_volume(&body, &tol).expect("volume");
    let want = 1_000.008_726_493_952_7;
    assert!((v - want).abs() < 1e-9 * want, "volume {v} for {want}");
}

/// A box `[0, 10]² × [0, 1]` united with a block beside it whose top
/// rises from the box's top edge by `rise` over its 10, the scene
/// turned by `place`.
fn creased_union(rise: f64, place: &Transform3) -> Result<BRep, popkern_bool::BoolError> {
    let tol = Tolerance::default();
    let quad = |x0: f64, x1: f64, z0: f64, z1: f64| {
        vec![
            Point3::new(x0, 0.0, z0),
            Point3::new(x1, 0.0, z1),
            Point3::new(x1, 10.0, z1),
            Point3::new(x0, 10.0, z0),
        ]
    };
    let solid = |bottom: Vec<Point3>, top: Vec<Point3>| {
        BRep::polygon_frustum(&bottom, &top, tol)
            .and_then(|b| b.transform(place))
            .expect("solid")
    };
    let a = solid(quad(0.0, 10.0, 0.0, 0.0), quad(0.0, 10.0, 1.0, 1.0));
    let b = solid(
        quad(10.0, 20.0, 0.0, 0.0),
        quad(10.0, 20.0, 1.0, 1.0 + rise),
    );
    boolean(BoolOperation::Union, a, b, tol, None)
}

/// Two tops meeting at a crease too shallow for the surface comparison
/// (normals within the angular tolerance, one plane's origin on the
/// other) stay two faces; they merged into one, its far edge `rise` off
/// its plane, validate_all failing on a pcurve 10 off (PK-1148). A rise
/// under the tolerance still merges.
#[test]
fn a_shallow_crease_keeps_both_faces() {
    let tol = Tolerance::default();
    let turned =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.7, 0.3), 0.6)
            .expect("axis");
    for place in [Transform3::identity(), turned] {
        for (rise, faces) in [(1e-3, 7), (0.01, 7), (0.1, 7), (1e-8, 6)] {
            let body = creased_union(rise, &place).expect("union");
            assert!(
                body.validate_all().is_empty(),
                "rise {rise}: {:?}",
                body.validate_all()
            );
            assert_eq!(body.faces.len(), faces, "rise {rise}");
            let v = popkern_measure::body_volume(&body, &tol).expect("volume");
            let want = 200.0 + 50.0 * rise;
            // Exact but for a merged top: within the tolerance times
            // its area.
            let bound = 1e-9 * want + tol.linear() * 200.0;
            assert!((v - want).abs() < bound, "rise {rise}: volume {v}");
        }
    }
}
