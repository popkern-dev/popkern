//! Two spherical sectors of one ball, about axes 60 degrees apart,
//! sharing the ball's centre as their cones' apex (PK-1138).
//!
//! The cones meet in two rulings through the apex. Their intersection
//! arm clipped each ruling to both cones' parameter windows, padded
//! past the apex onto the other nappe, so each ruling ran from the rim
//! through the apex vertex to the padding; the face has no edge at its
//! apex for the clip to cut a record at, the vertex table refused the
//! far ends (5 off the face) and the fallback's graph walked the
//! sector's cone wall twice. Built from the Rust constructors (each
//! sector a ball less a cone), five of seven self-turns declined the
//! common on a loop gap. A cone face lies on one nappe, and the
//! window is cut back to it.

use std::f64::consts::{FRAC_1_SQRT_2, PI};

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// The sector of the ball of radius 10 inside the cone of half-angle
/// 45 degrees about z from its centre, turned `spin` about z.
fn sector(spin: f64) -> BRep {
    let tol = Tolerance::default();
    let ball = BRep::sphere(10.0, tol).expect("ball");
    let cone = BRep::cone(10.0, 10.0, tol)
        .and_then(|c| c.transform(&Transform3::rotation_x(PI)))
        .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, 10.0)))
        .expect("cone");
    boolean(BoolOperation::Intersection, ball, cone, tol, None)
        .and_then(|s| Ok(s.transform(&Transform3::rotation_z(spin))?))
        .expect("sector")
}

/// Each operation of the sector turned `spin_a` and the one turned
/// `spin_b` then tilted 60 degrees about y, at `scale`: its volume.
fn assert_sectors(spin_a: f64, spin_b: f64, scale: f64) {
    let base = Tolerance::default();
    let tol = base.with_linear(base.linear() * scale);
    let tilt = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        PI / 3.0,
    )
    .expect("tilt");
    let place = Transform3::scale_uniform(scale);
    // R^3/3 times the solid angles: a 45 degree cap's, and the lens of
    // two of them 60 degrees apart (the bench row's quadrature).
    let sector_v = 2000.0 * PI / 3.0 * (1.0 - FRAC_1_SQRT_2);
    let common = 146.516_406_055_211_18;
    for (op, want) in [
        (BoolOperation::Intersection, common),
        (BoolOperation::Union, 2.0 * sector_v - common),
        (BoolOperation::Subtraction, sector_v - common),
    ] {
        let a = sector(spin_a).transform(&place).expect("place");
        let b = sector(spin_b)
            .transform(&tilt)
            .and_then(|b| b.transform(&place))
            .expect("place");
        let body = boolean(op, a, b, tol, None)
            .unwrap_or_else(|e| panic!("{op:?} at {spin_a}, {spin_b}, x{scale}: {e:?}"));
        assert!(
            body.validate_all().is_empty(),
            "{op:?}: {:?}",
            body.validate_all()
        );
        let v = popkern_measure::body_volume(&body, &tol).expect("volume");
        let want = want * scale.powi(3);
        assert!(
            (v - want).abs() < 1e-6 * want,
            "{op:?} at {spin_a}, {spin_b}, x{scale}: volume {v} for {want}"
        );
    }
}

#[test]
fn sectors_sharing_an_apex_cut_unite_and_meet() {
    for (a, b) in [(0.0, 0.0), (0.0, 0.7), (0.7, 0.7), (PI, 0.0), (0.3, 1.9)] {
        assert_sectors(a, b, 1.0);
    }
}

#[test]
fn kilometre_sectors_sharing_an_apex() {
    assert_sectors(0.0, 0.0, 1e3);
    assert_sectors(PI, 0.0, 1e3);
}
