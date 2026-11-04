//! A record piece outside the unused face's trim is a phantom split,
//! not a leak: the keyway slot cut into a bored block serves through
//! the complex.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::BoolOperation;
use crate::complex::casting_scenes::Casting;

/// The axis-aligned box over `[x0, x1] × [y0, y1] × [z0, z1]`.
fn cuboid(x: (f64, f64), y: (f64, f64), z: (f64, f64), tol: Tolerance) -> BRep {
    let profile = [
        Point3::new(x.0, y.0, z.0),
        Point3::new(x.1, y.0, z.0),
        Point3::new(x.1, y.1, z.0),
        Point3::new(x.0, y.1, z.0),
    ];
    BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), z.1 - z.0, tol).expect("cuboid")
}

/// A 1 × 1 × 0.3 block centred on the z axis, bored through by a
/// 0.15-radius hole.
fn bored_block(tol: Tolerance) -> BRep {
    let block = cuboid((-0.5, 0.5), (-0.5, 0.5), (0.0, 0.3), tol);
    let bore = BRep::cylinder(0.15, 0.5, tol)
        .expect("bore")
        .transform(&Transform3::translation(0.0, 0.0, -0.1))
        .expect("placed bore");
    crate::boolean(BoolOperation::Subtraction, block, bore, tol, None).expect("bored block")
}

/// A 0.08 × 0.06 slot over y ∈ [0.14, 0.20], through the block: its
/// inner wall y = 0.14 crosses the bore, so the records where that wall
/// meets the block's top and bottom faces run, between x = ±0.04,
/// across those faces' holes.
fn key(tol: Tolerance) -> BRep {
    cuboid((-0.04, 0.04), (0.14, 0.2), (-0.14, 0.44), tol)
}

#[test]
fn a_record_across_a_hole_does_not_decline_the_complex() {
    let tol = Tolerance::default();
    let block = bored_block(tol);
    let key = key(tol);
    let ix = crate::intersect::find_intersections(&block, &key, &tol).expect("intersections");
    let served = crate::complex::serve::rescue(&block, &key, &ix, BoolOperation::Subtraction, &tol)
        .expect("manifold");
    let (body, _) = served.expect("the complex serves the keyway cut");
    body.validate().expect("valid");
    // Material where the block stays, none in the bore or the slot.
    let boxes = crate::intersect::compute_face_aabbs(&body);
    for (p, inside) in [
        ((0.3, 0.0, 0.15), true),
        ((0.0, -0.3, 0.15), true),
        ((0.1, 0.17, 0.15), true),
        ((0.0, 0.0, 0.15), false),
        ((0.0, 0.17, 0.15), false),
        ((0.0, 0.17, 0.02), false),
        ((0.0, 0.17, -0.05), false),
    ] {
        let point = Point3::new(p.0, p.1, p.2);
        let got =
            crate::classify::point_in_solid_final(&point, &body, &boxes, &tol).expect("classifies");
        let want = if inside {
            crate::PointClassification::Inside
        } else {
            crate::PointClassification::Outside
        };
        assert_eq!(got, want, "at {p:?}");
    }
}

/// The tangent boss casting ([`Casting::tangent_boss`]): the boss's top
/// generator is tangent to the web's top face y = -17 and touches the
/// ring's extreme there. The boss's records with the ring wall run
/// inside the web, through the web's notch in the ring wall, and end on
/// the notch rim at that tangent point. An end on the rim is where any
/// record leaving the face lies, no evidence the record splits it
/// (measured: the six-axis arm's end casting declined in every step of
/// its tangency sweep). The boss stops short of the web's end face,
/// whose coplanar cap is a separate class.
#[test]
fn a_record_through_a_notch_ending_on_its_rim_does_not_decline_the_complex() {
    Casting::tangent_boss(-50.0, 29.0).complex_union_serves();
}

/// The boss swept through its tangency with the web's top face (the
/// arm casting's own sweep): within the tolerance band its record with
/// the ring wall, running through the notch, touches the notch rim at
/// the tangent point between its ends. A touch is not a ride along the
/// rim, and the record still misses the ring wall.
#[test]
fn a_record_touching_the_notch_rim_between_its_ends_does_not_decline_the_complex() {
    for delta in [
        0.0, -1e-6, -1e-7, 1e-7, 5e-7, 1e-6, 1e-4, 1e-3, 3e-3, 1e-2, -1e-3, -0.1, 0.1,
    ] {
        Casting::tangent_boss(-50.0 + delta, 29.0).complex_union_serves();
    }
}

/// The boss reaching the web's end face x = -740: its far cap is
/// coplanar with the web's end there, and the boss wall's record with
/// that face is the cap's rim, riding the wall's boundary all round.
/// With the wall's chart seam at the top the record ran two periods from
/// the trim's window and the riding test, trying one period each way,
/// kept half of it (PK-1018): the web's end was never cut along the
/// rest, the two copies of the coincident patch differed and stayed
/// separate complex faces (an ambiguous radial order on the web's edge
/// inside the cap; with the boss 0.1 low, a half-registered record).
/// Every seam placement serves.
#[test]
fn a_boss_cap_flush_with_the_web_end_serves() {
    use std::f64::consts::{FRAC_PI_2, PI};
    for boss_seam in [-FRAC_PI_2, 0.0, FRAC_PI_2, PI, 2.0] {
        for boss_y in [-50.0, -50.1, -49.9] {
            Casting {
                boss_seam,
                ..Casting::tangent_boss(boss_y, 30.0)
            }
            .complex_union_serves();
        }
    }
}
