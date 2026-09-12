//! Unit gates for the tilt-aware entry ring (PK-480): probe
//! selection, plunge/closure math, and the steep-entry decline.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::ProfileShape;

use super::{EntryProbe, entry_ring_profile, probe_entry};
use crate::FastenerError;
use crate::engine::extrude_shape;

/// A `w × d × h` block with its top face on z = 0 (material below,
/// so a hole axis of −Z̄... here +Z̄ points into material means the
/// canonical hole frame drills downward through the top at z = h).
fn block(w: f64, d: f64, h: f64) -> BRep {
    extrude_shape(
        ProfileShape::Rect {
            width: w,
            height: d,
        },
        0.0,
        h,
        "test block",
        Tolerance::default(),
    )
    .expect("test block")
}

#[test]
fn perpendicular_entry_probes_flat_and_hull_free() {
    let body = block(20.0, 20.0, 8.0);
    let probe = probe_entry(
        &body,
        Point3::new(1.0, -0.5, 8.0),
        Vec3::new(0.0, 0.0, -1.0),
        4.0,
    );
    assert!(probe.tilt_tan.abs() < 1e-12, "tilt {}", probe.tilt_tan);
    assert!(
        probe.hull_height.abs() < 1e-12,
        "hull {}",
        probe.hull_height
    );
}

#[test]
fn tilted_entry_probes_the_face_tilt_and_hull_rise() {
    let angle = 22.0_f64.to_radians();
    let body = block(20.0, 20.0, 8.0)
        .transform(&Transform3::rotation_y(angle))
        .expect("tilt block");
    // Drill straight down (−Z): the rotated top face is `angle` off
    // perpendicular. Origin on the rotated top plane at x-ish 0.
    let origin = Point3::new(0.0, 0.0, 8.0 / angle.cos());
    let probe = probe_entry(&body, origin, Vec3::new(0.0, 0.0, -1.0), 4.0);
    assert!(
        (probe.tilt_tan - angle.tan()).abs() < 1e-9,
        "tilt {} want {}",
        probe.tilt_tan,
        angle.tan()
    );
    // The rotated block's highest vertex rises above the origin.
    assert!(
        probe.hull_height > 1.0,
        "hull {} should see the uphill corners",
        probe.hull_height
    );
}

#[test]
fn side_wall_near_the_mouth_is_not_the_entry_face() {
    let body = block(20.0, 20.0, 8.0);
    // Hole 1.5 from the x = 10 wall: the wall plane passes well
    // inside the pierce band as measured by DISTANCE, but the axis
    // ray never pierces it — the probe must still read the flat top.
    let probe = probe_entry(
        &body,
        Point3::new(8.5, 0.0, 8.0),
        Vec3::new(0.0, 0.0, -1.0),
        4.0,
    );
    assert!(probe.tilt_tan.abs() < 1e-12, "tilt {}", probe.tilt_tan);
}

#[test]
fn missing_entry_face_falls_back_to_perpendicular() {
    let body = block(20.0, 20.0, 8.0);
    // Origin far above the body: no plane pierces within the band.
    let probe = probe_entry(
        &body,
        Point3::new(0.0, 0.0, 30.0),
        Vec3::new(0.0, 0.0, -1.0),
        4.0,
    );
    assert_eq!(probe.tilt_tan, 0.0);
}

#[test]
fn perpendicular_profile_matches_the_design_ring() {
    let probe = EntryProbe::perpendicular();
    let poly = entry_ring_profile(2.5, 3.3, 0.25, &probe, "test").expect("profile");
    // Zero plunge: cone still crosses z = 0 at the design rim; lid a
    // whisker above the mouth; wall one eps outside the rim.
    let want = [(2.5, -0.5), (3.55, -0.5), (3.55, -0.25), (2.5, 0.8)];
    for (got, want) in poly.iter().zip(want) {
        assert!((got.0 - want.0).abs() < 1e-12 && (got.1 - want.1).abs() < 1e-12);
    }
}

#[test]
fn tilted_profile_plunges_by_rim_times_tan_and_lids_above_hull() {
    let probe = EntryProbe {
        tilt_tan: 0.4,
        hull_height: 5.0,
    };
    let poly = entry_ring_profile(2.5, 3.0, 0.25, &probe, "test").expect("profile");
    let plunge = 3.0 * 0.4;
    let apex_sum = 3.0 + plunge;
    assert_eq!(
        poly,
        vec![
            (2.5, -5.5),
            (apex_sum + 0.25, -5.5),
            (apex_sum + 0.25, -0.25),
            (2.5, apex_sum - 2.5),
        ]
    );
    // The cone keeps its 45° slope through the extension.
    let (r0, z0) = poly[2];
    let (r1, z1) = poly[3];
    assert!(((r0 + z0) - (r1 + z1)).abs() < 1e-12, "cone not 45°");
}

#[test]
fn steep_entry_declines_with_the_measured_tilt() {
    let probe = EntryProbe {
        tilt_tan: 1.2,
        hull_height: 0.0,
    };
    let err = entry_ring_profile(2.5, 3.0, 0.25, &probe, "tapped hole")
        .expect_err("must decline past the tilt bound");
    match err {
        FastenerError::EntryTooSteep {
            stage,
            tilt_deg,
            limit_deg,
        } => {
            assert_eq!(stage, "tapped hole");
            assert!((tilt_deg - 1.2_f64.atan().to_degrees()).abs() < 1e-9);
            assert!(limit_deg < tilt_deg);
        }
        other => panic!("wrong decline: {other:?}"),
    }
}
