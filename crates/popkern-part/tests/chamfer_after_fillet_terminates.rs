//! PK-390 — a chamfer of the boss-top rim on the G1-filleted flagship
//! boss plate must serve (or decline) promptly.
//!
//! Before the fix the chamfer's edge-reconcile pass intersected each
//! straight-edge fillet cylinder with its corner sphere — the same
//! rolling-ball radius, centre on the strip axis, so the pair is
//! tangent along a full circle. The marching solver seeded ~2100
//! tangent-classified points and marched every one to its 20 001-step
//! cap (26 ms each) without a crossing certificate: ~55 s per pair
//! call, three calls per pass, and the test harness reported a hang
//! (>7 min at 97 % CPU). The coaxial cylinder×sphere analytic arm now
//! answers that pair exactly (tangent → silent), and the marching
//! driver marks a certificate-failed demoted trace consumed so no
//! other tangent band can repeat the identical march.
use popkern_part::prelude::*;
use std::time::Instant;

#[test]
fn boss_top_chamfer_after_g1_fillet_serves_within_the_wall_clock_gate() {
    let plate = Profile::on(Plane::xy_at(0.0))
        .rect(2.4, 1.6)
        .extrude(0.4)
        .expect("plate");
    let boss = Profile::on(Plane::xy_at(0.2))
        .circle(0.35)
        .extrude(0.8)
        .expect("boss");
    let body = plate.union(&boss).expect("union");
    let rounded = body
        .fillet(0.08, |e| {
            e.between(SurfaceKind::Cylinder, SurfaceKind::Plane)
                .at_height(0.4)
                .or(|q| q.straight().at_height(0.4))
                .or(|q| q.straight().direction(Dir::Z))
        })
        .expect("G1 fillet serves");
    let faces_before = rounded.face_count();
    let t = Instant::now();
    let chamfered = rounded
        .chamfer(0.05, |e| e.circular().at_height(1.0))
        .expect("boss-top chamfer serves");
    let elapsed = t.elapsed();
    assert!(
        elapsed.as_secs_f64() < 20.0,
        "chamfer took {elapsed:?} (gate 20 s; pre-fix >7 min)"
    );
    // One bevel face added between the boss top and the boss wall.
    assert_eq!(chamfered.face_count(), faces_before + 1);
    let issues = chamfered.brep().validate_all();
    assert!(issues.is_empty(), "validate_all: {issues:?}");
    assert_eq!(chamfered.open_edge_count().expect("mesh"), 0);
}
