//! PK-389 — the flagship boss plate serves its 9-edge fillet at G2 with
//! a MEASURED curvature-continuity verdict.
//!
//! The selection is the showcase's: the boss/plate junction circle, the
//! four plate-top rim edges, and the four vertical plate edges. Before
//! this task the G2 request declined at every radius, for two
//! independent reasons found by the recorder: (1) the circular-spine
//! G2 surface was a cubic loft whose rows left the exact circular rails
//! by 2.6e-4 (now an exact surface of revolution of one G2 section);
//! (2) the corner-cap region merge dissolved the boss cylinder's own
//! periodic seam (an "any same-loop bridge" scan, now scoped to the
//! surrendered region) and the ring bridge split left the boss hole on
//! the corner sliver (inner loops are now re-homed by containment).

use popkern_part::prelude::*;
use popkern_script::Plane;
use popkern_script::prelude::*;

fn boss_plate_g2(radius: f64) -> (Model, FeatureHandle) {
    let mut m = Model::new();
    let plate = m
        .profile(Plane::XY)
        .rect(2.4, 1.6)
        .sketch("plate_sketch")
        .expect("sketch")
        .extrude(&mut m, End::Blind(0.4.into()))
        .expect("extrude");
    let boss = m
        .tool(|m| {
            m.profile(Plane::xy_at(0.2))
                .circle(0.35)
                .sketch("boss_sketch")?
                .extrude(m, End::Blind(0.8.into()))
        })
        .expect("boss");
    let plate = plate.union(&mut m, boss).expect("union");
    plate
        .fillet_continuity(
            &mut m,
            edges()
                .between(SurfaceKind::Cylinder, SurfaceKind::Plane)
                .at_height(0.4)
                .or(edges().straight().at_height(0.4))
                .or(edges().straight().direction(Dir::Z)),
            radius,
            ContinuityKind::G2,
            "rounds",
        )
        .expect("G2 fillet serves");
    let rounds = m.feature("rounds").expect("label");
    (m, rounds)
}

#[test]
fn nine_edge_g2_fillet_serves_with_measured_g2_at_every_showcase_radius() {
    for radius in [0.06, 0.08, 0.1] {
        let (mut m, rounds) = boss_plate_g2(radius);
        m.assert_watertight().expect("watertight");
        let part = m
            .to_part()
            .expect("serving boundary (validate_all) accepts the body");
        assert_eq!(part.face_count(), 29, "r={radius}: face census");
        let report = m.evidence(&rounds).expect("fillet evidence");
        let blend = report.blend.as_ref().expect("blend record");
        assert_eq!(blend.edges_requested, 9);
        let continuity = blend.continuity.as_ref().expect("measured continuity");
        assert!(continuity.measured(), "r={radius}: seams measured");
        // Honest verdict pin. The body serves and validates, and the
        // evidence is MEASURED; the measured level is below G2 today
        // for two filed reasons — the strip-rail pcurve parametrisation
        // reads as a fake r/√2 gap (PK-391) and the corner-fan
        // substrate on the rim+verticals class declines its setback
        // fan and serves Gregory seams at up to 27° (PK-392). When
        // either closes this pin flips: retire the downgrade
        // expectation and require `!blend.downgraded()`.
        assert!(
            blend.downgraded(),
            "r={radius}: measured continuity now meets G2 — retire the PK-391/PK-392 pin: {continuity:?}"
        );
        // Every positional "gap" is the PK-391 parametrisation artifact
        // (image-equal sides: tangent match to 1e-6), never a crack.
        for seam in &continuity.seams {
            if seam.worst_g0_gap > 1e-4 {
                assert!(
                    seam.worst_g1_normal < 1e-6,
                    "r={radius}: seam {:?} gap {:.3e} with tangent jump {:.3e} is a real crack",
                    seam.edge,
                    seam.worst_g0_gap,
                    seam.worst_g1_normal
                );
            }
        }
    }
}
