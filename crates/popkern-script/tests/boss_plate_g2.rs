//! PK-389/PK-398 — the flagship boss plate serves its 9-edge fillet at
//! G2 with a MEASURED curvature-continuity verdict that meets the request.
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
        // Corner method decides the census: the KP20 setback fans
        // serve at EVERY showcase radius. r=0.1 joined when the
        // crossing detector became transversality-structural (PK-392:
        // the setback arc is UV-coincident with the pre-carve apex
        // boundary by construction; the old count was a
        // projection-noise knife-edge).
        let expected_faces = 65;
        assert_eq!(part.face_count(), expected_faces, "r={radius}: face census");
        let report = m.evidence(&rounds).expect("fillet evidence");
        let blend = report.blend.as_ref().expect("blend record");
        assert_eq!(blend.edges_requested, 9);
        let continuity = blend.continuity.as_ref().expect("measured continuity");
        assert!(continuity.measured(), "r={radius}: seams measured");
        // Honest verdict pin: the body serves, validates, and the
        // MEASURED continuity meets the request. Until PK-398 the fans
        // served as a recorded G1 downgrade (cap↔strip seam class
        // 6.85e-2, then a transverse-curvature class of 1.4): the
        // border rows were a cubic interpolant marched from control
        // points; the exact section rows with corner-consistent
        // Hermite blocks retired both.
        assert!(
            !blend.downgraded(),
            "r={radius}: measured continuity below the G2 request: {continuity:?}"
        );
        assert_eq!(
            continuity.achieved,
            Continuity::G2,
            "r={radius}: achieved level"
        );
        // No seam gap within orders of magnitude of the retired rail
        // artifact (0.64–0.73·r) or the cubic-sagitta class (4.2e-4).
        assert!(
            continuity.worst_g0_gap < 1e-8,
            "r={radius}: worst seam gap {:.3e}",
            continuity.worst_g0_gap
        );
    }
}
