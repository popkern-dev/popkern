//! The six-axis arm's gripper (design §4.2, PK-724): two instances of
//! one finger part slide on the gripper body, their jaw faces ±grip/2
//! from the tool axis; closed, the jaws touch, and opening over the
//! whole travel meets nothing.

mod arm_support;

use arm_support::*;
use popkern_assembly::Assembly;

/// Grip at the script's default (the fingers are authored there).
const GRIP0: f64 = 40.0;

#[test]
fn the_gripper_fingers_open_by_half_the_grip_and_meet_when_closed() {
    use popkern_assembly::{InterferenceOptions, PairStatus, SwingOptions};
    let mut asm = build();
    let placement = |asm: &Assembly, name: &str| {
        asm.tree()
            .iter()
            .find(|r| r.path.to_string() == name)
            .map(|r| r.world.transform().matrix)
            .unwrap_or_else(|| panic!("no {name}"))
    };
    let (a0, b0) = (placement(&asm, "finger_a"), placement(&asm, "finger_b"));
    for (q, grip) in [
        ([0.0; 6], 0.0),
        ([0.0; 6], 80.0),
        ([37.0, -71.0, 104.0, -122.0, 0.0, 215.0], 23.5),
    ] {
        pose(&mut asm, &q);
        asm.set_assembly_param("grip", grip).expect("grip");
        // Each finger moves along the jaw axis (the tool frame's x) by
        // half the grip change, finger B the other way; nothing else.
        let tool = {
            let t = frames(&table(), &q);
            mul(&t[6], &inverse(&frames(&table(), &[0.0; 6])[6]))
        };
        for (name, home, sign) in [("finger_a", a0, 1.0), ("finger_b", b0, -1.0)] {
            let shift = sign * (grip - GRIP0) / 2.0;
            let mut slid = home;
            for i in 0..3 {
                slid[i][3] += shift * [1.0, 0.0, 0.0][i];
            }
            let want = mul(&tool, &slid);
            let got = placement(&asm, name);
            for i in 0..3 {
                for j in 0..4 {
                    assert!(
                        (got[i][j] - want[i][j]).abs() < 1e-9,
                        "{name} at grip {grip}, pose {q:?}: [{i}][{j}] {} vs {}",
                        got[i][j],
                        want[i][j]
                    );
                }
            }
        }
    }
    pose(&mut asm, &[0.0; 6]);
    // Closed, the jaws meet face to face: touching, not interfering.
    asm.set_assembly_param("grip", 0.0).expect("grip");
    let report = asm
        .interference(&InterferenceOptions {
            among: vec!["finger_a".into(), "finger_b".into()],
            ..InterferenceOptions::default()
        })
        .expect("interference");
    assert!(
        report
            .pairs
            .iter()
            .all(|p| matches!(p.status, PairStatus::Touching { .. }))
            && !report.pairs.is_empty(),
        "closed jaws: {:?}",
        report.pairs
    );
    // Opening over the full travel meets nothing new.
    let swing = asm
        .swing_with(
            "grip",
            &SwingOptions {
                from: 0.0,
                to: 80.0,
                step: 5.0,
                resolution: 0.5,
                clearance: Some(0.5),
                exclude: Vec::new(),
            },
        )
        .expect("grip swing");
    assert_eq!(
        swing.first_contact, None,
        "grip swing: {:?}",
        swing.first_contact_pair
    );
}
