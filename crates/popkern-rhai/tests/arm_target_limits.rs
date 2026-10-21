//! Released joints hold the arm's limits under target posing (PK-808,
//! six-axis-arm §11 E4): on the skeleton of `arm_target_posing` with
//! its revolutes as released joints, every seeded target solves within
//! the limits or is refused with nothing moved, and a folded elbow is
//! refused.

mod arm_support;
mod arm_target_support;

use arm_support::*;
use arm_target_support::*;
use popkern_assembly::AssemblyError;

#[test]
fn released_joints_hold_their_limits_over_the_seeded_targets_and_the_walk() {
    // The same targets, the revolutes now released joints with the
    // arm's limits: every solve either keeps every joint within them or
    // is refused (a limit, or the mates) with nothing moved.
    let home = col(&frames(&table(), &[0.0; 6])[6], 3);
    let within = |q: &[f64; 6]| {
        q.iter()
            .zip(&LIMITS)
            .all(|(a, l)| *a >= l[0] - 1e-6 && *a <= l[1] + 1e-6)
    };
    let (mut solved, mut limited) = (Vec::new(), 0);
    for (i, p) in targets().into_iter().enumerate() {
        let mut asm = skeleton_with(offset(p, home), true);
        let start = placement(&asm, "link_6");
        match asm.set_assembly_param("s", 1.0) {
            Ok(_) => {
                let q = joints(&asm);
                assert!(within(&q), "target {i}: joints {q:?} past their limits");
                // The released parameters read back the solved angles.
                let j3 = asm
                    .params()
                    .iter()
                    .find(|x| x.name == "j3")
                    .expect("j3")
                    .display_value();
                assert!(
                    (j3 - q[2]).abs() < 1e-6,
                    "target {i}: j3 reads {j3} vs {}",
                    q[2]
                );
                solved.push((i, p));
            }
            Err(AssemblyError::JointLimit { joint, angle, .. }) => {
                assert_eq!(
                    placement(&asm, "link_6"),
                    start,
                    "target {i} moved on refusal"
                );
                assert!(joint == "j3", "target {i}: {joint} at {angle}");
                limited += 1;
            }
            Err(error) => {
                assert!(
                    matches!(error, AssemblyError::MateConflict { .. }),
                    "target {i}: {error}"
                );
                assert_eq!(
                    placement(&asm, "link_6"),
                    start,
                    "target {i} moved on refusal"
                );
            }
        }
    }
    eprintln!(
        "[pk808] {} solved within limits, {limited} refused at a limit",
        solved.len()
    );
    // A folded elbow, built at j3 = 178°, past its 170° limit. Another
    // branch reaches it at j3 = 153° (PK-902: refused before, the
    // nearest branch from home lying past the limit), checked against
    // the DH chain; with j3 held to ±100° no branch fits and it is
    // refused, nothing moved.
    let folded = [0.0, 90.0, -90.0, 180.0]
        .iter()
        .flat_map(|&q4| [0.0, 90.0, -90.0, 180.0].map(|q5| (q4, q5)))
        .find_map(|(q4, q5)| reachable_from([30.0, -40.0, 178.0, q4, q5, 0.0]))
        .expect("a folded target");
    let mut asm = skeleton_with(offset(folded, home), true);
    asm.set_assembly_param("s", 1.0)
        .expect("a branch within the limits");
    let q = joints(&asm);
    assert!(within(&q), "folded: joints {q:?}");
    assert!(
        dist(col(&frames(&table(), &q)[6], 3), folded) < 1e-6,
        "folded: {q:?}"
    );
    let mut tight = LIMITS;
    tight[2] = [-100.0, 100.0];
    let mut asm = skeleton_with_limits(offset(folded, home), true, &tight);
    let start = placement(&asm, "link_6");
    let refused = asm.set_assembly_param("s", 1.0);
    eprintln!("[pk808] folded elbow, j3 within ±100°: {refused:?}");
    assert!(
        matches!(refused, Err(AssemblyError::JointLimit { ref joint, .. }) if joint == "j3")
            || matches!(refused, Err(AssemblyError::MateConflict { .. })),
        "{refused:?}"
    );
    assert_eq!(placement(&asm, "link_6"), start);
    assert!(solved.len() + limited >= 95, "{} + {limited}", solved.len());
    for &(i, p) in solved.iter().take(30) {
        let delta = offset(p, home);
        let length = dist(p, home);
        let mut asm = skeleton_with(delta, true);
        asm.set_assembly_param("s", 1.0).expect("solved before");
        for k in 1..=20 {
            let s = 1.0 + 5.0 * k as f64 / length;
            if asm.set_assembly_param("s", s).is_err() {
                break;
            }
            let q = joints(&asm);
            assert!(
                within(&q),
                "target {i} step {k}: joints {q:?} past their limits"
            );
        }
    }
}

#[test]
fn a_target_with_a_branch_within_the_limits_is_served_on_it() {
    // PK-902: targets built from configurations within the limits, the
    // elbow near its 170° stop. From home the nearest branch of three of
    // them lay past it (j3 at 170.3°, −179.9°, 171.4°) and the target
    // was refused though the configuration it came from is a solution;
    // the joint-space rescue now ranks the branches by their limits.
    let home = col(&frames(&table(), &[0.0; 6])[6], 3);
    let within = |q: &[f64; 6]| {
        q.iter()
            .zip(&LIMITS)
            .all(|(a, l)| *a >= l[0] - 1e-6 && *a <= l[1] + 1e-6)
    };
    let mut served = Vec::new();
    for q1 in [-120.0, -60.0, 0.0, 60.0, 120.0] {
        for q2 in [-80.0, -40.0, 0.0, 40.0, 80.0] {
            for q3 in [150.0, 160.0, 165.0, 169.0, -150.0, -160.0, -165.0, -169.0] {
                let Some(p) = [0.0, 90.0]
                    .iter()
                    .find_map(|&q4| reachable_from([q1, q2, q3, q4, 0.0, 0.0]))
                else {
                    continue;
                };
                let mut asm = skeleton_with(offset(p, home), true);
                asm.set_assembly_param("s", 1.0)
                    .unwrap_or_else(|e| panic!("{:?}: {e}", [q1, q2, q3]));
                let q = joints(&asm);
                assert!(within(&q), "{:?}: joints {q:?}", [q1, q2, q3]);
                served.push((p, asm));
            }
        }
    }
    assert!(served.len() >= 30, "{} near-limit targets", served.len());
    // A continuous edit from there keeps its branch: 5 mm steps stay
    // within the limits or are refused, never flipped past them.
    for (p, asm) in served.iter_mut().take(10) {
        let length = dist(*p, home);
        for k in 1..=4 {
            let s = 1.0 + 5.0 * k as f64 / length;
            if asm.set_assembly_param("s", s).is_err() {
                break;
            }
            let q = joints(asm);
            assert!(within(&q), "step {k}: joints {q:?}");
        }
    }
}
