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
    // A folded elbow, j3 at 178 (both elbow branches past its 170
    // limit): refused, never served past the limit.
    let folded = [0.0, 90.0, -90.0, 180.0]
        .iter()
        .flat_map(|&q4| [0.0, 90.0, -90.0, 180.0].map(|q5| (q4, q5)))
        .find_map(|(q4, q5)| reachable_from([30.0, -40.0, 178.0, q4, q5, 0.0]))
        .expect("a folded target");
    let mut asm = skeleton_with(offset(folded, home), true);
    let start = placement(&asm, "link_6");
    let refused = asm.set_assembly_param("s", 1.0);
    eprintln!("[pk808] folded elbow: {refused:?}");
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
