//! Target posing of the arm by mates (six-axis-arm §11 E4, PK-726):
//! the arm's kinematic skeleton — seven links on the §4.1 DH axes,
//! free revolutes (concentric + axial seat, no angle), the flange mated
//! concentric + coincident to a target that three distance mates place,
//! all driven by one parameter `s` (home → target in one jump).
//!
//! The measured classes are the acceptance (E4's ≥ 99% bar):
//!
//! - **Converged** from home in one jump: ≥ 99 of 100 seeded reachable
//!   targets, each on its target to 1e-9 mm. Measured 100: 85 with
//!   Gauss–Newton alone, 90 with the Levenberg–Marquardt retry, 95 with
//!   the revolute chain's joint-space solve (PK-763), 100 with its
//!   branch seeds for a jump (each joint a half turn over).
//! - **Solver-basin misses**: targets a damped least-squares oracle in
//!   joint space reaches from home but the mate solve does not: none.
//! - **Branch misses**: the remainder, reachable only on another IK
//!   branch.
//! - Every miss is refused transactionally; an unreachable target is a
//!   `MateConflict` with its residual and nothing moves.
//! - **Continuity**: 5 mm target moves never flip a branch (no joint
//!   turns over 10° per step): a request under 2% of the model is a
//!   continuous edit and branch seeds stay off (with them on, a 5 mm
//!   step at a branch's reach turned a joint 112°).

mod arm_support;
mod arm_target_support;

use arm_support::*;
use arm_target_support::*;
use popkern_assembly::AssemblyError;

#[test]
fn the_flange_reaches_seeded_targets_by_mates_in_the_measured_classes() {
    let home = col(&frames(&table(), &[0.0; 6])[6], 3);
    let (mut converged, mut basin_misses, mut branch_misses) = (0, Vec::new(), Vec::new());
    let mut solved = Vec::new();
    for (i, p) in targets().into_iter().enumerate() {
        let mut asm = skeleton(offset(p, home));
        let start = placement(&asm, "link_6");
        match asm.set_assembly_param("s", 1.0) {
            Ok(_) => {
                let at = apply(&placement(&asm, "link_6"), home);
                assert!(dist(at, p) < 1e-9, "target {i}: flange {at:?} vs {p:?}");
                converged += 1;
                solved.push((i, p));
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
                if oracle_reaches(&[0.0; 6], p) {
                    basin_misses.push(i);
                } else {
                    branch_misses.push(i);
                }
            }
        }
    }
    eprintln!(
        "[pk726] converged {converged}/100; solver-basin misses {basin_misses:?}; \
         branch misses {branch_misses:?}"
    );
    assert!(converged >= 99, "converged {converged}/100");
    assert!(
        basin_misses.is_empty(),
        "solver-basin misses {basin_misses:?}"
    );
    // The Levenberg–Marquardt retry's fixture: Gauss–Newton alone
    // stalls on target 26 (a high target over the base).
    assert!(solved.iter().any(|(i, _)| *i == 26), "target 26 regressed");

    // Continuity: 5 mm moves outward along the approach from the first
    // 30 that converged; each step either solves without a branch flip
    // or leaves the workspace (the oracle cannot reach it either).
    for &(i, p) in solved.iter().take(30) {
        let delta = offset(p, home);
        let length = dist(p, home);
        let mut asm = skeleton(delta);
        asm.set_assembly_param("s", 1.0).expect("solved before");
        let mut q0 = joints(&asm);
        for k in 1..=20 {
            let s = 1.0 + 5.0 * k as f64 / length;
            if asm.set_assembly_param("s", s).is_err() {
                let next = [0, 1, 2].map(|c| home[c] + s * delta[c]);
                assert!(
                    !oracle_reaches(&q0, next),
                    "target {i} step {k}: the solve lost a reachable target"
                );
                break;
            }
            let q1 = joints(&asm);
            let turn = q0
                .iter()
                .zip(&q1)
                .map(|(a, b)| ((b - a + 540.0).rem_euclid(360.0) - 180.0).abs())
                .fold(0.0, f64::max);
            assert!(turn < 10.0, "target {i} step {k}: a joint turned {turn}°");
            q0 = q1;
        }
    }

    // Unreachable: two metres out is refused with its residual, and
    // nothing moves.
    let mut far = skeleton([2000.0, 0.0, 0.0]);
    let before = placement(&far, "link_6");
    match far.set_assembly_param("s", 1.0) {
        Err(AssemblyError::MateConflict { report }) => {
            assert!(!report.converged && report.residual > 1e-3, "{report:?}");
            assert!(!report.conflicting.is_empty());
        }
        other => panic!("an unreachable target: {other:?}"),
    }
    assert_eq!(placement(&far, "link_6"), before);
}
