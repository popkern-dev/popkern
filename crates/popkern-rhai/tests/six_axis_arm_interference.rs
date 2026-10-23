//! Whole-arm interference (PK-847): every occurrence pair of the
//! six-axis arm is touching or clear. It reported 129 interfering pairs
//! in 626 s: ~126 screws whose real threads overlapped their tap-drill
//! bores by design (the bores are real tapped holes now and every screw
//! phased to mate, PK-876; turned the wrong way, all 8 of j4's mount
//! screws interfere), the
//! base against the j1 housing (a dome cap read inside out by the trim
//! containment), and two tubes butted into their end castings (a closed
//! record along both faces' rims; web-plane trim residues counted as
//! crossings) whose overlap booleans declined.

mod arm_support;

use arm_support::*;
use popkern_assembly::{InterferenceOptions, PairStatus};

#[test]
fn the_whole_arm_has_no_interfering_pair() {
    let mut asm = build();
    let _ = popkern_bool::chart_walk_fallback::take_chart_walk_fallback_events();
    let watch = std::time::Instant::now();
    let report = asm
        .interference(&InterferenceOptions::default())
        .expect("interference");
    let took = watch.elapsed().as_secs_f64();
    // A link face whose 42.5 rim lies on a coaxial bore wall of the same
    // radius walked the 3D-Newell fallback on every pass: the records on
    // the rim registered in both directions (PK-920). Each interference
    // boolean now walks every face on its chart.
    let fallbacks: Vec<String> =
        popkern_bool::chart_walk_fallback::take_chart_walk_fallback_events()
            .into_iter()
            .map(|e| format!("{:?} a={}: {}", e.face, e.operand_a, e.reason))
            .collect();
    assert!(fallbacks.is_empty(), "3D-Newell fallbacks: {fallbacks:#?}");
    let interfering: Vec<String> = report
        .pairs
        .iter()
        .filter(|p| matches!(p.status, PairStatus::Interfering { .. }))
        .map(|p| format!("{} x {}: {:?}", p.a, p.b, p.status))
        .collect();
    assert!(
        interfering.is_empty(),
        "{} interfering: {interfering:#?}",
        interfering.len()
    );
    assert!(
        report.occurrences > 700,
        "{} occurrences",
        report.occurrences
    );
    // 18.3 s quiet at the test profile with cosmetic threads (626 s at
    // filing). Real mated threads measured 135.9 s, then 42.7 s once
    // certification re-meshed only each pair's near faces (56.5 s on
    // whole bodies) and the tap flanks meshed within their chord (a
    // flank sampled by phase luck sagged 16x it, and the gripper screws
    // descended to the finest level, 23 s each; PK-900).
    assert!(took < 60.0, "whole-arm interference took {took:.1} s");
}
