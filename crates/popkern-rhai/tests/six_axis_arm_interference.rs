//! Whole-arm interference (PK-847): every occurrence pair of the
//! six-axis arm is touching or clear. It reported 129 interfering pairs
//! in 626 s: ~126 screws whose real threads overlapped their tap-drill
//! bores by design (now cosmetic threads, PK-876 tracks real ones), the
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
    let watch = std::time::Instant::now();
    let report = asm
        .interference(&InterferenceOptions::default())
        .expect("interference");
    let took = watch.elapsed().as_secs_f64();
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
    // 18.3 s quiet at the test profile (626 s at filing).
    assert!(took < 60.0, "whole-arm interference took {took:.1} s");
}
