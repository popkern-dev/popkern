//! PK-600: the KP20 setback cap's per-corner assembly cost at n = 8 on
//! the production border (the lifted Catmull–Clark ring, genuine guide
//! misfit), warm: tau and the KKT pattern are memoised per valence, so
//! the first corner of a valence pays more.
//!
//! A release-only perf budget in the PK-272 pattern (it self-skips in
//! debug; `tools/quality-check.sh --perf` runs it):
//!
//! ```text
//! cargo test --release -p popkern-core --test kp20_cap_perf_budget -- --nocapture
//! ```
//!
//! # Budget
//!
//! `measured × 5 headroom`, best of 5, 16-core linux dev reference
//! machine, 2026-10-21: 0.34 s release (0.38 s in the test profile,
//! against the ticket's 0.5 s). PK-604 had it at 0.54 s in the test
//! profile; PK-676's deterministic linear algebra took it to 0.75 s
//! unnoticed, as no gate covered it. Three mechanisms brought it back:
//! a supernodal `L D Lᵀ` for the coupled curvature KKT (0.45 s ->
//! 0.10 s; dense panels, relative scatter), one elimination for the
//! guide's three coordinates (0.17 s -> 0.05 s, bit-identical), and
//! the frozen curvature Grams formed in parallel (0.18 s -> 0.03 s,
//! identical).

use std::time::{Duration, Instant};

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Kp20PaperGuideBorder, Kp20PaperGuideBorderPiece, assemble_kp20_paper_cap,
    cc_lifted_ring_border_pieces,
};

fn production_border(n: usize) -> Kp20PaperGuideBorder {
    let pieces =
        cc_lifted_ring_border_pieces(n, |x, y| 0.4 * (x * x + y * y)).expect("lifted ring");
    Kp20PaperGuideBorder {
        sectors: pieces
            .into_iter()
            .map(|[near_c, near_m]| {
                [
                    Kp20PaperGuideBorderPiece::from_bi3(near_c),
                    Kp20PaperGuideBorderPiece::from_bi3(near_m),
                ]
            })
            .collect(),
    }
}

#[test]
fn kp20_cap_at_valence_eight_within_budget() {
    let border = production_border(8);
    let assemble = || assemble_kp20_paper_cap(&border, Point3::new(0.0, 0.0, 0.0));
    // Every build: the cap assembles and the bi-9 resampling serves.
    let report = assemble().expect("assemble");
    let resample = report.bi9.as_ref().expect("bi-9 resample");
    assert_eq!(resample.pieces.len(), 32);
    // Both fair-interior candidates solve (the cap keeps the fairer of
    // them and the resampled start).
    assert!(
        resample
            .fair
            .attempts
            .iter()
            .flatten()
            .all(|a| a.max_curvature_after.is_some()),
        "{:?}",
        resample.fair.attempts
    );
    if cfg!(debug_assertions) {
        eprintln!("KP20 n = 8 cap budget skipped in debug: enforced by the release perf job");
        return;
    }
    let mut best = Duration::MAX;
    for _ in 0..5 {
        let t0 = Instant::now();
        assemble().expect("assemble");
        best = best.min(t0.elapsed());
    }
    let budget = Duration::from_millis(1700);
    eprintln!("KP20 n = 8 cap: best of 5 {best:?} (budget {budget:?})");
    assert!(
        best <= budget,
        "KP20 n = 8 cap assembly took {best:?}, budget {budget:?} (measured 0.34 s × 5)"
    );
}
