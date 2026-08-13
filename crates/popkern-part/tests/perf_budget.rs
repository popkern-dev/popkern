//! CI performance budget gate for the smooth-mesh path (PK-190).
//!
//! Companion to `popkern-fillet/tests/perf_budget.rs` — see its header
//! for the budget-justification method (best-of-N release wall clock
//! on the dev reference machine × 5 CI headroom). This file gates the
//! PnS mesh-to-BRep smoke tier: the cube cage (direct serve, the
//! minimal full pipeline) and the truncated-corner cage (one
//! Catmull–Clark normalization step + retry — the PK-198 mechanism).
//! The heavy fixtures (suzanne/hand/airplane) are budgeted separately
//! as the nightly workflow's job timeout, not here.
//!
//! Budgets are meaningless in debug builds: every test self-skips
//! outside release. CI runs this in the release `perf` job; locally
//! use `tools/quality-check.sh --perf`.

use popkern_core::tolerance::Tolerance;
use popkern_part::{SmoothMeshOptions, smooth_solid_from_obj};
use std::time::{Duration, Instant};

/// Closed unit-cube control cage (8 vertices, 6 quads) — serves
/// directly with no refinement step.
const CUBE_CAGE: &str = "\
v -1 -1 -1
v  1 -1 -1
v  1  1 -1
v -1  1 -1
v -1 -1  1
v  1 -1  1
v  1  1  1
v -1  1  1
f 1 4 3 2
f 5 6 7 8
f 1 2 6 5
f 2 3 7 6
f 3 4 8 7
f 4 1 5 8
";

/// Cube with corner (1,1,1) truncated to a triangle through the three
/// adjacent edge midpoints — pentagons and mixed valences decline the
/// direct serve, so it goes through one Catmull–Clark normalization
/// step (the PK-198 retry path). Same cage as
/// `smooth_mesh_from_cage.rs::truncated_corner_cage`, spelled as OBJ.
const TRUNCATED_CORNER_CAGE: &str = "\
v -1 -1 -1
v  1 -1 -1
v  1  1 -1
v -1  1 -1
v -1 -1  1
v  1 -1  1
v -1  1  1
v  1  0  1
v  1  1  0
v  0  1  1
f 1 4 3 2
f 5 6 8 10 7
f 1 2 6 5
f 2 3 9 8 6
f 3 4 7 10 9
f 4 1 5 7
f 8 9 10
";

fn debug_skip(name: &str) -> bool {
    if cfg!(debug_assertions) {
        eprintln!("perf budget [{name}]: skipped in debug — enforced by the CI release perf job");
        return true;
    }
    false
}

fn best_of<F: FnMut()>(n: usize, mut op: F) -> Duration {
    op();
    let mut best = Duration::MAX;
    for _ in 0..n {
        let t0 = Instant::now();
        op();
        best = best.min(t0.elapsed());
    }
    best
}

fn assert_budget(name: &str, measured: Duration, budget_ms: u64) {
    let budget = Duration::from_millis(budget_ms);
    eprintln!("perf budget [{name}]: best-of-run {measured:?} (budget {budget:?})");
    assert!(
        measured <= budget,
        "performance regression: {name} took {measured:?}, budget {budget:?} \
         (= measured baseline × 5 headroom; see file header)"
    );
}

/// Cube cage through the full smooth-mesh pipeline: construction,
/// weld, transactional repair, validation, `Part` wrap.
#[test]
fn smooth_mesh_cube_cage_within_budget() {
    if debug_skip("smooth-mesh cube cage") {
        return;
    }
    let measured = best_of(3, || {
        let (_, report) = smooth_solid_from_obj(
            CUBE_CAGE,
            Tolerance::default(),
            SmoothMeshOptions::default(),
        )
        .expect("cube cage serves");
        assert_eq!(report.refinement_steps, 0, "cube must serve directly");
    });
    // Measured 3.86 s best-of-3 release on the dev reference machine
    // (2026-08-13, PK-190 close) × 5 headroom ≈ 20 s.
    assert_budget("smooth-mesh cube cage", measured, 20000);
}

/// Truncated-corner cage: one Catmull–Clark refinement plus the full
/// retry — the cost a caller pays for a cage outside the direct
/// serving families.
#[test]
fn smooth_mesh_normalized_cage_within_budget() {
    if debug_skip("smooth-mesh truncated-corner cage (1 CC step)") {
        return;
    }
    let measured = best_of(3, || {
        let (_, report) = smooth_solid_from_obj(
            TRUNCATED_CORNER_CAGE,
            Tolerance::default(),
            SmoothMeshOptions::default(),
        )
        .expect("truncated-corner cage serves after normalization");
        assert_eq!(report.refinement_steps, 1, "must serve via one CC step");
    });
    // Measured 16.26 s best-of-3 release on the dev reference machine
    // (2026-08-13, PK-190 close) × 5 headroom ≈ 85 s.
    assert_budget(
        "smooth-mesh truncated-corner cage (1 CC step)",
        measured,
        85000,
    );
}
