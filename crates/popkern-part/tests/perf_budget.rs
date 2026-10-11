//! CI performance budget gate for the smooth-mesh path (PK-190).
//!
//! Companion to `popkern-blend/tests/perf_budget.rs` — see its header
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

/// The budgets in this file run one at a time. libtest runs them in
/// parallel by default, and the CI runner has two cores: on CI the
/// cage tests overlapped both thread-cut builds unevenly, so the ratio
/// measured contention as well as scaling. PK-786 measured 3.54 on CI
/// (len 40 = 25.4 s, len 80 = 90.0 s). After its fixes, serialized
/// runs measured 2.77-2.80 locally (len 40 = 9.0 s, len 80 = 25.1 s);
/// after PK-794's, 1.90-1.92 (len 40 = 4.6 s, len 80 = 8.8 s).
static SERIAL: std::sync::Mutex<()> = std::sync::Mutex::new(());

fn serial() -> std::sync::MutexGuard<'static, ()> {
    SERIAL
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

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
    if !toms_tables_or_skip() {
        return;
    }
    if debug_skip("smooth-mesh cube cage") {
        return;
    }
    let _serial = serial();
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
    if !toms_tables_or_skip() {
        return;
    }
    if debug_skip("smooth-mesh truncated-corner cage (1 CC step)") {
        return;
    }
    let _serial = serial();
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

/// Multi-wrap thread-cut serve cost stays SUB-QUADRATIC in tool
/// turns (PK-496). The cellular merge stage classified O(turns)
/// wrap-band fragments with independent point-in-solid rays over
/// multi-wrap walls — N² — until A-side verdict propagation
/// (one audited seed ray per component) and the span-hull projection
/// pre-gate landed. This gate is a RATIO, not an absolute budget:
/// doubling the turn count may at most triple the build plus mesh.
/// Measured 2.5 after PK-511's pierce-inventory cache (2026-09-23).
/// PK-794 (2026-10-11) made the remaining stages linear (Delaunay
/// insertion order, pinch search, orientation votes, Bezier
/// decomposition, one face build per tessellation), and it measures
/// 1.90-1.92. The gate moved from 3.5 to 3.0: at 3.5 the mesh stage
/// had grown 3.3x per doubling underneath it unnoticed. A quadratic
/// stage measures >= 3.9. Ratios stay meaningful under uniform CI
/// contention, unlike wall-clock budgets.
#[test]
fn thread_cut_serve_cost_stays_subquadratic() {
    use popkern_part::fastener::{ThreadSpec, ThreadedRod};
    if debug_skip("thread-cut serve cost ratio (40 vs 80)") {
        return;
    }
    let _serial = serial();
    let build = |len: f64| {
        let t0 = Instant::now();
        let spec = ThreadSpec::parse("M6").expect("M6");
        let rod = ThreadedRod::new(spec, len)
            .tip_chamfer(false)
            .build(Tolerance::default())
            .expect("rod");
        assert_eq!(rod.open_edge_count().expect("mesh"), 0);
        t0.elapsed()
    };
    let t40 = build(40.0);
    let t80 = build(80.0);
    let ratio = t80.as_secs_f64() / t40.as_secs_f64().max(1e-9);
    eprintln!(
        "perf budget [thread-cut ratio]: len 40 = {t40:?}, len 80 = {t80:?}, ratio {ratio:.2} (gate 3.0)"
    );
    assert!(
        ratio <= 3.0,
        "thread-cut cost curve regressed toward quadratic: \
         len 40 = {t40:?}, len 80 = {t80:?}, ratio {ratio:.2} > 3.0 \
         (measured 1.9 at PK-794 close; quadratic >= 3.9)"
    );
}

/// The PnS mask tables are optional (ACM ToMS package, untracked —
/// PK-314/PK-315). Tests that build smooth-mesh bodies need them at
/// runtime; on checkouts without the package they skip with a visible
/// notice instead of failing.
fn toms_tables_or_skip() -> bool {
    if popkern_core::geometry::polyhedral_spline::toms_tables_available() {
        return true;
    }
    eprintln!("[toms] ACM ToMS 1032 package absent — smooth-mesh test skipped (PK-314/PK-315)");
    false
}
