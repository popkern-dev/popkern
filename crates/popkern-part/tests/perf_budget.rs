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

/// A lone M5x12 socket cap screw, tip chamfer on (PK-799). The chamfer
/// cone crosses the thread flanks; before the radial level-curve arm
/// served cone x flank pairs the marcher paid a multi-start projection
/// per landing past the flanks' clamped ends, and the chamfer alone
/// cost 0.8 s of a 2.05 s screw (the rod gate above runs without it,
/// so the class was never measured). Two gates: the chamfer's
/// marginal cost as a RATIO against the same screw without it (robust
/// to uniform CI contention; build-only best-of-3 measured 1.12-1.13
/// with the arm, 1.84-1.94 with the cone arm disabled), and an
/// absolute budget of the measured best-of-3 x 5.
#[test]
fn lone_socket_screw_tip_chamfer_within_budget() {
    use popkern_part::fastener::{SocketCapScrew, ThreadSpec};
    if debug_skip("lone M5x12 socket cap screw") {
        return;
    }
    let _serial = serial();
    // Build time only; watertightness is checked once, untimed.
    let build = |chamfer: bool| {
        let mut screw = None;
        let best = best_of(3, || {
            let spec = ThreadSpec::parse("M5").expect("M5");
            screw = Some(
                SocketCapScrew::new(spec, 12.0)
                    .tip_chamfer(chamfer)
                    .build(Tolerance::default())
                    .expect("screw"),
            );
        });
        let screw = screw.expect("built");
        assert_eq!(screw.open_edge_count().expect("mesh"), 0);
        best
    };
    let plain = build(false);
    let chamfered = build(true);
    let ratio = chamfered.as_secs_f64() / plain.as_secs_f64().max(1e-9);
    eprintln!(
        "perf budget [tip chamfer ratio]: plain {plain:?}, chamfered {chamfered:?}, ratio {ratio:.2} (gate 1.4)"
    );
    assert!(
        ratio <= 1.4,
        "the tip chamfer's cone x flank pairs regressed: plain {plain:?}, \
         chamfered {chamfered:?}, ratio {ratio:.2} > 1.4 (1.9 on the marcher)"
    );
    // Measured 0.94-0.97 s best-of-3 release (2026-10-12, PK-799) x 5.
    assert_budget("lone M5x12 socket cap screw", chamfered, 5000);
}

/// Twelve M5 taps on a housing's bolt circle (the six-axis arm's R5
/// fixture, PK-744): the twelfth tap must cost about what it costs on
/// the bare housing. Every tap's booleans used to re-measure the
/// earlier taps' threads: check 9 ran the cold multisheet projection
/// on every sample of every edge before trying the edge's own pcurve,
/// and edge overlap detection solved closest points for every pair of
/// operand edges. The gate times the SAME tap (best of 3) on the bare
/// housing and on the housing holding the other eleven, so only the
/// accumulated body differs: a RATIO, robust to uniform CI
/// contention. Measured 1.61-1.63; 2.29-2.74 with check 9's cold
/// search always run (the gate's negative control). Without the
/// edge-pair prescreen it reads 1.68-1.70, the twelfth tap 2.06-2.22 s
/// against 1.87 s: real, but under the ratio's resolution, so that
/// mechanism is pinned by its conservative-box tests instead.
#[test]
fn tapped_housing_taps_cost_what_the_first_did() {
    use popkern_part::fastener::{HoleDepth, TapOptions, ThreadSpec};
    use popkern_part::prelude::*;
    if debug_skip("tapped housing per-tap growth") {
        return;
    }
    let _serial = serial();
    let spec = ThreadSpec::parse("M5").expect("M5");
    let housing = Profile::on(Plane::XY)
        .circle(55.0)
        .extrude(30.0)
        .expect("housing");
    let tap = |part: &popkern_part::Part, k: usize| {
        let a = std::f64::consts::TAU * (k as f64 + 0.5) / 12.0;
        part.tapped_hole(
            &spec,
            popkern_core::geometry::Point3::new(42.0 * a.cos(), 42.0 * a.sin(), 30.0),
            popkern_core::geometry::Vec3::new(0.0, 0.0, -1.0),
            HoleDepth::Blind(10.0),
            &TapOptions::default(),
        )
        .expect("tap")
    };
    let mut eleven = housing.clone();
    for k in 0..11 {
        eleven = tap(&eleven, k);
    }
    let bare = best_of(3, || {
        tap(&housing, 11);
    });
    let mut twelve = None;
    let last = best_of(3, || {
        twelve = Some(tap(&eleven, 11));
    });
    let ratio = last.as_secs_f64() / bare.as_secs_f64().max(1e-9);
    eprintln!(
        "perf budget [tapped housing growth]: bare {bare:?}, twelfth {last:?}, ratio {ratio:.2} (gate 2.0)"
    );
    assert!(
        ratio <= 2.0,
        "the twelfth tap costs {ratio:.2}x the same tap on the bare housing \
         ({last:?} vs {bare:?}) > 2.0: tap cost grows with the body"
    );
    let twelve = twelve.expect("built");
    assert_eq!(twelve.open_edge_count().expect("mesh"), 0, "watertight");
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
