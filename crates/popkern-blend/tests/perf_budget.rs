//! §18 CI performance budget gate (PK-190).
//!
//! Absolute wall-clock budgets for the fillet pipeline's mainstream
//! paths. A regression past a budget **fails the build** — this is the
//! gate the 2026-04-13 state assessment asked for ("track in CI, fail
//! on regression beyond a margin"), complementing:
//!
//! - `tests/perf_acceptance.rs` — the §18 analytic ≥ 5× *ratio* gate
//!   (noise-robust, runs in every debug test sweep), and
//! - `benches/fillet_perf.rs` — criterion baselines for local,
//!   fine-grained tracking (`cargo bench -p popkern-blend`).
//!
//! # Budget justification (acceptance row 276)
//!
//! Each budget is `measured × headroom`, not a round figure. Measured
//! values are the **best-of-5** release-mode wall clock on the dev
//! reference machine (16-core linux, 2026-08-13, commit at PK-190
//! close); best-of-N is used because the minimum is the statistic
//! least sensitive to scheduler noise. The headroom multiplier is 5×:
//! GitHub-hosted runners measure 2–4× slower than the reference
//! machine on compute-bound work, and the gate must catch order-of-
//! magnitude regressions (a fallback path silently replacing a fast
//! path), not 20% drift — criterion owns fine-grained drift.
//!
//! Budgets are meaningless in debug builds (10–50× slower, different
//! optimizer): every test here self-skips outside release. CI runs
//! this file in the release-mode `perf` job; locally use
//! `tools/quality-check.sh --perf`.

use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::time::{Duration, Instant};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Skip guard: budgets only bind in release builds.
fn debug_skip(name: &str) -> bool {
    if cfg!(debug_assertions) {
        eprintln!(
            "§18 perf budget [{name}]: skipped in debug — enforced by the CI release perf job"
        );
        return true;
    }
    false
}

/// Best-of-N wall clock for `op`. A warmup call precedes the timed
/// runs except at `n == 1`, where the op is long enough (tens of
/// seconds) that cache warmup is noise and doubling the runtime would
/// cost more CI time than the precision buys.
fn best_of<F: FnMut()>(n: usize, mut op: F) -> Duration {
    if n > 1 {
        op();
    }
    let mut best = Duration::MAX;
    for _ in 0..n {
        let t0 = Instant::now();
        op();
        best = best.min(t0.elapsed());
    }
    best
}

/// Assert `measured` is within `budget_ms`, printing the number either
/// way so CI logs double as a coarse performance history.
fn assert_budget(name: &str, measured: Duration, budget_ms: u64) {
    let budget = Duration::from_millis(budget_ms);
    eprintln!("§18 perf budget [{name}]: best-of-run {measured:?} (budget {budget:?})");
    assert!(
        measured <= budget,
        "§18 performance regression: {name} took {measured:?}, budget {budget:?} \
         (= measured baseline × 5 headroom; see file header). A jump past this \
         budget means a fast path stopped applying or a phase went super-linear — \
         profile with `cargo bench -p popkern-blend` before touching the budget."
    );
}

fn cube_specs(cube: &BRep, count: usize) -> Vec<FilletSpec> {
    cube.edges
        .keys()
        .take(count)
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect()
}

/// Single-edge constant-radius cube fillet — the §11.1 analytic-torus
/// fast path plus full pipeline overhead. The cheapest end-to-end op.
#[test]
fn single_edge_cube_fillet_within_budget() {
    if debug_skip("single-edge cube fillet") {
        return;
    }
    let cube = BRep::unit_cube(tol()).expect("unit cube builds");
    let specs = cube_specs(&cube, 1);
    let measured = best_of(5, || {
        fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None)
            .expect("single-edge cube fillet succeeds");
    });
    // Measured 17.0 ms best-of-5 release (2026-08-13) × 5 ≈ 85 ms.
    assert_budget("single-edge cube fillet", measured, 85);
}

/// Full 12-edge cube fillet at the default G1 continuity — batch
/// analysis, chain detection, 8 three-edge corner patches, topology
/// apply. The mainstream "fillet the whole part" operation.
#[test]
fn twelve_edge_cube_fillet_within_budget() {
    if debug_skip("12-edge cube fillet (G1)") {
        return;
    }
    let cube = BRep::unit_cube(tol()).expect("unit cube builds");
    let specs = cube_specs(&cube, 12);
    let measured = best_of(5, || {
        fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None)
            .expect("12-edge cube fillet succeeds");
    });
    // Measured 186 ms best-of-5 release (2026-08-17, PK-298) × 5 ≈
    // 930 ms. Re-baselined from 85.4 ms (2026-08-14, PK-281): the
    // PK-294 endpoint-window fix un-crashed the n=3 PnS dispatch, so
    // each of the 8 corners now honestly constructs and measures a
    // PnS cap candidate (~10 ms) before declining to the sphere —
    // real acceptance work, where before the attempt died free on a
    // split error. The other 330 ms of that regression was the
    // PK-203 curve-off-surface audit running unconditionally; PK-298
    // gated it behind the pns-curve-audit trace topic (decision-neutral —
    // nothing reads its output).
    assert_budget("12-edge cube fillet (G1)", measured, 930);
}

/// Full 12-edge cube fillet through the G2 construction ladder
/// (rational Gregory rungs + KP setback/polar caps) — the heaviest
/// mainstream path and the one most likely to regress when ladder
/// rungs are reordered.
#[test]
fn twelve_edge_cube_g2_fillet_within_budget() {
    if debug_skip("12-edge cube fillet (G2 ladder)") {
        return;
    }
    let cube = BRep::unit_cube(tol()).expect("unit cube builds");
    let specs = cube_specs(&cube, 12);
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    // Single timed run: the op is ~10 s, so best-of-N would dominate
    // the whole perf job for precision the 5× headroom doesn't need.
    let measured = best_of(1, || {
        fillet_edges(&cube, &specs, &options, tol(), None)
            .expect("12-edge G2 cube fillet succeeds");
    });
    // Measured 2.7-2.8 s single-run release (2026-08-20, PK-278
    // parallel corner presolve + PK-301 strip-binder locality +
    // PK-302 parallel tail measurement; was 5.9 s after PK-289/291,
    // 10.0 s at PK-281, 38.5 s before that) × 5 ≈ 14 s.
    assert_budget("12-edge cube fillet (G2 ladder)", measured, 14000);
}

/// Cylinder base-circle fillet — the curved-support seam-crossing
/// path (plane×cylinder analytic torus on a periodic edge),
/// structurally different from the cube's plane-plane analytics.
#[test]
fn cylinder_base_fillet_within_budget() {
    if debug_skip("cylinder base fillet (curved support)") {
        return;
    }
    let cyl = BRep::cylinder(1.0, 2.0, tol()).expect("cylinder builds");
    let edge = cyl
        .edges
        .keys()
        .find(|&e| {
            let he = cyl.edges[e].half_edge;
            let twin = cyl.half_edges[he].twin;
            cyl.loops[cyl.half_edges[he].loop_].face != cyl.loops[cyl.half_edges[twin].loop_].face
        })
        .expect("cylinder has a base edge with distinct faces");
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let measured = best_of(5, || {
        let result = fillet_edges(&cyl, &specs, &FilletOptions::default(), tol(), None)
            .expect("cylinder base fillet succeeds");
        assert!(result.brep.validate().is_ok(), "result must validate");
    });
    // Measured 11–20 ms release (2026-08-13) × 5 ≈ 100 ms.
    assert_budget("cylinder base fillet (curved support)", measured, 100);
}

/// Cone base-circle fillet — the phase-1b walker path (no analytic
/// offset exists for cones, so the rolling-ball contact walker
/// produces the center curve; the §11.1 plane-cone arm then promotes
/// it to an exact torus). Restored by PK-264: this op failed at every
/// radius until the convexity mirror-ball and footprint-clip defects
/// were fixed at the root, and its previous bench entry silently
/// timed the failure path.
#[test]
fn cone_base_fillet_within_budget() {
    if debug_skip("cone base fillet (walker)") {
        return;
    }
    let cone = BRep::cone(1.0, 2.0, tol()).expect("cone builds");
    let edge = cone
        .edges
        .keys()
        .find(|&e| {
            let he = cone.edges[e].half_edge;
            let twin = cone.half_edges[he].twin;
            cone.loops[cone.half_edges[he].loop_].face
                != cone.loops[cone.half_edges[twin].loop_].face
        })
        .expect("cone has a base edge with distinct faces");
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.03),
    }];
    let measured = best_of(5, || {
        let result = fillet_edges(&cone, &specs, &FilletOptions::default(), tol(), None)
            .expect("cone base fillet succeeds (PK-264)");
        assert!(result.brep.validate().is_ok(), "result must validate");
    });
    // Measured 15.6–16.6 ms release (2026-08-13) × 5 ≈ 85 ms.
    assert_budget("cone base fillet (walker)", measured, 85);
}
