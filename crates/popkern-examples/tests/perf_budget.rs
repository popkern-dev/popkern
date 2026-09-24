//! Offset/shell performance budgets (PK-329 row 485, PK-190 pattern).
//!
//! Each budget is `measured × headroom`, not a round figure. Measured
//! values are the best-of-5 release-mode wall clock on the dev
//! reference machine (16-core linux, idle). The headroom multiplier is
//! 5×: GitHub-hosted runners measure 2–4× slower, and the gate must
//! catch order-of-magnitude regressions (a fallback path silently
//! replacing a fast path), not 20% drift. Each row also asserts a
//! correctness invariant so the timing can never be of a failure path
//! (PK-190 honesty rule).
//!
//! Baselines (2026-08-27, post-PK-329 campaign): banded least-squares
//! curve fitting (240 ms → 4 ms per marched SSI fit), analytic pair
//! separation certificates (concentric-skin pairs certified O(1)),
//! and rayon ordered-collect dispatch across the pipeline. The
//! pre-campaign corpus numbers these rows protect: wavy 5405 ms → 170
//! ms, rod 256 ms → 7 ms, bracket 53 ms → 10 ms, filleted cube 216 ms
//! → 57 ms.

use popkern_core::tolerance::Tolerance;
use popkern_examples::offset_corpus::{
    corpus, filleted_cube, l_bracket, sphere_capped_rod, wavy_plate,
};
use popkern_shell::thicken::thicken;
use popkern_shell::{ShellDirection, ShellOptions, max_safe_offset};
use std::time::{Duration, Instant};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Perf gates only mean something in release mode; a debug run
/// self-skips (the workspace sweep runs debug, CI runs this file with
/// `--release` in the perf job).
fn debug_skip(name: &str) -> bool {
    if cfg!(debug_assertions) {
        eprintln!("[perf-budget] skipping '{name}' (debug build)");
        return true;
    }
    false
}

fn best_of<F: FnMut()>(n: usize, mut op: F) -> Duration {
    op(); // warmup
    (0..n)
        .map(|_| {
            let start = Instant::now();
            op();
            start.elapsed()
        })
        .min()
        .expect("at least one run")
}

/// Load-robust cost ratio `time(b) / time(a)`: the median over `n`
/// back-to-back `(a, b)` pairs of their pair ratio, plus each side's
/// best time for the report. Each pair sees the same machine state,
/// and a preemption burst spoils only the pairs it lands on, which the
/// median ignores while they stay a minority. Timing one side's
/// samples after the other's let load changes skew the ratio
/// (measured under a full-workspace nextest run: cube 4.9 ms then
/// l_bracket 12.2 ms, ratio 2.50 against a 1.3 nominal), and a
/// best-of-15 interleave still read 2.84 once with every core busy
/// (PK-526).
fn median_pair_ratio<A: FnMut(), B: FnMut()>(
    n: usize,
    mut op_a: A,
    mut op_b: B,
) -> (f64, Duration, Duration) {
    op_a(); // warmup
    op_b();
    let (mut best_a, mut best_b) = (Duration::MAX, Duration::MAX);
    let mut ratios: Vec<f64> = (0..n)
        .map(|_| {
            let start = Instant::now();
            op_a();
            let a = start.elapsed();
            let start = Instant::now();
            op_b();
            let b = start.elapsed();
            best_a = best_a.min(a);
            best_b = best_b.min(b);
            b.as_secs_f64() / a.as_secs_f64().max(f64::MIN_POSITIVE)
        })
        .collect();
    ratios.sort_by(f64::total_cmp);
    (ratios[ratios.len() / 2], best_a, best_b)
}

fn assert_budget(name: &str, measured: Duration, budget_ms: u64) {
    eprintln!(
        "[perf-budget] {name}: {:.1} ms (budget {budget_ms} ms)",
        measured.as_secs_f64() * 1e3
    );
    assert!(
        measured <= Duration::from_millis(budget_ms),
        "{name} took {:.1} ms, budget {budget_ms} ms",
        measured.as_secs_f64() * 1e3
    );
}

#[test]
fn plane_only_cube_shell_within_budget() {
    if debug_skip("plane-only cube shell") {
        return;
    }
    let cube = popkern_core::topology::BRep::unit_cube(tol()).expect("cube");
    let top = cube.faces.keys().next().expect("cube has faces");
    let measured = best_of(5, || {
        let r = popkern_shell::shell(
            &cube,
            0.05,
            ShellDirection::Inward,
            &[top],
            &[],
            &ShellOptions::default(),
            tol(),
            None,
        )
        .expect("cube shell serves");
        assert!(r.brep.faces.len() > 6, "shell must add walls and skin");
    });
    // Measured 1.3 ms best-of-5 release (2026-08-27) × 5 ≈ 7 ms.
    assert_budget("plane-only cube shell", measured, 7);
}

#[test]
fn plane_dominant_thicken_within_2x_of_plane_only() {
    // Design §8: plane-dominant inputs stay within 2× of the pure
    // plane pipeline. Ratio gate (machine-speed independent): the
    // l_bracket adds concave joins and a bigger face count over the
    // cube's pure plane ladder; measured single-thread median pair
    // ratio 1.62 idle, 1.32..1.81 with every core busy, and 3.23 with
    // the l_bracket side's work doubled (2026-09-24, PK-526).
    if debug_skip("plane-dominant ratio") {
        return;
    }
    let cube = popkern_core::topology::BRep::unit_cube(tol()).expect("cube");
    let lb = l_bracket(tol()).expect("fixture builds");
    // Both sides run in ONE dedicated single-thread pool: the gate
    // bounds algorithmic cost, and thicken's rayon fan-out otherwise
    // loses a different share of its speedup per fixture under machine
    // load (measured ratio 1.28..2.70 with every core busy, PK-526).
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .expect("single-thread pool");
    let (ratio, cube_time, lb_time) = pool.install(|| {
        median_pair_ratio(
            21,
            || {
                let r = thicken(&cube, 0.02, &ShellOptions::default(), tol()).expect("serves");
                // Minkowski-rounded outer skin: 6 planes each side, 12 arc
                // strips, 8 octant corners.
                assert_eq!(r.brep.faces.len(), 32, "wall body");
            },
            || {
                let r = thicken(&lb.brep, 0.02, &ShellOptions::default(), tol()).expect("serves");
                assert!(r.brep.validate().is_ok());
            },
        )
    });
    let (cube_ms, lb_ms) = (cube_time.as_secs_f64(), lb_time.as_secs_f64());
    eprintln!("[perf-budget] plane-dominant ratio: {ratio:.2} (l_bracket/cube thicken)");
    assert!(
        ratio <= 2.0,
        "plane-dominant thicken ratio {ratio:.2} exceeds the §8 2× target \
         (best l_bracket {:.2} ms vs cube {:.2} ms)",
        lb_ms * 1e3,
        cube_ms * 1e3
    );
}

#[test]
fn filleted_cube_thicken_within_budget() {
    if debug_skip("filleted-cube thicken") {
        return;
    }
    let f = filleted_cube(tol()).expect("fixture builds");
    let measured = best_of(3, || {
        let r = thicken(&f.brep, 0.02, &ShellOptions::default(), tol()).expect("serves");
        assert_eq!(r.brep.faces.len(), 2 * f.expected_faces, "wall body");
    });
    // Measured 57 ms best-of-3 release (2026-08-27) × 5 ≈ 285 ms.
    assert_budget("filleted-cube thicken", measured, 285);
}

#[test]
fn wavy_plate_thicken_within_budget() {
    if debug_skip("wavy-plate thicken") {
        return;
    }
    let f = wavy_plate(tol()).expect("fixture builds");
    let measured = best_of(3, || {
        let r = thicken(&f.brep, 0.02, &ShellOptions::default(), tol()).expect("serves");
        assert!(r.brep.validate().is_ok());
        assert!(
            r.warnings.iter().any(|w| matches!(
                w,
                popkern_shell::ShellWarning::FreeFormFoldScan {
                    truncated: false,
                    ..
                }
            )),
            "tier-3 scan must run clean on the NURBS plate"
        );
    });
    // Measured 170 ms best-of-3 release (2026-08-27; was 5405 ms
    // pre-campaign — banded fitting + parallel edge recompute) × 5 ≈
    // 850 ms.
    assert_budget("wavy-plate thicken", measured, 850);
}

#[test]
fn sphere_capped_rod_thicken_within_budget() {
    if debug_skip("sphere-capped-rod thicken") {
        return;
    }
    let f = sphere_capped_rod(tol()).expect("fixture builds");
    let measured = best_of(5, || {
        let r = thicken(&f.brep, 0.02, &ShellOptions::default(), tol()).expect("serves");
        assert!(r.brep.validate().is_ok());
    });
    // Measured 7 ms best-of-5 release (2026-08-27; was 256 ms — the
    // concentric-cap pair burned the subdivision seeder's budget until
    // the analytic separation certificate) × 5 ≈ 35 ms.
    assert_budget("sphere-capped-rod thicken", measured, 35);
}

#[test]
fn max_safe_offset_corpus_within_budget() {
    if debug_skip("max_safe_offset corpus sweep") {
        return;
    }
    let fixtures = corpus(tol()).expect("corpus builds");
    let measured = best_of(3, || {
        for f in &fixtures {
            for dir in [ShellDirection::Inward, ShellDirection::Outward] {
                let rep = max_safe_offset(&f.brep, dir, tol()).expect("report");
                assert!(rep.max_distance > 0.0, "{}: degenerate bound", f.name);
            }
        }
    });
    // Whole-corpus sweep (10 bodies × 2 directions) measured 45 ms
    // best-of-3 release (2026-08-27) × 5 ≈ 225 ms. The per-query
    // <100 ms interactive gate lives in max_safe_offset_gates.
    assert_budget("max_safe_offset corpus sweep", measured, 225);
}
