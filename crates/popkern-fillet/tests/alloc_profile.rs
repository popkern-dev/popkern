//! §18 memory-high-water / allocation-profile measurement.
//!
//! Integration test that installs a `CountingAlloc` as the
//! `#[global_allocator]` for this binary only, runs representative
//! fillet operations, and prints per-operation allocation counts
//! and bytes. Doubles as a regression guard: a hard assertion on
//! a generous upper bound catches unbounded-allocation regressions
//! (e.g. an accidental `O(n²)` vec push in a hot loop) without
//! requiring external tools like Valgrind's DHAT or heaptrack.
//!
//! ## Why a custom allocator instead of `dhat` / `heaptrack`
//!
//! - **No external toolchain dependency.** Counting atomic adds
//!   on every alloc/dealloc is cheap (≤ 10ns per op) and runs on
//!   every target that supports Rust's stable `GlobalAlloc` trait.
//! - **Queryable at runtime.** §18 requires telemetry to be
//!   "queryable" — a static atomic counter we read before and
//!   after a pipeline call IS queryable, no profiling UI required.
//! - **Per-binary scope.** `#[global_allocator]` inside an
//!   integration test only affects that test's binary. No risk
//!   of coupling allocation measurement to the library's runtime
//!   behaviour.
//!
//! ## What this tests ships
//!
//! 1. **Baseline measurement** for three representative fillet
//!    operations, printed to stderr via `--nocapture`.
//! 2. **Hard upper bounds** on alloc count and byte total,
//!    generous enough to absorb normal variance but tight enough
//!    to catch regressions like an `O(n²)` explosion.
//!
//! ## Baseline numbers (2026-04-14, debug build, single-threaded)
//!
//! | Operation                     | Allocs  | Bytes    | Notes                       |
//! |-------------------------------|---------|----------|-----------------------------|
//! | Cone base (walker path)       |      56 |   7.8 KB | Pure walker, no overlap pair|
//! | Single-edge cube (legacy)     |   2,131 |    251 KB| No overlap pairs (1 edge)   |
//! | 12-edge cube (legacy, batch)  |  40,943 |   6.5 MB | Post-de_boor + AABB filter  |
//!
//! **Per-phase breakdown for the 12-edge cube** (same workload):
//!
//! | Phase            | Allocs  | Bytes    | % of total |
//! |------------------|---------|----------|------------|
//! | Validate         |      24 |     624  |   0.06%    |
//! | Analyze          |      85 |  14,796  |   0.21%    |
//! | ComputeGeometry  |     639 |   160 KB |   1.6%     |
//! | DetectOverlaps   |  13,724 |   2.4 MB |  33.5%     |
//! | ApplyTopology    |   6,720 |   537 KB |  16.4%     |
//! | RepairPcurves    |  19,744 |   3.2 MB |  48.2%     |
//!
//! **Optimisation history:**
//!
//! 1. **Baseline (commit 5b136e1):** 342,107 allocs / 103 MB. Flat
//!    profile was deceptive — investigation showed DetectOverlaps
//!    was 89% of the volume.
//!
//! 2. **AABB pre-filter (commit dcc0f4a):** 243,802 allocs / 75 MB
//!    (−29% / −31%). `overlap::curves_intersect` now checks curve
//!    bounding boxes before wrapping curves in `CurveGeometry::Nurbs`
//!    and calling the full NURBS intersector. Dropped DetectOverlaps
//!    from 303,788 → 205,484 allocs (−32%).
//!
//! 3. **Allocation-free de Boor (this commit):** 40,943 allocs / 6.5 MB
//!    (−83% / −91% from AABB-filtered baseline). `NurbsCurve::de_boor`
//!    used to build a fresh `Vec<Homogeneous>` on every `point_at` call
//!    and `de_boor_1d` kept a Vec scratch buffer for the recursion.
//!    Both now use a fixed-size stack array (max degree 31). Total
//!    improvement from 5b136e1 baseline: **−88% allocs, −94% bytes**.
//!
//! The dominant phase is now RepairPcurves (pcurve recomputation),
//! which is inherently allocation-heavy because it builds new NURBS
//! pcurves per half-edge. Further reductions would require changing
//! how pcurves are stored/computed — significantly larger scope.
//!
//! Run with:
//! ```sh
//! cargo test -p popkern-fillet --test alloc_profile -- --nocapture
//! ```

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    FilletOptions, FilletRadius, FilletSpec, ProgressEvent, ProgressReporter, fillet_edges,
    fillet_edges_with_progress,
};
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

static ALLOC_COUNT: AtomicUsize = AtomicUsize::new(0);
static ALLOC_BYTES: AtomicUsize = AtomicUsize::new(0);
static DEALLOC_COUNT: AtomicUsize = AtomicUsize::new(0);

/// Serialises the profiled tests so parallel `cargo test` execution
/// doesn't blend alloc counts from concurrent tests into one another.
/// Every test that calls `profile()` or installs a `PhaseAllocProfiler`
/// must `let _guard = PROFILE_LOCK.lock();` at the top.
static PROFILE_LOCK: Mutex<()> = Mutex::new(());

struct CountingAlloc;

unsafe impl GlobalAlloc for CountingAlloc {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ALLOC_COUNT.fetch_add(1, Ordering::Relaxed);
        ALLOC_BYTES.fetch_add(layout.size(), Ordering::Relaxed);
        unsafe { System.alloc(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        DEALLOC_COUNT.fetch_add(1, Ordering::Relaxed);
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        ALLOC_COUNT.fetch_add(1, Ordering::Relaxed);
        ALLOC_BYTES.fetch_add(layout.size(), Ordering::Relaxed);
        unsafe { System.alloc_zeroed(layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        // Realloc counts as a new allocation for our purposes —
        // it either shrinks in place (rare) or moves + frees the
        // old buffer. Either way the size delta is worth tracking.
        ALLOC_COUNT.fetch_add(1, Ordering::Relaxed);
        if new_size > layout.size() {
            ALLOC_BYTES.fetch_add(new_size - layout.size(), Ordering::Relaxed);
        }
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOC: CountingAlloc = CountingAlloc;

/// Snapshot of the allocation counters at a point in time.
#[derive(Debug, Clone, Copy)]
struct AllocSnapshot {
    count: usize,
    bytes: usize,
    deallocs: usize,
}

impl AllocSnapshot {
    fn capture() -> Self {
        Self {
            count: ALLOC_COUNT.load(Ordering::Relaxed),
            bytes: ALLOC_BYTES.load(Ordering::Relaxed),
            deallocs: DEALLOC_COUNT.load(Ordering::Relaxed),
        }
    }
    fn delta(self, start: Self) -> Self {
        Self {
            count: self.count.saturating_sub(start.count),
            bytes: self.bytes.saturating_sub(start.bytes),
            deallocs: self.deallocs.saturating_sub(start.deallocs),
        }
    }
}

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Run `body`, report the allocation delta it produced, and
/// return the delta for further assertions.
fn profile<F: FnOnce()>(label: &str, body: F) -> AllocSnapshot {
    let start = AllocSnapshot::capture();
    body();
    let end = AllocSnapshot::capture();
    let delta = end.delta(start);
    eprintln!(
        "§18 {label}: {} allocs, {} bytes, {} deallocs (net {} live)",
        delta.count,
        delta.bytes,
        delta.deallocs,
        delta.count.saturating_sub(delta.deallocs)
    );
    delta
}

// ---------------------------------------------------------------------------
// Per-phase allocation profiler
// ---------------------------------------------------------------------------

/// `ProgressReporter` impl that snapshots the alloc counters at
/// every `on_phase` event and prints the delta between consecutive
/// phases, giving a per-phase allocation breakdown.
///
/// Important caveat: phases are contiguous in wall time but
/// ALSO include nested calls — e.g. `ComputeGeometry` covers the
/// whole per-edge loop including `offset_surface`,
/// `intersect_surfaces`, `build_fillet_surface`, and trim-curve
/// fitting. The per-phase numbers tell us which phase of the
/// pipeline dominates, not which inner function.
struct PhaseAllocProfiler {
    last_snapshot: AllocSnapshot,
    first: bool,
    label: &'static str,
}

impl PhaseAllocProfiler {
    fn new(label: &'static str) -> Self {
        Self {
            last_snapshot: AllocSnapshot::capture(),
            first: true,
            label,
        }
    }
}

impl ProgressReporter for PhaseAllocProfiler {
    fn on_phase(&mut self, event: ProgressEvent) {
        let now = AllocSnapshot::capture();
        let delta = now.delta(self.last_snapshot);
        if self.first {
            eprintln!("§18 [{label}] per-phase breakdown:", label = self.label);
            self.first = false;
        }
        eprintln!(
            "  {phase:?}: {allocs} allocs, {bytes} bytes",
            phase = event.phase,
            allocs = delta.count,
            bytes = delta.bytes
        );
        self.last_snapshot = now;
    }
}

#[test]
fn per_phase_allocation_breakdown_12_edge_cube() {
    let _guard = PROFILE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    // Warmup the pipeline outside the measurement window.
    let cube = BRep::unit_cube(tol()).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let _ = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None);

    // Measured run with a per-phase profiler reporter.
    let mut profiler = PhaseAllocProfiler::new("12-edge cube");
    let _ = fillet_edges_with_progress(
        &cube,
        &specs,
        &FilletOptions::default(),
        tol(),
        None,
        &mut profiler,
    );
    // No assertion — this is a diagnostic printer. The next
    // commit will use this output to identify which phase needs
    // arena adoption first.
}

#[test]
fn per_phase_allocation_breakdown_cone_walker() {
    let _guard = PROFILE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let cone = BRep::cone(1.0, 2.0, tol()).expect("cone BRep should build");
    let chosen = cone
        .edges
        .keys()
        .find(|&e| {
            let he = cone.edges[e].half_edge;
            let twin = cone.half_edges[he].twin;
            let face_a = cone.loops[cone.half_edges[he].loop_].face;
            let face_b = cone.loops[cone.half_edges[twin].loop_].face;
            face_a != face_b
        })
        .unwrap_or_else(|| cone.edges.keys().next().unwrap());
    let specs = vec![FilletSpec {
        edge: chosen,
        radius: FilletRadius::Constant(0.005),
    }];
    let _ = fillet_edges(&cone, &specs, &FilletOptions::default(), tol(), None);

    let mut profiler = PhaseAllocProfiler::new("cone base walker");
    let _ = fillet_edges_with_progress(
        &cone,
        &specs,
        &FilletOptions::default(),
        tol(),
        None,
        &mut profiler,
    );
}

// ---------------------------------------------------------------------------
// Baseline measurements (gross totals, no phase breakdown)
// ---------------------------------------------------------------------------

#[test]
fn baseline_12_edge_cube_fillet_allocation_profile() {
    let _guard = PROFILE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    // Setup: build the cube and specs outside the profiled block
    // so setup costs don't distort the pipeline's own allocation
    // numbers. We only profile the fillet call itself.
    let cube = BRep::unit_cube(tol()).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();

    // Warmup: run once to pull caches/JIT into steady state.
    let _ = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None);

    let delta = profile("12-edge cube fillet", || {
        let _ = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None);
    });

    // Upper bound set ~50% above the current ~41k alloc baseline
    // to absorb normal variance while catching unbounded regressions
    // (e.g. an accidental O(n²) push slipping into a hot loop). If
    // this trips, first check the per-phase breakdown to see which
    // phase regressed rather than assuming the fix is "bump the
    // threshold".
    assert!(
        delta.count < 60_000,
        "§18 regression: 12-edge cube fillet made {} allocations \
         (threshold 60,000; baseline ~41k) — run the per-phase \
         profiler to locate the regression",
        delta.count
    );
    assert!(
        delta.bytes < 10 * 1024 * 1024,
        "§18 regression: 12-edge cube fillet allocated {} bytes \
         (threshold 10 MiB; baseline ~6.5 MB) — check for \
         unbounded buffer growth",
        delta.bytes
    );
}

#[test]
fn baseline_single_edge_cube_fillet_allocation_profile() {
    let _guard = PROFILE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];

    // Warmup.
    let _ = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None);

    let delta = profile("single-edge cube fillet", || {
        let _ = fillet_edges(&cube, &specs, &FilletOptions::default(), tol(), None);
    });

    // Single-edge is a single spec with no overlap pairs, so should
    // be tiny relative to the 12-edge case. Threshold ~5× baseline
    // (~2.1k allocs) catches per-edge allocation explosion without
    // false-positiving on variance.
    assert!(
        delta.count < 10_000,
        "§18 regression: single-edge cube fillet made {} allocations \
         (threshold 10,000; baseline ~2.1k) — check for regression \
         on the fast path",
        delta.count
    );
}

#[test]
fn baseline_cone_base_fillet_allocation_profile() {
    let _guard = PROFILE_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    // The cone base fillet takes the phase-1b walker path, which
    // is expected to allocate more than the legacy offset-intersect
    // path because of per-sample contact solutions. This baseline
    // captures the current state so future walker optimisations
    // (arena pooling, Vec reuse) can be measured against it.
    let cone = BRep::cone(1.0, 2.0, tol()).expect("cone BRep should build");
    let chosen = cone
        .edges
        .keys()
        .find(|&e| {
            let he = cone.edges[e].half_edge;
            let twin = cone.half_edges[he].twin;
            let face_a = cone.loops[cone.half_edges[he].loop_].face;
            let face_b = cone.loops[cone.half_edges[twin].loop_].face;
            face_a != face_b
        })
        .unwrap_or_else(|| cone.edges.keys().next().unwrap());
    let specs = vec![FilletSpec {
        edge: chosen,
        radius: FilletRadius::Constant(0.005),
    }];

    // Warmup.
    let _ = fillet_edges(&cone, &specs, &FilletOptions::default(), tol(), None);

    let delta = profile("cone base fillet (walker path)", || {
        let _ = fillet_edges(&cone, &specs, &FilletOptions::default(), tol(), None);
    });

    // Walker path allocates per sample (default 20 samples) plus
    // per-Newton-step temporaries. Baseline ~56 allocs; threshold
    // set to 10× headroom to absorb any future change in sampling
    // or Newton iteration count without false-positiving.
    assert!(
        delta.count < 600,
        "§18 regression: cone walker fillet made {} allocations \
         (threshold 600; baseline ~56) — check for regression \
         on the walker path",
        delta.count
    );
}
