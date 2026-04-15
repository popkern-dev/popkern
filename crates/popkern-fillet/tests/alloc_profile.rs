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
//! | Cone base (walker path)       | 176     | 30 KiB   | Pure walker, no overlap pair|
//! | Single-edge cube (legacy)     | 3,047   | 322 KiB  | No overlap pairs (1 edge)   |
//! | 12-edge cube (legacy, batch)  | 243,802 | 71.7 MiB | 84% in DetectOverlaps phase |
//!
//! **Per-phase breakdown for the 12-edge cube** (same workload):
//!
//! | Phase            | Allocs  | Bytes    | % of total |
//! |------------------|---------|----------|------------|
//! | Validate         |      24 |     624  |   0.01%    |
//! | Analyze          |      85 |  14,796  |   0.03%    |
//! | ComputeGeometry  |   1,167 |   192 KB |   0.47%    |
//! | DetectOverlaps   | 205,484 |  66.7 MB |  84.3%     |
//! | ApplyTopology    |   9,037 |   743 KB |   3.7%     |
//! | RepairPcurves    |  28,000 |   4.0 MB |  11.5%     |
//!
//! **Key finding:** DetectOverlaps dominates (84% of allocs, 93% of
//! bytes). Round 1 of optimisation — an AABB pre-filter in
//! `overlap::curves_intersect` — already dropped DetectOverlaps from
//! 303,788 → 205,484 allocs (−32%) by skipping the full NURBS
//! intersector when two trim curves' bounding boxes don't overlap.
//! The remaining 205k allocs come from `intersect_curves` ⇒
//! `intersect_nurbs_nurbs` ⇒ `decompose_to_bezier` +
//! bezier clipping; further gains require cross-crate work on
//! `NurbsCurve::de_boor` (which allocates a fresh homogeneous Vec
//! on every point evaluation) or a tighter overlap heuristic.
//!
//! ## What this commit ships vs. defers
//!
//! **Ships:**
//! - Measurement infrastructure (`CountingAlloc` + `profile`)
//! - Baseline numbers for three representative operations
//! - Per-phase breakdown via `PhaseAllocProfiler`
//! - AABB pre-filter in `overlap::curves_intersect` (−98k allocs)
//! - Tightened regression guards around the new baselines
//!
//! **Defers:**
//! - Cross-crate `de_boor` refactor to eliminate per-call Vec alloc
//! - Arena/pool adoption in `decompose_to_bezier` and the bezier
//!   clipping recursion (popkern-core scope)
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

    // Upper bound set ~40% above the current ~244k alloc baseline
    // to absorb normal variance while catching unbounded regressions
    // (e.g. an accidental O(n²) push slipping into a hot loop). If
    // this trips, first check the per-phase breakdown to see which
    // phase regressed rather than assuming the fix is "bump the
    // threshold".
    assert!(
        delta.count < 350_000,
        "§18 regression: 12-edge cube fillet made {} allocations \
         (threshold 350,000; baseline ~244k) — run the per-phase \
         profiler to locate the regression",
        delta.count
    );
    assert!(
        delta.bytes < 110 * 1024 * 1024,
        "§18 regression: 12-edge cube fillet allocated {} bytes \
         (threshold 110 MiB; baseline ~72 MB) — check for \
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
    // (~3k allocs) catches per-edge allocation explosion without
    // false-positiving on variance.
    assert!(
        delta.count < 15_000,
        "§18 regression: single-edge cube fillet made {} allocations \
         (threshold 15,000; baseline ~3k) — check for regression on \
         the fast path",
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
    // per-Newton-step temporaries. Baseline ~176 allocs; threshold
    // set to 10× headroom to absorb any future change in sampling
    // or Newton iteration count without false-positiving.
    assert!(
        delta.count < 2_000,
        "§18 regression: cone walker fillet made {} allocations \
         (threshold 2,000; baseline ~176) — check for regression \
         on the walker path",
        delta.count
    );
}
