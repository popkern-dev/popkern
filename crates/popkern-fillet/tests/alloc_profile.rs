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
//! ## Baseline numbers (2026-04-15, debug build, local dev box)
//!
//! | Operation                     | Allocs  | Bytes    | Per-edge |
//! |-------------------------------|---------|----------|----------|
//! | Cone base (walker path)       | 621     | 99 KiB   | 621      |
//! | Single-edge cube (legacy)     | 18,025  | 5.2 MiB  | 18,025   |
//! | 12-edge cube (legacy, batch)  | 342,107 | 103.6 MiB| ~28,500  |
//!
//! **Key finding:** the phase-1b walker path is ~30× more
//! allocation-efficient per edge than the legacy offset-intersect
//! path. Future arena/pool work should target the legacy path
//! (`offset::analytic` / `center_curve::constant_radius_center_curve`
//! and the core `intersect_surfaces` they call into), NOT the
//! walker. Targeting the walker would optimise 0.2% of the
//! allocation volume; targeting the legacy path would hit 99%.
//!
//! ## What this commit ships vs. defers
//!
//! **Ships:**
//! - Measurement infrastructure (`CountingAlloc` + `profile`)
//! - Baseline numbers for three representative operations
//! - Hard regression guards with generous upper bounds
//!
//! **Defers:**
//! - Actual arena/pool adoption. The baselines above identify the
//!   hot paths; a follow-up commit can wire a `bumpalo::Bump` or
//!   `typed-arena::Arena` through the legacy offset-intersect
//!   pipeline. Cross-crate work: most of the allocation volume
//!   originates in popkern-core's `intersect_surfaces` and NURBS
//!   fitter, so arena adoption needs coordinated changes there.
//!
//! Run with:
//! ```sh
//! cargo test -p popkern-fillet --test alloc_profile -- --nocapture
//! ```

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};

static ALLOC_COUNT: AtomicUsize = AtomicUsize::new(0);
static ALLOC_BYTES: AtomicUsize = AtomicUsize::new(0);
static DEALLOC_COUNT: AtomicUsize = AtomicUsize::new(0);

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
// Baseline measurements
// ---------------------------------------------------------------------------

#[test]
fn baseline_12_edge_cube_fillet_allocation_profile() {
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

    // Generous upper bound: a 12-edge cube fillet producing a
    // 26-face result should allocate on the order of tens of
    // thousands of items, not millions. The absolute threshold
    // catches unbounded-allocation regressions (e.g. an O(n²)
    // vec push slipping into a hot loop) without false-positiving
    // on normal variance.
    assert!(
        delta.count < 1_000_000,
        "§18 regression: 12-edge cube fillet made {} allocations \
         (threshold 1,000,000) — check for unbounded allocation growth",
        delta.count
    );
    assert!(
        delta.bytes < 200 * 1024 * 1024,
        "§18 regression: 12-edge cube fillet allocated {} bytes \
         (threshold 200 MiB) — check for unbounded buffer growth",
        delta.bytes
    );
}

#[test]
fn baseline_single_edge_cube_fillet_allocation_profile() {
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

    // Single-edge should be dramatically smaller than the 12-edge
    // case. Upper bound set to catch per-edge allocation explosion.
    assert!(
        delta.count < 200_000,
        "§18 regression: single-edge cube fillet made {} allocations \
         (threshold 200,000) — check for regression on the fast path",
        delta.count
    );
}

#[test]
fn baseline_cone_base_fillet_allocation_profile() {
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
    // per-Newton-step temporaries. Upper bound set generously to
    // catch regressions without false-positiving on the current
    // baseline, which is unknown until the first profiling run.
    assert!(
        delta.count < 2_000_000,
        "§18 regression: cone walker fillet made {} allocations \
         (threshold 2,000,000) — check for regression on the walker path",
        delta.count
    );
}
