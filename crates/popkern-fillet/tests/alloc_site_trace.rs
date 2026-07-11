//! Path-D diagnostic: locate the downstream cascade amplifier
//! triggered by authoring linear pcurves on both twins of each
//! corner-split MEF edge (Stage C1+, commit `eb2f672d`).
//!
//! Approach: install a `BacktraceCountingAlloc` as the
//! `#[global_allocator]` for this binary. Sample every `SAMPLE_RATE`
//! allocation and record an unresolved backtrace (instruction
//! pointers only — cheap). At the end, resolve symbols once and
//! print the top-N stack signatures by sample count.
//!
//! Two runs:
//! 1. `cascade_baseline_12_edge_cube` — full Stage C1+ as-of-main:
//!    both twins of each MEF corner-split edge have linear pcurves.
//!    Expected: ~612k allocs, ~6.1k samples at rate 100.
//! 2. `cascade_inhibited_12_edge_cube` — same flow but with the
//!    parent-side half-edge pcurve write deleted (via the env
//!    `POPKERN_FILLET_INHIBIT_PARENT_PCURVE=1` hook the prober
//!    installs). Expected: ~310k allocs, ~3.1k samples.
//!
//! Diff: the backtraces over-represented in run 1 vs run 2 identify
//! the downstream pass that amplifies the second pcurve into the
//! cascade.

#![cfg(test)]

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::collections::HashMap;
use std::sync::Mutex;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

/// Sample every Nth allocation. 100 keeps the run time bounded;
/// backtrace capture is ~1µs per call, so 6000 samples adds ~6ms
/// to a 12-edge cube fillet that already runs in ~30ms.
const SAMPLE_RATE: usize = 100;

/// Skip the top N frames of every backtrace — the allocator hook
/// itself and the `System` allocator passthrough are stable noise.
const SKIP_TOP_FRAMES: usize = 4;

/// Keep only this many leading frames after the skip. Deeper frames
/// are usually framework noise; the immediate callers are where
/// per-pcurve work lives.
const KEEP_FRAMES: usize = 12;

thread_local! {
    /// Re-entry guard. `backtrace::trace` allocates internally;
    /// without this guard, every capture would re-enter the
    /// allocator and either deadlock or spiral into unbounded
    /// recursion.
    static IN_CAPTURE: Cell<bool> = const { Cell::new(false) };
}

/// Sampled allocation sites: backtrace hash → (sample count, frame
/// instruction pointers).
type AllocSiteSamples = HashMap<u64, (usize, Vec<usize>)>;

static ALLOC_COUNT: AtomicUsize = AtomicUsize::new(0);
static SAMPLES: Mutex<Option<AllocSiteSamples>> = Mutex::new(None);
static TRACING_ENABLED: AtomicBool = AtomicBool::new(false);

struct BacktraceCountingAlloc;

unsafe impl GlobalAlloc for BacktraceCountingAlloc {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let n = ALLOC_COUNT.fetch_add(1, Ordering::Relaxed);
        if TRACING_ENABLED.load(Ordering::Relaxed) && n.is_multiple_of(SAMPLE_RATE) {
            capture_sample();
        }
        unsafe { System.alloc(layout) }
    }
    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let n = ALLOC_COUNT.fetch_add(1, Ordering::Relaxed);
        if TRACING_ENABLED.load(Ordering::Relaxed) && n.is_multiple_of(SAMPLE_RATE) {
            capture_sample();
        }
        unsafe { System.alloc_zeroed(layout) }
    }
    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let n = ALLOC_COUNT.fetch_add(1, Ordering::Relaxed);
        if TRACING_ENABLED.load(Ordering::Relaxed) && n.is_multiple_of(SAMPLE_RATE) {
            capture_sample();
        }
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

#[global_allocator]
static ALLOC: BacktraceCountingAlloc = BacktraceCountingAlloc;

fn capture_sample() {
    let already = IN_CAPTURE.with(|f| f.replace(true));
    if already {
        return;
    }
    let mut frames: Vec<usize> = Vec::with_capacity(KEEP_FRAMES);
    let mut seen = 0usize;
    backtrace::trace(|frame| {
        if seen < SKIP_TOP_FRAMES {
            seen += 1;
            return true;
        }
        if frames.len() == KEEP_FRAMES {
            return false;
        }
        frames.push(frame.ip() as usize);
        true
    });
    let mut hash: u64 = 1469598103934665603u64;
    for &ip in &frames {
        hash ^= ip as u64;
        hash = hash.wrapping_mul(1099511628211);
    }
    let mut guard = SAMPLES.lock().unwrap();
    let map = guard.get_or_insert_with(HashMap::new);
    map.entry(hash).or_insert_with(|| (0, frames.clone())).0 += 1;
    drop(guard);
    IN_CAPTURE.with(|f| f.set(false));
}

fn reset_tracing() {
    ALLOC_COUNT.store(0, Ordering::Relaxed);
    *SAMPLES.lock().unwrap() = Some(HashMap::new());
}

fn enable_tracing() {
    TRACING_ENABLED.store(true, Ordering::Relaxed);
}

fn disable_tracing() {
    TRACING_ENABLED.store(false, Ordering::Relaxed);
}

fn finish_and_dump(label: &str, top_n: usize) -> Vec<(u64, usize, Vec<String>)> {
    let map = SAMPLES.lock().unwrap().take().unwrap_or_default();
    let mut entries: Vec<(u64, usize, Vec<usize>)> = map
        .into_iter()
        .map(|(h, (count, frames))| (h, count, frames))
        .collect();
    entries.sort_by_key(|entry| std::cmp::Reverse(entry.1));
    entries.truncate(top_n);
    let total_samples: usize = entries.iter().map(|(_, c, _)| *c).sum();
    eprintln!(
        "\n## {label} — top {} backtrace groups (out of {} sampled allocations, rate 1/{})\n",
        entries.len(),
        total_samples,
        SAMPLE_RATE
    );
    let mut resolved = Vec::with_capacity(entries.len());
    for (hash, count, frames) in entries {
        let symbols: Vec<String> = frames
            .iter()
            .map(|&ip| {
                let mut name = format!("0x{ip:x}");
                backtrace::resolve(ip as *mut _, |symbol| {
                    if let Some(symname) = symbol.name() {
                        name = format!("{symname}");
                        if let (Some(filename), Some(line)) = (symbol.filename(), symbol.lineno()) {
                            name = format!(
                                "{symname} ({}:{line})",
                                filename.file_name().and_then(|s| s.to_str()).unwrap_or("?")
                            );
                        }
                    }
                });
                name
            })
            .collect();
        eprintln!("--- hash={hash:016x} count={count} ---");
        for (i, sym) in symbols.iter().enumerate() {
            eprintln!("  [{i:2}] {sym}");
        }
        resolved.push((hash, count, symbols));
    }
    resolved
}

fn run_cube_fillet() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    // Warmup outside measurement.
    let _ = fillet_edges(
        &cube,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    );
    reset_tracing();
    enable_tracing();
    let _ = fillet_edges(
        &cube,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    );
    disable_tracing();
}

#[test]
#[ignore = "diagnostic — run with --ignored --nocapture"]
fn cascade_baseline_12_edge_cube_backtrace_histogram() {
    run_cube_fillet();
    let total = ALLOC_COUNT.load(Ordering::Relaxed);
    eprintln!("\nTotal allocations recorded: {total}");
    let _ = finish_and_dump("12-edge cube fillet (Stage C1+ both-twin pcurves)", 30);
}
