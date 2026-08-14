//! Content-keyed memoization for the Gregory self-intersection
//! report (PK-273).
//!
//! `gregory_self_intersection_report` is a pure function of the patch
//! content and the tolerance, and its slow tier (bidegree-9 NURBS
//! conversion + two sub-patch surface-surface intersections) was
//! measured at ~25 ms per patch — the binding term of whole-model
//! fast validation on Gregory-laden bodies (96 patches ≈ 2.4 s per
//! sweep). Construction fold gates and the LM polish re-check
//! unchanged patches repeatedly; memoizing by content turns every
//! repeat into a hash + lookup.
//!
//! Design (docs/architecture/region-scoped-validation.md):
//! - **Key**: patch serde byte stream hashed by two independently
//!   seeded 64-bit hashers, plus the byte length and the tolerance
//!   bits. Content identity is exact (bit patterns, so NaN/-0.0 are
//!   distinct keys, never value-compared); collision probability is
//!   ~2⁻¹²⁸ per pair — far below geometric noise.
//! - **Storage**: thread-local, so no locks and no cross-thread
//!   coherence; each thread warms its own cache.
//! - **Bound**: the map is cleared when it exceeds
//!   [`CACHE_CAPACITY`] entries. Clearing (not evicting) keeps the
//!   policy deterministic; a re-miss only re-pays the old cost.
//! - **No invalidation**: keys are content-derived, so a mutated
//!   patch is simply a different key.

use std::cell::RefCell;
use std::collections::HashMap;
use std::hash::Hasher;

use super::GregorySelfIntersectionReport;
use crate::geometry::gregory::GregoryPatch;

/// Entry bound for the per-thread verdict map. Keys are ~40 bytes and
/// reports ~40 bytes, so the worst case is ~1 MB per thread.
const CACHE_CAPACITY: usize = 8192;

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct CacheKey {
    content_len: u64,
    content_hash_a: u64,
    content_hash_b: u64,
    tolerance_bits: u64,
}

#[derive(Default)]
struct CacheState {
    map: HashMap<CacheKey, GregorySelfIntersectionReport>,
    hits: u64,
    misses: u64,
}

thread_local! {
    static CACHE: RefCell<CacheState> = RefCell::new(CacheState::default());
}

/// Feeds one byte stream into two independently seeded hashers while
/// counting length, so a single serialization pass yields the whole
/// key.
struct DualHashWriter {
    hash_a: std::hash::DefaultHasher,
    hash_b: std::hash::DefaultHasher,
    len: u64,
}

impl DualHashWriter {
    fn new() -> Self {
        let mut hash_a = std::hash::DefaultHasher::new();
        let mut hash_b = std::hash::DefaultHasher::new();
        // Distinct seed prefixes make the two 64-bit digests
        // effectively independent for non-adversarial input.
        hash_a.write_u64(0x9e37_79b9_7f4a_7c15);
        hash_b.write_u64(0x2545_f491_4f6c_dd1d);
        Self {
            hash_a,
            hash_b,
            len: 0,
        }
    }
}

impl std::io::Write for DualHashWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.hash_a.write(buf);
        self.hash_b.write(buf);
        self.len += buf.len() as u64;
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

fn content_key(patch: &GregoryPatch, tolerance: f64) -> Option<CacheKey> {
    let mut writer = DualHashWriter::new();
    serde_json::to_writer(&mut writer, patch).ok()?;
    Some(CacheKey {
        content_len: writer.len,
        content_hash_a: writer.hash_a.finish(),
        content_hash_b: writer.hash_b.finish(),
        tolerance_bits: tolerance.to_bits(),
    })
}

/// Memoized front for the self-intersection report. Falls back to a
/// direct computation if the patch cannot be serialized (no such
/// variant exists today; the fallback keeps the cache strictly an
/// optimization).
pub(super) fn cached_report(
    patch: &GregoryPatch,
    tolerance: f64,
    compute: impl FnOnce() -> GregorySelfIntersectionReport,
) -> GregorySelfIntersectionReport {
    let Some(key) = content_key(patch, tolerance) else {
        return compute();
    };
    if let Some(report) = CACHE.with(|cache| {
        let mut cache = cache.borrow_mut();
        let found = cache.map.get(&key).cloned();
        match found {
            Some(report) => {
                cache.hits += 1;
                Some(report)
            }
            None => {
                cache.misses += 1;
                None
            }
        }
    }) {
        return report;
    }
    let report = compute();
    CACHE.with(|cache| {
        let mut cache = cache.borrow_mut();
        if cache.map.len() >= CACHE_CAPACITY {
            cache.map.clear();
        }
        cache.map.insert(key, report.clone());
    });
    report
}

impl GregoryPatch {
    /// (hits, misses) of the current thread's self-intersection
    /// verdict cache — structured observability for tests asserting
    /// cache behavior (PK-273).
    pub fn self_intersection_cache_stats() -> (u64, u64) {
        CACHE.with(|cache| {
            let cache = cache.borrow();
            (cache.hits, cache.misses)
        })
    }

    /// Clear the current thread's self-intersection verdict cache and
    /// its counters. Test isolation only — production code never
    /// needs invalidation because keys are content-derived.
    pub fn self_intersection_cache_reset() {
        CACHE.with(|cache| {
            *cache.borrow_mut() = CacheState::default();
        });
    }
}
