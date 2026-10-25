//! Order-preserving parallel map over an index range, for independent
//! per-item work whose results are then combined in index order.
//!
//! Natively the range splits into one fixed, contiguous block per
//! worker; each item's result is computed exactly as the sequential
//! map would, and the results come back in index order, so whatever
//! combines them sees the same values in the same order. wasm has no
//! threads: there it is the sequential map.

use std::ops::Range;

use crate::linalg::sparse::available_workers;

/// `range.map(f).collect()`, computed on [`available_workers`] threads
/// when the range has at least `min_items` items; in index order.
pub fn map_range<R, F>(range: Range<usize>, min_items: usize, f: F) -> Vec<R>
where
    R: Send,
    F: Fn(usize) -> R + Sync,
{
    let workers = if range.len() >= min_items.max(2) {
        available_workers()
    } else {
        1
    };
    map_range_on(range, workers, f)
}

/// [`map_range`] on `workers` threads.
pub fn map_range_on<R, F>(range: Range<usize>, workers: usize, f: F) -> Vec<R>
where
    R: Send,
    F: Fn(usize) -> R + Sync,
{
    #[cfg(not(target_arch = "wasm32"))]
    if workers > 1 && range.len() > workers {
        let block = range.len().div_ceil(workers);
        let f = &f;
        return std::thread::scope(|scope| {
            let handles: Vec<_> = (range.start..range.end)
                .step_by(block)
                .map(|start| {
                    let end = (start + block).min(range.end);
                    scope.spawn(move || (start..end).map(f).collect::<Vec<R>>())
                })
                .collect();
            handles
                .into_iter()
                .flat_map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                .collect()
        });
    }
    let _ = workers;
    range.map(f).collect()
}

#[cfg(test)]
#[path = "parallel_tests.rs"]
mod tests;
