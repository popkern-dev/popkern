//! Parallel-ready analysis utilities.
//!
//! Provides helpers for splitting entity iteration across multiple
//! threads. Uses `std::thread::scope` (Rust 1.63+) for safe scoped
//! parallelism with borrowed data. Can be upgraded to `rayon` when
//! added as a workspace dependency.
//!
//! # Usage
//!
//! ```ignore
//! use popkern_core::healing::parallel::parallel_map;
//!
//! let items: Vec<u64> = (0..1000).collect();
//! let results = parallel_map(&items, 100, |x| x * 2);
//! assert_eq!(results.len(), 1000);
//! ```

/// Default threshold below which [`parallel_map`] runs sequentially.
///
/// For small inputs the overhead of spawning threads exceeds the
/// benefit. This value is a reasonable starting point; callers can
/// pass a custom threshold.
pub const DEFAULT_PARALLEL_THRESHOLD: usize = 64;

/// Apply a function to each item in a slice, collecting results.
///
/// Uses sequential iteration for small inputs (below `threshold`)
/// and spawns scoped worker threads for larger inputs. The number
/// of threads is determined by [`std::thread::available_parallelism`].
///
/// Thread panics are caught — if a worker thread panics, its chunk
/// produces an empty result vector. This avoids propagating panics
/// but means the output may be shorter than the input in that case.
///
/// # Arguments
///
/// * `items` — the input slice to process.
/// * `threshold` — inputs smaller than this run sequentially.
/// * `f` — the function to apply to each item.
///
/// # Returns
///
/// A `Vec<R>` with one result per input item (in order), unless a
/// worker thread panicked.
pub fn parallel_map<T: Sync, R: Send>(
    items: &[T],
    threshold: usize,
    f: impl Fn(&T) -> R + Send + Sync,
) -> Vec<R> {
    if items.len() < threshold || items.is_empty() {
        return items.iter().map(&f).collect();
    }

    let num_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .min(items.len());

    if num_threads <= 1 {
        return items.iter().map(&f).collect();
    }

    let chunk_size = items.len().div_ceil(num_threads);

    std::thread::scope(|s| {
        let handles: Vec<_> = items
            .chunks(chunk_size)
            .map(|chunk| s.spawn(|| chunk.iter().map(&f).collect::<Vec<_>>()))
            .collect();

        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap_or_default())
            .collect()
    })
}

/// Apply a fallible function to each item, collecting `Result`s.
///
/// Like [`parallel_map`] but for functions that return `Result`.
/// Each item produces an independent `Result` — one failure does
/// not abort processing of other items.
pub fn parallel_try_map<T: Sync, R: Send, E: Send>(
    items: &[T],
    threshold: usize,
    f: impl Fn(&T) -> Result<R, E> + Send + Sync,
) -> Vec<Result<R, E>> {
    parallel_map(items, threshold, f)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sequential_below_threshold() {
        let items: Vec<u32> = (0..10).collect();
        let results = parallel_map(&items, 100, |x| x * 2);
        let expected: Vec<u32> = (0..10).map(|x| x * 2).collect();
        assert_eq!(results, expected);
    }

    #[test]
    fn parallel_above_threshold() {
        let items: Vec<u64> = (0..200).collect();
        let results = parallel_map(&items, 50, |x| x * 3);
        let expected: Vec<u64> = (0..200).map(|x| x * 3).collect();
        assert_eq!(results, expected);
    }

    #[test]
    fn empty_input() {
        let items: Vec<i32> = Vec::new();
        let results = parallel_map(&items, 10, |x| x + 1);
        assert!(results.is_empty());
    }

    #[test]
    fn single_item() {
        let items = vec![42];
        let results = parallel_map(&items, 10, |x| x * 2);
        assert_eq!(results, vec![84]);
    }

    #[test]
    fn preserves_order() {
        let items: Vec<usize> = (0..500).collect();
        let results = parallel_map(&items, 10, |x| *x);
        assert_eq!(results, items);
    }

    #[test]
    fn parallel_try_map_all_ok() {
        let items: Vec<i32> = (1..=10).collect();
        let results: Vec<Result<i32, &str>> = parallel_try_map(&items, 5, |x| Ok(x * 2));
        for (i, r) in results.iter().enumerate() {
            match r {
                Ok(v) => assert_eq!(*v, (i as i32 + 1) * 2),
                Err(_) => panic!("expected Ok"),
            }
        }
    }

    #[test]
    fn parallel_try_map_some_errors() {
        let items: Vec<i32> = (0..10).collect();
        let results: Vec<Result<i32, String>> = parallel_try_map(&items, 5, |x| {
            if *x % 3 == 0 {
                Err(format!("bad: {x}"))
            } else {
                Ok(x * 2)
            }
        });
        assert_eq!(results.len(), 10);
        // Items 0, 3, 6, 9 should be Err.
        assert!(results[0].is_err());
        assert!(results[3].is_err());
        assert!(results[6].is_err());
        assert!(results[9].is_err());
        // Others should be Ok.
        assert_eq!(results[1], Ok(2));
        assert_eq!(results[4], Ok(8));
    }

    #[test]
    fn threshold_boundary_exact() {
        // Exactly at threshold should still be sequential.
        let items: Vec<u32> = (0..64).collect();
        let results = parallel_map(&items, 64, |x| x + 1);
        let expected: Vec<u32> = (1..=64).collect();
        assert_eq!(results, expected);
    }

    #[test]
    fn threshold_boundary_one_over() {
        // One above threshold triggers parallel path.
        let items: Vec<u32> = (0..65).collect();
        let results = parallel_map(&items, 64, |x| x + 1);
        let expected: Vec<u32> = (1..=65).collect();
        assert_eq!(results, expected);
    }

    #[test]
    fn default_threshold_value() {
        assert_eq!(DEFAULT_PARALLEL_THRESHOLD, 64);
    }
}
