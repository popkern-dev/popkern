//! Order-preserving parallel map for independent per-item work (the
//! interference check's candidate pairs, PK-749; the watertightness
//! gate's definitions, PK-725).
//!
//! Native builds run one scoped thread per available core; each pulls
//! the next unclaimed item (dynamic scheduling: a few expensive items —
//! the arm's threaded screws — must not land in one contiguous chunk)
//! and the results are put back in item order, so the outcome is
//! identical to the sequential map. wasm has no threads: there it is
//! the sequential map.

/// `items.iter().map(f).collect()`, computed in parallel where threads
/// exist; the result is in item order.
pub(crate) fn ordered<T, R, F>(items: &[T], f: F) -> Vec<R>
where
    T: Sync,
    R: Send,
    F: Fn(&T) -> R + Sync,
{
    #[cfg(not(target_arch = "wasm32"))]
    {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let workers = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
            .min(items.len());
        if workers > 1 {
            let next = AtomicUsize::new(0);
            let (f, next) = (&f, &next);
            let mut done: Vec<(usize, R)> = std::thread::scope(|scope| {
                let handles: Vec<_> = (0..workers)
                    .map(|_| {
                        scope.spawn(move || {
                            let mut out = Vec::new();
                            loop {
                                let i = next.fetch_add(1, Ordering::Relaxed);
                                let Some(item) = items.get(i) else {
                                    break;
                                };
                                out.push((i, f(item)));
                            }
                            out
                        })
                    })
                    .collect();
                handles
                    .into_iter()
                    .flat_map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                    .collect()
            });
            done.sort_unstable_by_key(|(i, _)| *i);
            return done.into_iter().map(|(_, r)| r).collect();
        }
    }
    items.iter().map(f).collect()
}
