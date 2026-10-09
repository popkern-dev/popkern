//! Order-preserving parallel map for independent per-item work
//! (the interference check's candidate pairs, PK-749).
//!
//! Native builds split the items into contiguous chunks over scoped
//! threads (one per available core) and concatenate the results in
//! item order, so the outcome is identical to the sequential map.
//! wasm has no threads: there it is the sequential map.

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
        let workers = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
            .min(items.len());
        if workers > 1 {
            let chunk = items.len().div_ceil(workers);
            let f = &f;
            return std::thread::scope(|scope| {
                let handles: Vec<_> = items
                    .chunks(chunk)
                    .map(|part| scope.spawn(move || part.iter().map(f).collect::<Vec<R>>()))
                    .collect();
                handles
                    .into_iter()
                    .flat_map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                    .collect()
            });
        }
    }
    items.iter().map(f).collect()
}
