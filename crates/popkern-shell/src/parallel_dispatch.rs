//! Ordered parallel dispatch for pure per-item pipeline stages
//! (PK-329, PK-277 discipline).
//!
//! Every parallel site in this crate follows the same contract: the
//! task list is deterministic, each task reads shared state immutably
//! (plain `SlotMap` storage, no interior mutability), and results are
//! collected positionally — so output is bit-identical to the
//! sequential loop at any thread count, and opting out of `parallel`
//! is a resource decision, not a correctness one. Stages that honor an
//! env-gated trace pass `serial: true` while tracing so interleaved
//! thread output cannot garble the dump.

/// Map `f` over `items`, collecting results in item order. Runs on the
/// rayon pool unless the `parallel` feature is off or `serial` is set.
pub(crate) fn dispatch_ordered<I: Sync, T: Send>(
    items: &[I],
    serial: bool,
    f: impl Fn(&I) -> T + Send + Sync,
) -> Vec<T> {
    #[cfg(feature = "parallel")]
    {
        if !serial {
            use rayon::prelude::*;
            return items.par_iter().map(f).collect();
        }
    }
    let _ = serial;
    items.iter().map(f).collect()
}
