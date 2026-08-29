//! Topic-based diagnostic tracing.
//!
//! One env var replaces the historical per-subsystem trace flags:
//! `POPKERN_TRACE=topic1,topic2` enables the named topics,
//! `POPKERN_TRACE=all` enables everything. Topic names are the
//! semantic subsystem names (`boolean`, `sphere-fit`, `fragment-walk`,
//! …) — see `docs/reference/debug-flags.md` for the registry.
//!
//! The topic set is parsed once per process (`OnceLock`): reads are a
//! lock-free set lookup, so hot loops can gate on `trace_enabled`
//! without taking the process env lock per item (which measurably
//! serialized parallel pipeline stages when traces were raw env
//! reads). Tests that need a topic must set the env var before the
//! process's first `trace_enabled` call.

use std::collections::BTreeSet;
use std::sync::OnceLock;

fn topics() -> &'static Option<BTreeSet<String>> {
    static TOPICS: OnceLock<Option<BTreeSet<String>>> = OnceLock::new();
    TOPICS.get_or_init(|| {
        let raw = std::env::var("POPKERN_TRACE").ok()?;
        Some(
            raw.split(',')
                .map(|t| t.trim().to_ascii_lowercase())
                .filter(|t| !t.is_empty())
                .collect(),
        )
    })
}

/// Whether diagnostic tracing is enabled for `topic`.
///
/// `POPKERN_TRACE=all` enables every topic. Topic matching is exact
/// (after ASCII-lowercasing the env list); topics are kebab-case
/// subsystem names.
pub fn trace_enabled(topic: &str) -> bool {
    match topics() {
        Some(set) => set.contains("all") || set.contains(topic),
        None => false,
    }
}
