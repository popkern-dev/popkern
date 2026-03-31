//! Per-step measurement cache for the replay engine.
//!
//! Caches measurement results during a single replay step. Invalidated
//! when the BRep changes (after any non-Reference feature executes).

use std::collections::HashMap;

use crate::error::ReplayError;
use crate::measure_query::{MeasureContext, MeasureQuery};

/// A cache for measurement results during replay.
///
/// Results are keyed by a string hash of the query. The cache is invalidated
/// by calling [`invalidate`](Self::invalidate) after any BRep-modifying feature.
pub(crate) struct MeasureCache {
    entries: HashMap<String, f64>,
    generation: u64,
}

impl MeasureCache {
    /// Create a new empty cache.
    pub fn new() -> Self {
        Self {
            entries: HashMap::new(),
            generation: 0,
        }
    }

    /// Invalidate all cached entries (call after BRep changes).
    pub fn invalidate(&mut self) {
        self.entries.clear();
        self.generation += 1;
    }

    /// Look up or compute a measurement result.
    ///
    /// If the query is cached, returns the cached value. Otherwise calls
    /// `resolve_and_evaluate`, caches the result, and returns it.
    pub fn get_or_evaluate(
        &mut self,
        query: &MeasureQuery,
        ctx: &MeasureContext<'_>,
    ) -> Result<f64, ReplayError> {
        let key = format!("{query:?}");
        if let Some(&value) = self.entries.get(&key) {
            return Ok(value);
        }
        let value = crate::measure_resolve::resolve_and_evaluate(query, ctx)?;
        self.entries.insert(key, value);
        Ok(value)
    }
}
