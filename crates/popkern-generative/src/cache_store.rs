//! Solved designs carried between processes: the cache's entries as
//! JSON (each its key, report, sources and the surface it was faceted
//! from) and back. Importing facets each surface again, which makes the
//! same body (the faceting is deterministic and moves no vertex).

use std::sync::Arc;

use popkern_core::tolerance::Tolerance;
use serde::{Deserialize, Serialize};

use super::{CAPACITY, Entry, Key, StoredSurface, TOUCHED, cache};
use crate::facet::{FacetOptions, facet};
use crate::generate::GenerativeReport;

/// Format of the exported text.
const FORMAT: u32 = 1;

#[derive(Serialize, Deserialize)]
struct StoredEntry {
    /// The key, as two hex digests.
    key: [String; 2],
    report: GenerativeReport,
    sources: Vec<Option<u32>>,
    surface: StoredSurface,
    tolerance: Tolerance,
}

#[derive(Serialize, Deserialize)]
struct Stored {
    format: u32,
    entries: Vec<StoredEntry>,
}

/// Every design the cache holds, as JSON (empty `entries` when none).
pub fn export_cache() -> String {
    export(&|_| true)
}

/// The designs served on this thread since the last call, as JSON
/// (`None` when there were none): what one build used, for a host that
/// hands it to another process (a gallery bake).
pub fn export_touched() -> Option<String> {
    let touched = TOUCHED.with(|t| std::mem::take(&mut *t.borrow_mut()));
    (!touched.is_empty()).then(|| export(&|k| touched.contains(k)))
}

fn export(keep: &dyn Fn(&Key) -> bool) -> String {
    let entries = cache()
        .lock()
        .map(|c| {
            let mut entries: Vec<StoredEntry> = c
                .iter()
                .filter(|(k, _)| keep(k))
                .map(|(k, e)| StoredEntry {
                    key: k.map(|d| format!("{d:016x}")),
                    report: e.report.clone(),
                    sources: e.sources.clone(),
                    surface: e.surface.clone(),
                    tolerance: e.tolerance,
                })
                .collect();
            entries.sort_by(|a, b| a.key.cmp(&b.key));
            entries
        })
        .unwrap_or_default();
    serde_json::to_string(&Stored {
        format: FORMAT,
        entries,
    })
    .unwrap_or_default()
}

/// Adds the designs `text` holds (an [`export_cache`]) to this process's
/// cache. Returns how many were added.
///
/// # Errors
///
/// A message for text of another format, or a surface that does not
/// facet again (never expected: it faceted when exported).
pub fn import_cache(text: &str) -> Result<usize, String> {
    let stored: Stored = serde_json::from_str(text).map_err(|e| e.to_string())?;
    if stored.format != FORMAT {
        return Err(format!(
            "generative cache format {} (expected {FORMAT})",
            stored.format
        ));
    }
    let mut added = Vec::new();
    for e in stored.entries {
        let key: Key = [
            u64::from_str_radix(&e.key[0], 16).map_err(|x| x.to_string())?,
            u64::from_str_radix(&e.key[1], 16).map_err(|x| x.to_string())?,
        ];
        let s = &e.surface;
        let faceted = facet(
            &s.vertices,
            &s.triangles,
            &s.owners,
            &FacetOptions::for_tolerance(&e.tolerance),
            e.tolerance,
        )
        .map_err(|x| x.to_string())?;
        added.push((
            key,
            Arc::new(Entry {
                brep: faceted.brep,
                sources: e.sources,
                report: e.report,
                surface: e.surface,
                tolerance: e.tolerance,
            }),
        ));
    }
    let count = added.len();
    if let Ok(mut c) = cache().lock() {
        if c.len() + count > CAPACITY {
            c.clear();
        }
        c.extend(added);
    }
    Ok(count)
}
