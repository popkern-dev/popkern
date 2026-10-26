//! The histories of solved designs (PK-959), kept beside the cache by
//! design id: recorded by every solve on this process, or imported from
//! another's (a gallery bake's), and read by a host showing a design
//! converging.

use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

use super::{CAPACITY, Key, TOUCHED, design_id, parse_design_id};
use crate::history::DesignHistory;

type Histories = Mutex<HashMap<Key, Arc<DesignHistory>>>;

fn histories() -> &'static Histories {
    static HISTORIES: OnceLock<Histories> = OnceLock::new();
    HISTORIES.get_or_init(|| Mutex::new(HashMap::new()))
}

/// Keeps `history` for the design `key` (a few are kept, as designs
/// are).
pub(super) fn keep(key: Key, history: DesignHistory) {
    if let Ok(mut h) = histories().lock() {
        if h.len() >= CAPACITY && !h.contains_key(&key) {
            h.clear();
        }
        h.insert(key, Arc::new(history));
    }
}

/// The history of the design `design` (its id), when this process
/// solved it or imported its history.
pub fn design_history(design: &str) -> Option<Arc<DesignHistory>> {
    let key = parse_design_id(design)?;
    histories().lock().ok()?.get(&key).cloned()
}

/// The ids of the designs served on this thread since the host last
/// took them ([`crate::export_touched`]), without taking them.
pub fn touched_designs() -> Vec<String> {
    TOUCHED.with(|t| t.borrow().iter().map(|&k| design_id(k)).collect())
}

/// The histories of `designs` this process holds, as bytes
/// (`history_store` layout); `None` when it holds none of them.
pub fn export_histories(designs: &[String]) -> Option<Vec<u8>> {
    let held: Vec<(String, Arc<DesignHistory>)> = designs
        .iter()
        .filter_map(|d| Some((d.clone(), design_history(d)?)))
        .collect();
    if held.is_empty() {
        return None;
    }
    let refs: Vec<(String, &DesignHistory)> =
        held.iter().map(|(d, h)| (d.clone(), h.as_ref())).collect();
    Some(crate::history_store::encode(&refs))
}

/// Keeps the histories `bytes` hold (an [`export_histories`]). Returns
/// how many.
///
/// # Errors
///
/// A message for bytes that are not a design history, or a design id
/// that is not one.
pub fn import_histories(bytes: &[u8]) -> Result<usize, String> {
    let decoded = crate::history_store::decode(bytes)?;
    let count = decoded.len();
    for (design, history) in decoded {
        let key =
            parse_design_id(&design).ok_or_else(|| format!("'{design}' is not a design id"))?;
        keep(key, history);
    }
    Ok(count)
}
