//! Designs already solved, and the host's view of one being solved.
//!
//! A script re-runs from the top on every edit, so a design whose
//! inputs did not change is served from here: the key is a 128-bit
//! digest of everything the solve reads (the design space's analysis
//! mesh to the bit, the setup, the units and tolerance), so a hit is the
//! same design the solve would make (it is deterministic, GD6, and the
//! same in wasm as natively, DT1), never a stale one. Sources are kept
//! as analysis face ids and mapped back to the design space's handles on
//! each hit. A few designs are kept; past that the cache empties (a bake
//! or a session edits one design at a time). Each entry keeps the
//! surface it was faceted from, so a host can carry solved designs to
//! another process ([`export_cache`], [`import_cache`]: a gallery bake
//! hands the page the designs it solved natively).
//!
//! Hosts that show progress install an observer for the duration of a
//! build ([`observe_with`]); it sees every update of every solve on its
//! thread and may stop the solve (`Cancelled`). A host where a solve is
//! too slow to run on every edit (popkern live in the browser: a crisp design
//! runs 250+ updates, minutes in wasm, PK-956) builds with solving deferred ([`defer_with`]):
//! a design not in the cache then declines as `NeedsSolve`, and the
//! host asks for an explicit solve.

use std::cell::RefCell;
use std::collections::HashMap;
use std::ops::ControlFlow;
use std::sync::{Arc, Mutex, OnceLock};

use popkern_analysis::structural::generative::Progress;
use popkern_core::topology::BRep;

use crate::generate::{
    GenerativeError, GenerativeReport, GenerativeRequest, GenerativeResult, Prepared, prepare,
    solve,
};

/// Designs kept.
const CAPACITY: usize = 4;

/// The surface a design body is faceted from, in the body's units:
/// triangles of the kept pieces, each with its owner (an analysis face
/// id, or `faces` and above for the design's own surface).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct StoredSurface {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangles.
    pub triangles: Vec<[u32; 3]>,
    /// Each triangle's owner.
    pub owners: Vec<u32>,
    /// The design space's analysis faces (owners at or above are the
    /// design's own surface).
    pub faces: u32,
}

/// A solved design, sources as analysis face ids.
struct Entry {
    brep: BRep,
    sources: Vec<Option<u32>>,
    report: GenerativeReport,
    surface: StoredSurface,
    tolerance: popkern_core::tolerance::Tolerance,
}

/// A digest of everything a solve reads.
type Key = [u64; 2];

type Cache = Mutex<HashMap<Key, Arc<Entry>>>;

fn cache() -> &'static Cache {
    static CACHE: OnceLock<Cache> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(HashMap::new()))
}

/// One optimisation update as a host sees it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Update {
    /// Updates so far.
    pub count: usize,
    /// The most the solve will run.
    pub max: usize,
    /// Compliance before the update, N·mm.
    pub compliance: f64,
    /// Material volume over the design space's.
    pub volume_fraction: f64,
    /// The largest density change of the update (the solve stops when
    /// it falls under 0.01).
    pub change: f64,
}

type Observer = Box<dyn FnMut(&Update) -> ControlFlow<()>>;

thread_local! {
    static OBSERVER: RefCell<Option<Observer>> = const { RefCell::new(None) };
    /// Designs served on this thread since the host last asked
    /// ([`store::export_touched`]).
    static TOUCHED: RefCell<Vec<Key>> = const { RefCell::new(Vec::new()) };
    static DEFERRED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Runs `f` with solving deferred (`defer`) or allowed on this thread:
/// deferred, a design not in the cache declines as
/// [`GenerativeError::NeedsSolve`] instead of solving.
pub fn defer_with<R>(defer: bool, f: impl FnOnce() -> R) -> R {
    let previous = DEFERRED.with(|d| d.replace(defer));
    let out = f();
    DEFERRED.with(|d| d.set(previous));
    out
}

/// Runs `f` with `observer` seeing every update of every design solved
/// on this thread meanwhile (it may stop the solve with `Break`).
pub fn observe_with<R>(
    observer: impl FnMut(&Update) -> ControlFlow<()> + 'static,
    f: impl FnOnce() -> R,
) -> R {
    let previous = OBSERVER.with(|o| o.borrow_mut().replace(Box::new(observer)));
    let out = f();
    OBSERVER.with(|o| *o.borrow_mut() = previous);
    out
}

fn notify(update: &Update) -> ControlFlow<()> {
    OBSERVER.with(|o| match o.borrow_mut().as_mut() {
        Some(observer) => observer(update),
        None => ControlFlow::Continue(()),
    })
}

/// Two independent FNV-1a 64 digests of `bytes` (a 128-bit key: a
/// collision between the handful of designs a cache holds is beyond
/// any practical chance).
fn digest(bytes: &[u8]) -> Key {
    let mut a: u64 = 0xcbf2_9ce4_8422_2325;
    let mut b: u64 = 0x6c62_272e_07bb_0142;
    for &byte in bytes {
        a = (a ^ u64::from(byte)).wrapping_mul(0x0100_0000_01b3);
        b = (b ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3 ^ 0x9e37);
    }
    [a, b]
}

/// Everything a solve reads, digested.
fn key(prepared: &Prepared, request: &GenerativeRequest, space: &BRep) -> Key {
    let mut k = Vec::new();
    for solid in &prepared.named.model.solids {
        for v in &solid.vertices {
            for x in v {
                k.extend_from_slice(&x.to_bits().to_le_bytes());
            }
        }
        for t in &solid.triangles {
            for i in t {
                k.extend_from_slice(&i.to_le_bytes());
            }
        }
        for f in &solid.triangle_face {
            k.extend_from_slice(&f.0.to_le_bytes());
        }
    }
    // Debug renders every f64 to its shortest round-trip form.
    k.extend_from_slice(
        format!(
            "{:?}|{:?}|{:?}",
            prepared.setup, request.mm_per_unit, space.tolerance
        )
        .as_bytes(),
    );
    digest(&k)
}

/// [`crate::generate::generate`] through the cache: an unchanged
/// design is served without solving; a solved one is kept. Each update
/// of a solve goes to the thread's observer ([`observe_with`]).
///
/// # Errors
///
/// As [`crate::generate::generate`].
pub fn generate_cached(
    space: &BRep,
    request: &GenerativeRequest,
) -> Result<(GenerativeResult, bool), GenerativeError> {
    let prepared = prepare(space, request)?;
    let k = key(&prepared, request, space);
    let hit = cache().lock().ok().and_then(|c| c.get(&k).cloned());
    let entry = match hit {
        Some(entry) => (entry, true),
        None if DEFERRED.with(std::cell::Cell::get) => {
            return Err(GenerativeError::NeedsSolve);
        }
        None => {
            let max = request.max_iterations;
            let mut observer = |p: Progress<'_>| {
                notify(&Update {
                    count: p.count,
                    max,
                    compliance: p.iteration.compliance,
                    volume_fraction: p.iteration.volume_fraction,
                    change: p.iteration.change,
                })
            };
            let (solved, _) = solve(
                &prepared,
                request.mm_per_unit,
                space.tolerance,
                &mut observer,
            )?;
            let entry = Arc::new(Entry {
                brep: solved.brep,
                sources: solved.sources,
                report: solved.report,
                surface: solved.surface,
                tolerance: space.tolerance,
            });
            if let Ok(mut c) = cache().lock() {
                if c.len() >= CAPACITY {
                    c.clear();
                }
                c.insert(k, Arc::clone(&entry));
            }
            (entry, false)
        }
    };
    let (entry, hit) = entry;
    TOUCHED.with(|t| {
        let mut t = t.borrow_mut();
        if !t.contains(&k) {
            t.push(k);
        }
    });
    Ok((
        GenerativeResult {
            brep: entry.brep.clone(),
            sources: prepared.handles(&entry.sources),
            report: entry.report.clone(),
        },
        hit,
    ))
}

#[path = "cache_store.rs"]
mod store;
pub use store::{export_cache, export_touched, import_cache};
