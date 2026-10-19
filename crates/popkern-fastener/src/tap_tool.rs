//! The tapped hole's cutting tool, built once per shape (PK-876).
//!
//! A tapped hole removes the bore, the internal thread's groove and the
//! entry chamfer. Cut one after another, every hole paid the thread
//! boolean (~1.45 s for an M5 against its own fresh bore) and the
//! chamfer cone against the new thread (~0.49 s): about 2 s per hole,
//! 250 s for the ~126 screws of the six-axis arm. The three are one
//! solid in the hole's local frame — their union — and the hole is the
//! body minus that solid placed on the hole axis: one subtraction that
//! meets only the faces at the hole (measured ~55 ms for an M5 in a
//! plate, identical face and edge counts).
//!
//! The tool depends only on its construction inputs, so it is built
//! once per input set and reused, bit for bit, for every hole of that
//! shape: a cache hit returns the tool a miss would have built. Two
//! inputs depend on the body only beyond its material — the chamfer
//! lid height (the body's vertex hull above the entry) and a through
//! hole's reach past the far side. A taller lid or a longer reach
//! removes nothing more (there is no material past the hull), so both
//! are rounded UP to a power-of-two multiple of the pitch, and holes in
//! different bodies share one tool.
//!
//! Any decline (the tool's unions, or the placed subtraction) falls
//! back to the sequential cuts, recorded in [`take_tap_tool_fallbacks`].

use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::FastenerError;
use crate::ThreadSpec;
use crate::engine::{extrude_shape, helical_tool, plan_turns, revolved_tool};
use crate::thread_form::ToolProfile;
use popkern_sketch::ProfileShape;

/// Distinct tools kept before the cache is cleared: an assembly uses a
/// handful of sizes and depths.
const CACHE_CAPACITY: usize = 64;

/// Everything a tapped hole's local-frame tool is built from.
pub(crate) struct TapToolInputs<'a> {
    /// The thread.
    pub spec: &'a ThreadSpec,
    /// The internal thread's ridge (groove cutter) profile.
    pub ridge: &'a ToolProfile,
    /// Bore radius.
    pub r_bore: f64,
    /// Bore headroom above the entry.
    pub up: f64,
    /// Bore depth below the entry.
    pub bore_depth: f64,
    /// Thread start (internal phase) and end, along the axis.
    pub thread_span: (f64, f64),
    /// The entry chamfer's revolved profile, when chamfered.
    pub ring: Option<&'a [(f64, f64)]>,
}

impl TapToolInputs<'_> {
    /// The cache key: every input, f64s by their bits.
    fn key(&self, tolerance: Tolerance) -> String {
        let bits = |x: f64| x.to_bits();
        let ring: Vec<(u64, u64)> = self
            .ring
            .unwrap_or_default()
            .iter()
            .map(|&(r, z)| (bits(r), bits(z)))
            .collect();
        format!(
            "{:?}|{:?}|{}|{}|{}|{}|{}|{ring:?}|{}",
            self.spec,
            self.ridge.points,
            bits(self.r_bore),
            bits(self.up),
            bits(self.bore_depth),
            bits(self.thread_span.0),
            bits(self.thread_span.1),
            bits(tolerance.linear()),
        )
    }
}

/// `x` rounded up to `unit · 2^k` (the smallest such value ≥ `x`, at
/// least `unit`). Used for tool extents beyond the body's material.
pub(crate) fn round_up_pow2(x: f64, unit: f64) -> f64 {
    let mut v = unit;
    while v < x && v.is_finite() {
        v *= 2.0;
    }
    v
}

#[cfg(test)]
thread_local! {
    /// Test switch: cut tapped holes in sequence only (the equivalence
    /// reference for the one-tool cut).
    pub(crate) static SEQUENTIAL_ONLY: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    /// Test probe: tools built on this thread.
    pub(crate) static BUILDS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Whether tapped holes take the one-tool cut first.
pub(crate) fn one_tool_enabled() -> bool {
    #[cfg(test)]
    {
        !SEQUENTIAL_ONLY.with(std::cell::Cell::get)
    }
    #[cfg(not(test))]
    {
        true
    }
}

fn cache() -> &'static Mutex<HashMap<String, Arc<BRep>>> {
    static CACHE: OnceLock<Mutex<HashMap<String, Arc<BRep>>>> = OnceLock::new();
    CACHE.get_or_init(|| Mutex::new(HashMap::new()))
}

/// The tapped hole's tool in its local frame (entry at z = 0, material
/// toward +z), from the cache or built.
pub(crate) fn tap_tool(
    inputs: &TapToolInputs<'_>,
    stage: &str,
    tolerance: Tolerance,
) -> Result<Arc<BRep>, FastenerError> {
    let key = inputs.key(tolerance);
    if let Some(tool) = cache().lock().ok().and_then(|c| c.get(&key).cloned()) {
        return Ok(tool);
    }
    let tool = Arc::new(build(inputs, stage, tolerance)?);
    if let Ok(mut c) = cache().lock() {
        if c.len() >= CACHE_CAPACITY {
            c.clear();
        }
        c.insert(key, Arc::clone(&tool));
    }
    Ok(tool)
}

fn build(
    inputs: &TapToolInputs<'_>,
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    #[cfg(test)]
    BUILDS.with(|b| b.set(b.get() + 1));
    let spec = inputs.spec;
    let mut tool = extrude_shape(
        ProfileShape::Circle {
            radius: inputs.r_bore,
        },
        -inputs.up,
        inputs.up + inputs.bore_depth,
        "tap bore",
        tolerance,
    )?;
    let (z0, turns) = plan_turns(inputs.thread_span.0, inputs.thread_span.1, spec.lead());
    // Every start's groove; a multi-start thread's grooves are disjoint
    // coils, each joined to the bore.
    for start in 0..spec.starts() {
        let groove = helical_tool(spec, inputs.ridge, z0, turns, start, stage, tolerance)?;
        tool = union(&tool, &groove, stage, tolerance)?;
    }
    if let Some(ring) = inputs.ring {
        let ring = revolved_tool(ring, "tap entry chamfer", tolerance)?;
        tool = union(&tool, &ring, stage, tolerance)?;
    }
    Ok(tool)
}

fn union(a: &BRep, b: &BRep, stage: &str, tolerance: Tolerance) -> Result<BRep, FastenerError> {
    popkern_bool::boolean(
        popkern_bool::BoolOperation::Union,
        a.clone(),
        b.clone(),
        tolerance,
        None,
    )
    .map_err(|e| FastenerError::Construction {
        stage: stage.to_string(),
        detail: format!("tap tool union failed: {e}"),
    })
}

/// Why a tapped hole fell back from the one-tool cut to the
/// sequential bore, thread and chamfer cuts.
#[derive(Debug, Clone, PartialEq)]
pub struct TapToolFallback {
    /// The decline's display form.
    pub reason: String,
}

fn fallbacks() -> &'static Mutex<Vec<TapToolFallback>> {
    static LOG: OnceLock<Mutex<Vec<TapToolFallback>>> = OnceLock::new();
    LOG.get_or_init(|| Mutex::new(Vec::new()))
}

/// Record a fallback (and say so on the `fastener` trace).
pub(crate) fn record_fallback(reason: &FastenerError) {
    if popkern_core::trace::trace_enabled("fastener") {
        eprintln!("[fastener] tapped hole: one-tool cut declined ({reason}); cutting in sequence");
    }
    if let Ok(mut log) = fallbacks().lock() {
        log.push(TapToolFallback {
            reason: reason.to_string(),
        });
    }
}

/// Drain the tapped-hole fallbacks recorded so far (process-wide).
pub fn take_tap_tool_fallbacks() -> Vec<TapToolFallback> {
    fallbacks()
        .lock()
        .map(|mut log| std::mem::take(&mut *log))
        .unwrap_or_default()
}

#[cfg(test)]
#[path = "tap_tool_tests.rs"]
mod tests;
