//! Process-local dispatch overrides for tests and diagnostics.
//!
//! These replace the load-bearing env flags the test suite used to
//! set (`POPKERN_SPHERE_FIT_BOUNDARY`, `POPKERN_KP13_BI5`, …): a
//! typed tri-state per dispatch decision, set through an RAII guard
//! that restores the previous state on drop. Unlike env vars, the
//! overrides are data-race-free under parallel in-process tests and
//! cannot leak into child processes.
//!
//! Production code never sets these; every slot's `None` state maps
//! to the production default at its single read site.

use std::sync::atomic::{AtomicU8, Ordering};

const UNSET: u8 = 0;
const ON: u8 = 1;
const OFF: u8 = 2;

/// One tri-state dispatch override slot (`None` = production default).
pub struct DispatchOverride {
    state: AtomicU8,
}

impl DispatchOverride {
    const fn new() -> Self {
        Self {
            state: AtomicU8::new(UNSET),
        }
    }

    /// Current override: `None` when unset (production default).
    pub fn get(&self) -> Option<bool> {
        match self.state.load(Ordering::Relaxed) {
            ON => Some(true),
            OFF => Some(false),
            _ => None,
        }
    }

    fn raw_set(&self, value: Option<bool>) -> u8 {
        let next = match value {
            None => UNSET,
            Some(true) => ON,
            Some(false) => OFF,
        };
        self.state.swap(next, Ordering::Relaxed)
    }

    /// Set the override with no restore — for test binaries whose
    /// every test wants the same value (a per-test RAII restore would
    /// unset it under a still-running parallel sibling).
    pub fn set(&self, value: Option<bool>) {
        self.raw_set(value);
    }

    /// Set the override for the guard's lifetime; the previous state
    /// is restored on drop.
    pub fn set_guard(&'static self, value: bool) -> DispatchOverrideGuard {
        let previous = self.raw_set(Some(value));
        DispatchOverrideGuard {
            slot: self,
            previous,
        }
    }
}

/// RAII restore for [`DispatchOverride::set_guard`].
pub struct DispatchOverrideGuard {
    slot: &'static DispatchOverride,
    previous: u8,
}

impl Drop for DispatchOverrideGuard {
    fn drop(&mut self) {
        self.slot.state.store(self.previous, Ordering::Relaxed);
    }
}

/// V-free repeated-apex contract: `Some(false)` disables the V-free
/// requirement (mutating-builder acceptance path), `Some(true)`
/// forces it at every valence, `None` keeps the DefaultOn valence
/// scoping.
pub static V_FREE_REPEATED_APEX: DispatchOverride = DispatchOverride::new();

/// Sphere-fit boundary substitution: `Some(true)` arms the sphere-fit
/// termination substitution globally (production default: off unless
/// the per-call `BlendBoundaryOptions` asks for it).
pub static SPHERE_FIT_BOUNDARY: DispatchOverride = DispatchOverride::new();

/// KP13 bi-5 macro-cap dispatch: `Some(true)` routes eligible caps
/// through the bi-5 promotion (production default: off — the bi-6
/// polar path serves N ≥ 5).
pub static KP13_BI5_DISPATCH: DispatchOverride = DispatchOverride::new();

/// Outer-fillet G1 rung: `Some(true)` extends the binonic optimiser's
/// L-M residual with cap-vs-fillet G1 samples along the outer
/// boundary (production default: off, PK-275 armed-rung contract).
pub static OUTER_FILLET_G1: DispatchOverride = DispatchOverride::new();

/// Setback surgery force-off: `Some(true)` disables the setback
/// vertex blend so tests can pin the fallback dispatch
/// (production default: setback serves eligible apexes).
pub static SETBACK_DISABLE: DispatchOverride = DispatchOverride::new();

/// Cylinder off-surface-vertex topology repair: `Some(true)` runs the
/// repair everywhere, `Some(false)` disables it even for apex
/// boundary contexts (production default: apex-boundary contexts
/// only).
pub static TOPO_REPAIR: DispatchOverride = DispatchOverride::new();

/// Different-planes arm of the topology repair: `Some(true)` enables
/// the cross-plane chord-pair repair (production default: off — it
/// once regressed N ≥ 8 pyramids).
pub static TOPO_REPAIR_DIFF: DispatchOverride = DispatchOverride::new();
