//! Thread-local in-place-rewrite switch for the face-family bridge
//! operator: lets the dispatch-level serving-shape gate
//! (PK-205/PK-206) redo one build through the aggregate bridge
//! path without process-global state.

thread_local! {
    static IN_PLACE_REWRITE_DISABLED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// RAII guard: disables the in-place support rewrite on this thread
/// for the guard's lifetime, forcing the aggregate bridge path.
/// Used by the dispatch-level serving-shape gate to redo a build
/// whose in-place handover left glued surfaceless membranes.
pub(crate) struct InPlaceRewriteDisabledGuard;

impl InPlaceRewriteDisabledGuard {
    pub(crate) fn activate() -> Self {
        IN_PLACE_REWRITE_DISABLED.with(|cell| cell.set(true));
        Self
    }
}

impl Drop for InPlaceRewriteDisabledGuard {
    fn drop(&mut self) {
        IN_PLACE_REWRITE_DISABLED.with(|cell| cell.set(false));
    }
}

pub(super) fn in_place_rewrite_disabled() -> bool {
    IN_PLACE_REWRITE_DISABLED.with(|cell| cell.get())
}
