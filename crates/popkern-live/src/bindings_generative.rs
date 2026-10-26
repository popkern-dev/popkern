//! Generative designs in the browser (PK-956): a design takes minutes
//! to solve in wasm (a crisp design runs 250+ updates: 21 s for 15 at 2 mm
//! cells on a bracket, against 11 s native),
//! so builds here never solve one on their own. A design not solved yet
//! declines (`GenerativeDeclined`, kind `NeedsSolve`) and the page keeps
//! the model it shows; the page's re-solve runs one build with solving
//! on, its updates reported to the page's progress callback.

use std::ops::ControlFlow;

use wasm_bindgen::prelude::*;

use crate::bindings::LiveSession;
use crate::session::Session;

#[wasm_bindgen]
impl LiveSession {
    /// Whether the next builds may solve generative designs (the page's
    /// re-solve turns it on for one build).
    pub fn set_generative_solve(&mut self, solve: bool) {
        self.generative_solve = solve;
    }

    /// Adds a bake's solved generative designs (its `generative`, as
    /// JSON text) so builds serve them without solving. Returns how many
    /// were added, or throws the reason.
    pub fn seed_generative(&mut self, designs: &str) -> Result<usize, JsValue> {
        self.session
            .seed_generative(designs)
            .map_err(|e| JsValue::from_str(&e))
    }

    /// The function called with each generative update during a build,
    /// as JSON `{"count", "max", "compliance_nmm", "volume_fraction",
    /// "change"}` (`null` clears it).
    pub fn set_progress(&mut self, callback: JsValue) {
        self.progress = callback.dyn_into::<js_sys::Function>().ok();
    }
}

impl LiveSession {
    /// Runs one build on the session with this host's generative policy
    /// and progress callback in force.
    pub(crate) fn building<R>(&mut self, build: impl FnOnce(&mut Session) -> R) -> R {
        let progress = self.progress.clone();
        let observer = move |u: &popkern_generative::Update| {
            if let Some(f) = &progress {
                let json = serde_json::json!({
                    "count": u.count,
                    "max": u.max,
                    "compliance_nmm": u.compliance,
                    "volume_fraction": u.volume_fraction,
                    "change": u.change,
                });
                let _ = f.call1(&JsValue::NULL, &JsValue::from_str(&json.to_string()));
            }
            ControlFlow::Continue(())
        };
        // An explicit re-solve runs as long as the design takes: the
        // page's budgets would stop it at the first script step after.
        let solve = self.generative_solve;
        let budgets = self.session.options();
        if solve {
            self.session.set_options(popkern_rhai::HostOptions {
                max_operations: None,
                timeout: None,
                ..budgets
            });
        }
        let session = &mut self.session;
        let out = popkern_generative::observe_with(observer, || {
            popkern_generative::defer_with(!solve, || build(session))
        });
        if solve {
            self.session.set_options(budgets);
        }
        out
    }
}
