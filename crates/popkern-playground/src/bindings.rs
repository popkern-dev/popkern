//! The JS boundary (browser target only). Everything crosses as JSON
//! text or typed arrays; no per-vertex calls (§9). The page's worker
//! module speaks exactly these verbs — `web/playground/PROTOCOL.md`
//! documents the message shapes built on them.

use wasm_bindgen::prelude::*;

use popkern_core::time::{monotonic_clock_available, register_monotonic_clock};

use crate::envelope::BodyMesh;
use crate::session::{RunOptions, Session};

/// The boundary's protocol version. Bumped when a verb's shape
/// changes; the page checks it before speaking.
pub const PROTOCOL_VERSION: &str = "1";

/// `performance.now()` through the global object (window, worker or
/// node all expose it); falls back to `Date.now()` — still monotonic
/// enough for a budget — when it is missing.
fn performance_now_ms() -> f64 {
    let global = js_sys::global();
    let performance = js_sys::Reflect::get(&global, &JsValue::from_str("performance"))
        .ok()
        .filter(|v| !v.is_undefined() && !v.is_null());
    let now_ms = performance.and_then(|performance| {
        let now = js_sys::Reflect::get(&performance, &JsValue::from_str("now")).ok()?;
        let now = now.dyn_ref::<js_sys::Function>()?;
        now.call0(&performance).ok()?.as_f64()
    });
    now_ms.unwrap_or_else(js_sys::Date::now)
}

/// Register the host clock so stage telemetry and the wall-clock
/// budget are live. Returns whether a monotonic clock is available
/// afterwards. Call once, before the first run.
#[wasm_bindgen]
pub fn init() -> bool {
    let _ = register_monotonic_clock(performance_now_ms);
    monotonic_clock_available()
}

/// The protocol version of this module.
#[wasm_bindgen]
pub fn protocol_version() -> String {
    PROTOCOL_VERSION.to_string()
}

/// A playground session as the page holds it.
#[wasm_bindgen]
pub struct Playground {
    session: Session,
    bodies: Vec<BodyMesh>,
}

#[wasm_bindgen]
impl Playground {
    /// A session under the budgets in `options_json`
    /// (`{"max_operations": 5000000, "timeout_ms": 10000}`; either
    /// may be null).
    #[wasm_bindgen(constructor)]
    pub fn new(options_json: &str) -> Result<Playground, JsValue> {
        let options: RunOptions = if options_json.trim().is_empty() {
            RunOptions::default()
        } else {
            serde_json::from_str(options_json)
                .map_err(|e| JsValue::from_str(&format!("options: {e}")))?
        };
        Ok(Playground {
            session: Session::new(options.host()),
            bodies: Vec::new(),
        })
    }

    /// Run a script; the envelope as JSON. Body buffers are then
    /// readable by index until the next run or parameter edit.
    pub fn run(&mut self, source: &str, name: &str) -> String {
        let envelope = self.session.run(source, name);
        let json = envelope.to_json().to_string();
        self.bodies = envelope.bodies;
        json
    }

    /// Set a declared parameter and replay; the envelope as JSON.
    pub fn set_param(&mut self, name: &str, value: f64) -> String {
        let envelope = self.session.set_param(name, value);
        let json = envelope.to_json().to_string();
        self.bodies = envelope.bodies;
        json
    }

    /// The document as canonical Rhai, or the outcome JSON as the error.
    pub fn to_rhai(&mut self) -> Result<String, JsValue> {
        self.session
            .to_rhai()
            .map_err(|outcome| outcome_js(&outcome))
    }

    /// The active body as STEP bytes (`Uint8Array`), or the outcome
    /// JSON as the error.
    pub fn export_step(&mut self, name: &str) -> Result<Vec<u8>, JsValue> {
        self.session
            .export_step(name)
            .map_err(|outcome| outcome_js(&outcome))
    }

    /// The document as `.pkn` bytes (`Uint8Array`), or the outcome
    /// JSON as the error.
    pub fn export_pkn(&self) -> Result<Vec<u8>, JsValue> {
        self.session
            .export_pkn()
            .map_err(|outcome| outcome_js(&outcome))
    }

    /// Bodies the last envelope described.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// `x y z` per vertex of body `index` (`Float32Array`).
    pub fn body_positions(&self, index: usize) -> Vec<f32> {
        self.bodies
            .get(index)
            .map(|b| b.positions.clone())
            .unwrap_or_default()
    }

    /// `nx ny nz` per vertex of body `index` (`Float32Array`).
    pub fn body_normals(&self, index: usize) -> Vec<f32> {
        self.bodies
            .get(index)
            .map(|b| b.normals.clone())
            .unwrap_or_default()
    }

    /// Triangle vertex indices of body `index` (`Uint32Array`).
    pub fn body_indices(&self, index: usize) -> Vec<u32> {
        self.bodies
            .get(index)
            .map(|b| b.indices.clone())
            .unwrap_or_default()
    }

    /// Face table index per triangle of body `index` (`Uint32Array`).
    pub fn body_face_ids(&self, index: usize) -> Vec<u32> {
        self.bodies
            .get(index)
            .map(|b| b.face_ids.clone())
            .unwrap_or_default()
    }
}

fn outcome_js(outcome: &popkern_rhai::Outcome) -> JsValue {
    JsValue::from_str(
        &serde_json::to_value(outcome)
            .unwrap_or(serde_json::Value::Null)
            .to_string(),
    )
}
