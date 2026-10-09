//! The JS boundary (browser target only). Everything crosses as JSON
//! text or typed arrays; no per-vertex calls (§9). The page's worker
//! module speaks exactly these verbs — `web/live/PROTOCOL.md`
//! documents the message shapes built on them.

use wasm_bindgen::prelude::*;

use popkern_core::time::{monotonic_clock_available, register_monotonic_clock};

use crate::envelope::BodyMesh;
use crate::session::{RunOptions, Session};

pub use crate::PROTOCOL_VERSION;

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

/// A popkern live session as the page holds it.
#[wasm_bindgen]
pub struct LiveSession {
    session: Session,
    bodies: Vec<BodyMesh>,
    /// The zebra analysis meshes of the current model (PK-577), built on
    /// the first `analysis_mesh` request after a run.
    analysis: Vec<BodyMesh>,
    /// The bodies of the last `step_meshes` answer, by row (`None`: the
    /// page already holds that key).
    step: Vec<Option<BodyMesh>>,
}

#[wasm_bindgen]
impl LiveSession {
    /// A session under the budgets in `options_json`
    /// (`{"max_operations": 5000000, "timeout_ms": 10000}`; either
    /// may be null).
    #[wasm_bindgen(constructor)]
    pub fn new(options_json: &str) -> Result<LiveSession, JsValue> {
        let options: RunOptions = if options_json.trim().is_empty() {
            RunOptions::default()
        } else {
            serde_json::from_str(options_json)
                .map_err(|e| JsValue::from_str(&format!("options: {e}")))?
        };
        Ok(LiveSession {
            session: Session::new(options.host()),
            bodies: Vec::new(),
            analysis: Vec::new(),
            step: Vec::new(),
        })
    }

    /// Replace the operation and wall-clock budgets for later runs
    /// (`{"max_operations": …, "timeout_ms": …}`, null = unlimited);
    /// the certification mode and the model are kept. Answers the
    /// budgets now in force as JSON, or the parse error.
    pub fn set_budget(&mut self, options_json: &str) -> Result<String, JsValue> {
        let options: RunOptions = serde_json::from_str(options_json)
            .map_err(|e| JsValue::from_str(&format!("options: {e}")))?;
        let host = options.host();
        self.session.set_options(popkern_rhai::HostOptions {
            max_operations: host.max_operations,
            timeout: host.timeout,
            continuity_certification: self.session.options().continuity_certification,
            track_sources: self.session.options().track_sources,
        });
        Ok(serde_json::json!({
            "max_operations": options.max_operations,
            "timeout_ms": options.timeout_ms,
        })
        .to_string())
    }

    /// Run a script; the envelope as JSON. Body buffers are then
    /// readable by index until the next run or parameter edit.
    pub fn run(&mut self, source: &str, name: &str) -> String {
        let envelope = self.session.run(source, name);
        let json = envelope.to_json().to_string();
        self.bodies = envelope.bodies;
        self.analysis.clear();
        self.step.clear();
        json
    }

    /// Open a STEP file (`Uint8Array`) as a view-only assembly; the
    /// envelope as JSON, buffers readable by index as after `run`.
    pub fn open_step(&mut self, bytes: &[u8], name: &str) -> String {
        let envelope = self.session.open_step(bytes, name);
        let json = envelope.to_json().to_string();
        self.bodies = envelope.bodies;
        self.analysis.clear();
        self.step.clear();
        json
    }

    /// Build a project from `files_json` (`{"files": {path: text},
    /// "base64": {path: data}}`); the envelope as JSON, with a
    /// `project` member holding one row per migration. Checkpoints
    /// persist across calls; buffers readable by index as after `run`.
    pub fn run_project(&mut self, files_json: &str) -> String {
        let envelope = match crate::session_project::decode_project_files(files_json) {
            Ok(files) => self.session.run_project(files),
            Err(detail) => {
                return serde_json::json!({
                    "outcome": "gate_failure", "step": "project", "detail": detail,
                })
                .to_string();
            }
        };
        let json = envelope.to_json().to_string();
        self.bodies = envelope.bodies;
        self.analysis.clear();
        self.step.clear();
        json
    }

    /// Set a declared parameter and replay; the envelope as JSON.
    pub fn set_param(&mut self, name: &str, value: f64) -> String {
        let envelope = self.session.set_param(name, value);
        let json = envelope.to_json().to_string();
        // A joint moves placements only: the meshes (and their analysis
        // and step meshes) stay as they are.
        if !envelope.bodies_unchanged {
            self.bodies = envelope.bodies;
            self.analysis.clear();
            self.step.clear();
        }
        json
    }

    /// Certify deferred blend seams for about `slice_ms` (at least one
    /// seam); the step as JSON: `{certified, remaining, report,
    /// failure}` with `report` set once `remaining` is zero. Body
    /// buffers are untouched — certification never changes geometry.
    pub fn certify(&mut self, slice_ms: f64) -> String {
        let slice = if slice_ms.is_finite() && slice_ms > 0.0 {
            std::time::Duration::from_secs_f64(slice_ms / 1000.0)
        } else {
            std::time::Duration::ZERO
        };
        serde_json::to_value(self.session.certify(slice))
            .unwrap_or(serde_json::Value::Null)
            .to_string()
    }

    /// The document as canonical Rhai, or the outcome JSON as the error.
    pub fn to_rhai(&mut self) -> Result<String, JsValue> {
        self.session
            .to_rhai()
            .map_err(|outcome| outcome_js(&outcome))
    }

    /// The whole model as STEP bytes (`Uint8Array`, every body its own
    /// solid), or the outcome JSON as the error.
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

    /// Show or hide the assembly occurrence at `path`
    /// (`stage/planet_1`); the outcome JSON as the error.
    pub fn set_visible(&mut self, path: &str, visible: bool) -> Result<(), JsValue> {
        self.session
            .set_visible(path, visible)
            .map_err(|outcome| outcome_js(&outcome))
    }

    /// Interference and clearance between assembly occurrences within
    /// `clearance`: the report JSON, or the outcome JSON as the error.
    pub fn interference(&mut self, clearance: f64) -> Result<String, JsValue> {
        self.session
            .interference(clearance)
            .map(|report| {
                serde_json::to_value(report)
                    .unwrap_or(serde_json::Value::Null)
                    .to_string()
            })
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

    /// Re-mesh the live bodies at the zebra analysis bound (6°, PK-577)
    /// and keep the buffers until the next run or parameter edit. The
    /// JSON reply: `{"ok": true, "facet_deg", "triangles", "mesh_ms"}`,
    /// or `{"ok": false, "outcome": ...}`.
    pub fn analysis_mesh(&mut self) -> String {
        let started = performance_now_ms();
        match self.session.analysis_meshes() {
            Ok(bodies) => {
                let triangles: usize = bodies.iter().map(|b| b.indices.len() / 3).sum();
                self.analysis = bodies;
                let elapsed = performance_now_ms() - started;
                let mesh_ms = if elapsed.is_finite() {
                    serde_json::json!(elapsed)
                } else {
                    serde_json::Value::Null
                };
                serde_json::json!({
                    "ok": true,
                    "facet_deg": crate::mesh_buffers::ANALYSIS_MAX_FACET_ANGLE_DEG,
                    "triangles": triangles,
                    "mesh_ms": mesh_ms,
                })
                .to_string()
            }
            Err(outcome) => serde_json::json!({ "ok": false, "outcome": outcome }).to_string(),
        }
    }

    /// Analysis-mesh buffers of body `index`, by `field`: `"positions"`
    /// and `"normals"` (`Float32Array`). Empty before `analysis_mesh`.
    pub fn analysis_f32(&self, index: usize, field: &str) -> Vec<f32> {
        self.analysis
            .get(index)
            .map(|b| match field {
                "positions" => b.positions.clone(),
                "normals" => b.normals.clone(),
                _ => Vec::new(),
            })
            .unwrap_or_default()
    }

    /// Analysis-mesh buffers of body `index`, by `field`: `"indices"`
    /// and `"face_ids"` (`Uint32Array`). Empty before `analysis_mesh`.
    pub fn analysis_u32(&self, index: usize, field: &str) -> Vec<u32> {
        self.analysis
            .get(index)
            .map(|b| match field {
                "indices" => b.indices.clone(),
                "face_ids" => b.face_ids.clone(),
                _ => Vec::new(),
            })
            .unwrap_or_default()
    }

    /// The document's history timeline as JSON: `[{part, steps: [{index,
    /// feature, label, operation, migration, suppressed, failed,
    /// bodies_changed, source}]}]`, or `{"outcome": …}` without a model.
    pub fn timeline(&mut self) -> String {
        match self.session.timeline() {
            Ok(parts) => serde_json::to_string(&parts).unwrap_or_default(),
            Err(outcome) => serde_json::json!({ "outcome": outcome }).to_string(),
        }
    }

    /// The assembly's construction timeline as JSON (`{steps: [{index,
    /// kind, label, instances, migration, source}], recorded}`), `null`
    /// for a part document, or `{"outcome": …}` without a model.
    pub fn assembly_timeline(&self) -> String {
        match self.session.assembly_timeline() {
            Ok(timeline) => serde_json::to_string(&timeline).unwrap_or_default(),
            Err(outcome) => serde_json::json!({ "outcome": outcome }).to_string(),
        }
    }

    /// The occurrences after the first `index` construction steps as
    /// JSON rows (the run's `assembly.instances` shape), or
    /// `{"outcome": …}`.
    pub fn assembly_step(&self, index: usize) -> String {
        match self.session.assembly_step(index) {
            Ok(rows) => serde_json::to_string(&rows).unwrap_or_default(),
            Err(outcome) => serde_json::json!({ "outcome": outcome }).to_string(),
        }
    }

    /// The bodies after the first `index` features of `part` (`""` for
    /// a part document) as JSON rows `{key, body, created_faces,
    /// buffers, mesh}`. `known_json` lists keys the page already holds
    /// (their rows carry no mesh). Buffers of row `i` are then readable
    /// through `step_f32` / `step_u32`.
    pub fn step_meshes(&mut self, part: &str, index: usize, known_json: &str) -> String {
        let known: Vec<String> = serde_json::from_str(known_json).unwrap_or_default();
        let part = (!part.is_empty()).then_some(part);
        match self.session.step_bodies(part, index, &known) {
            Ok(rows) => {
                self.step = rows.iter().map(|r| r.mesh.clone()).collect();
                serde_json::to_string(&rows).unwrap_or_default()
            }
            Err(outcome) => serde_json::json!({ "outcome": outcome }).to_string(),
        }
    }

    /// Step-mesh buffers of row `index`, by `field`: `"positions"`,
    /// `"normals"` (`Float32Array`). Empty for a row without a mesh.
    pub fn step_f32(&self, index: usize, field: &str) -> Vec<f32> {
        self.step
            .get(index)
            .and_then(Option::as_ref)
            .map(|b| match field {
                "positions" => b.positions.clone(),
                "normals" => b.normals.clone(),
                _ => Vec::new(),
            })
            .unwrap_or_default()
    }

    /// Step-mesh buffers of row `index`, by `field`: `"indices"`,
    /// `"face_ids"` (`Uint32Array`). Empty for a row without a mesh.
    pub fn step_u32(&self, index: usize, field: &str) -> Vec<u32> {
        self.step
            .get(index)
            .and_then(Option::as_ref)
            .map(|b| match field {
                "indices" => b.indices.clone(),
                "face_ids" => b.face_ids.clone(),
                _ => Vec::new(),
            })
            .unwrap_or_default()
    }

    /// Analytic per-vertex curvature of body `index`
    /// (`Float32Array`): `kind` is `"gaussian"` or `"mean"`. Empty
    /// for an unknown kind or missing body.
    pub fn body_curvature(&self, index: usize, kind: &str) -> Vec<f32> {
        self.bodies
            .get(index)
            .map(|b| match kind {
                "gaussian" => b.gaussian.clone(),
                "mean" => b.mean.clone(),
                _ => Vec::new(),
            })
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
