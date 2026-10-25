//! The `simulate` and `materials` requests over wasm-bindgen (PK-948):
//! the summary as JSON, then the field meshes' buffers by index, as the
//! body buffers are read after a run.

use wasm_bindgen::prelude::*;

use crate::bindings::LiveSession;

#[wasm_bindgen]
impl LiveSession {
    /// Solve `setup` (analysis setup file text) on the session's part.
    /// The JSON reply: `{"ok": true, "simulation": {status, analysis,
    /// result | error, ...}, "meshes": [{body, field_names}]}`, or
    /// `{"ok": false, "outcome": ...}`. Mesh buffers are then readable
    /// with `simulation_f32` / `simulation_u32` until the next solve.
    pub fn simulate(&mut self, setup: &str) -> String {
        match self.session.simulate(setup) {
            Ok(result) => {
                let meshes: Vec<serde_json::Value> = result
                    .meshes
                    .iter()
                    .map(|m| serde_json::json!({ "body": m.body, "field_names": m.fields.keys().collect::<Vec<_>>() }))
                    .collect();
                self.simulation = result.meshes;
                serde_json::json!({ "ok": true, "simulation": result.summary, "meshes": meshes })
                    .to_string()
            }
            Err(outcome) => serde_json::json!({ "ok": false, "outcome": outcome }).to_string(),
        }
    }

    /// Field-mesh buffer `field` of solid `index` (`Float32Array`):
    /// `"positions"`, `"displacement"`, or `"fields.<name>"`.
    pub fn simulation_f32(&self, index: usize, field: &str) -> Vec<f32> {
        let Some(m) = self.simulation.get(index) else {
            return Vec::new();
        };
        match field {
            "positions" => m.positions.clone(),
            "displacement" => m.displacement.clone(),
            other => other
                .strip_prefix("fields.")
                .and_then(|name| m.fields.get(name))
                .cloned()
                .unwrap_or_default(),
        }
    }

    /// Field-mesh `"indices"` of solid `index` (`Uint32Array`).
    pub fn simulation_u32(&self, index: usize, field: &str) -> Vec<u32> {
        match (self.simulation.get(index), field) {
            (Some(m), "indices") => m.indices.clone(),
            _ => Vec::new(),
        }
    }

    /// The material presets a setup may name, as JSON.
    pub fn materials(&self) -> String {
        crate::session_simulate::materials().to_string()
    }
}
