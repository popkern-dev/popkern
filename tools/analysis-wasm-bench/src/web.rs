//! The JavaScript entry point (wasm32 only).

use wasm_bindgen::prelude::wasm_bindgen;

/// Solves the case in `bytes` under `max_active_cells` and returns
/// `hash cells dof iterations` (hash in hex), or `error: …`.
#[wasm_bindgen]
pub fn run_case(bytes: &[u8], max_active_cells: usize) -> String {
    match super::run(bytes, max_active_cells) {
        Ok(o) => format!("{:016x} {} {} {}", o.hash, o.cells, o.dof, o.iterations),
        Err(e) => format!("error: {e}"),
    }
}
