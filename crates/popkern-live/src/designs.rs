//! Generative designs converging (PK-959): the histories a host holds
//! for its films. A design solved in this process keeps its history; a
//! bake's arrive as bytes (`bakes/<stem>.designs.bin`) and are seeded.
//! A film's design layer asks for the recorded designs and for a mesh
//! at a position among them (two neighbouring records blended), and
//! for the design space it is carved from.

use serde_json::{Value, json};

use popkern_generative::DisplayMesh;

/// Keeps the histories `bytes` hold (a bake's design file). Returns
/// how many.
///
/// # Errors
///
/// The bytes are not a design history file.
pub fn seed(bytes: &[u8]) -> Result<usize, String> {
    popkern_generative::import_histories(bytes)
}

/// What design `design`'s history holds, as the protocol's
/// `design_history` reply fields (`{design, updates, cell, frames:
/// [{update, compliance_nmm, volume_fraction}]}`); `None` when this
/// process holds no history for it.
pub fn summary(design: &str) -> Option<Value> {
    let history = popkern_generative::design_history(design)?;
    Some(json!({
        "design": design,
        "updates": history.updates,
        "cell": history.cell,
        "frames": history.frames.iter().map(|f| json!({
            "update": f.update,
            "compliance_nmm": f.compliance_nmm,
            "volume_fraction": f.volume_fraction,
        })).collect::<Vec<_>>(),
    }))
}

/// Design `design` at position `at` among its records (`None`: the
/// design space), in the design body's units; `None` when this process
/// holds no history for it.
pub fn mesh(design: &str, at: Option<f64>) -> Option<DisplayMesh> {
    let history = popkern_generative::design_history(design)?;
    Some(match at {
        Some(at) => history.frame(at),
        None => history.space.clone(),
    })
}

/// The reply to a design request for `design` that this process holds
/// no history for.
pub fn missing(design: &str) -> Value {
    json!({
        "outcome": "gate_failure",
        "step": "design_history",
        "detail": format!(
            "no history for design {design}: it was solved elsewhere and its bake's design file is not loaded"
        ),
    })
}
