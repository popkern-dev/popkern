//! What a run hands the page: outcome, diagnostics, report, parameters
//! and per-body triangle buffers. The JSON form (`serde`) carries
//! everything except the bulk buffers, which cross the boundary as
//! typed arrays fetched by body index (see `bindings`).

use serde::Serialize;
use serde_json::Value;

use popkern_rhai::Outcome;

/// Whether the wall-clock budget a host asked for is actually
/// enforceable. On the browser target a clock exists only after the
/// page registers `performance.now`; without one the budget cannot
/// fire, and the envelope says so instead of pretending.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum WallBudget {
    /// No wall-clock budget was requested.
    None,
    /// A budget was requested and a monotonic clock enforces it.
    Armed,
    /// A budget was requested but no monotonic clock is registered:
    /// the operation budget (and the page's worker cancel) are the
    /// only limits in force.
    Unavailable,
}

impl WallBudget {
    /// The state for a host option under the current clock.
    pub fn describe(timeout_requested: bool, clock_available: bool) -> WallBudget {
        match (timeout_requested, clock_available) {
            (false, _) => WallBudget::None,
            (true, true) => WallBudget::Armed,
            (true, false) => WallBudget::Unavailable,
        }
    }
}

/// One body's triangles, ready for a GPU buffer. Positions and normals
/// are interleaved per vertex as `[x, y, z]`; `indices` are triangle
/// triples into them; `face_ids[t]` indexes `face_labels` for triangle
/// `t` so the page can pick and name faces.
#[derive(Debug, Clone, Serialize)]
pub struct BodyMesh {
    /// The body's tree id.
    pub id: u64,
    /// The body's label, when one was given.
    pub label: Option<String>,
    /// Vertex count (`positions.len() / 3`).
    pub vertex_count: usize,
    /// Triangle count (`indices.len() / 3`).
    pub triangle_count: usize,
    /// Face count of the BRep.
    pub face_count: usize,
    /// One spelled entity name per face (`rounds/Generated/0`), or
    /// `"unnamed"` for a face the naming pass did not reach.
    pub face_labels: Vec<String>,
    /// Axis-aligned bounds `[min, max]` of the positions.
    pub bounds: [[f64; 3]; 2],
    /// Open edges in the welded mesh (0 = watertight).
    pub open_edges: usize,
    /// `x y z` per vertex.
    #[serde(skip)]
    pub positions: Vec<f32>,
    /// `nx ny nz` per vertex.
    #[serde(skip)]
    pub normals: Vec<f32>,
    /// Triangle vertex indices.
    #[serde(skip)]
    pub indices: Vec<u32>,
    /// Face table index per triangle.
    #[serde(skip)]
    pub face_ids: Vec<u32>,
    /// Analytic Gaussian curvature per vertex (PK-520 analysis view).
    #[serde(skip)]
    pub gaussian: Vec<f32>,
    /// Analytic mean curvature per vertex (PK-520 analysis view).
    #[serde(skip)]
    pub mean: Vec<f32>,
}

/// One declared parameter as the page's slider panel sees it.
#[derive(Debug, Clone, Serialize)]
pub struct ParamRow {
    /// Qualified name (`wall`, `extrude-1.distance`).
    pub name: String,
    /// Current value in `unit`.
    pub value: f64,
    /// Display unit.
    pub unit: Value,
    /// Driving expression, when any.
    pub expression: Option<String>,
    /// Declared bounds `[min, max]`.
    pub bounds: Option<[f64; 2]>,
    /// Published by a reference; not editable.
    pub read_only: bool,
}

/// The answer to every session verb that changes the model.
#[derive(Debug, Clone, Serialize)]
pub struct Envelope {
    /// How the run ended (flattened: `outcome`, `kind`, `payload`, …).
    #[serde(flatten)]
    pub outcome: Outcome,
    /// Lines the script printed, in order.
    pub printed: Vec<String>,
    /// Wall time of the script run when a clock exists.
    pub elapsed_ms: Option<f64>,
    /// Whether a monotonic clock is available in this process.
    pub clock_available: bool,
    /// Whether the wall-clock budget is enforceable.
    pub wall_budget: WallBudget,
    /// The replay report (spec §11) with `is_clean`, `failures`,
    /// `downgraded_blends` beside it, when a model exists.
    pub report: Option<Value>,
    /// Feature count of the document.
    pub features: usize,
    /// Declared parameters.
    pub params: Vec<ParamRow>,
    /// One triangle buffer set per alive body.
    pub bodies: Vec<BodyMesh>,
}

impl Envelope {
    /// An envelope carrying only an outcome (no model to describe).
    pub fn bare(outcome: Outcome, clock_available: bool, wall_budget: WallBudget) -> Envelope {
        Envelope {
            outcome,
            printed: Vec::new(),
            elapsed_ms: None,
            clock_available,
            wall_budget,
            report: None,
            features: 0,
            params: Vec::new(),
            bodies: Vec::new(),
        }
    }

    /// The JSON form (bulk buffers omitted).
    pub fn to_json(&self) -> Value {
        serde_json::to_value(self).unwrap_or(Value::Null)
    }
}
