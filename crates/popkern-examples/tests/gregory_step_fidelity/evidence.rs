//! Fidelity-state evidence records shared by the Gregory STEP
//! fidelity diagnostics.

#[derive(Debug, Clone, PartialEq)]
pub struct FidelityState {
    pub label: &'static str,
    pub edges: usize,
    pub faces: usize,
    pub trim_findings: usize,
    pub pcurve_blockers: usize,
    pub max_pcurve_deviation: f64,
    pub edge_surface_blockers: usize,
    pub max_edge_surface_deviation: f64,
    pub worst_pcurve_edges: Vec<EdgePcurveEvidence>,
    pub worst_edge_surface_edges: Vec<EdgeSurfaceEvidence>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct EdgePcurveEvidence {
    pub edge: String,
    pub max_deviation: f64,
    pub half_edges: Vec<HalfEdgePcurveEvidence>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct HalfEdgePcurveEvidence {
    pub half_edge: String,
    pub face: String,
    pub surface_kind: &'static str,
    pub max_deviation: f64,
    pub same_fraction_deviation: f64,
    pub trace_deviation: f64,
    pub worst_fraction: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EdgeSurfaceEvidence {
    pub edge: String,
    pub half_edge: String,
    pub face: String,
    pub surface_kind: &'static str,
    pub max_deviation: f64,
    pub max_curve_parameter: f64,
}
