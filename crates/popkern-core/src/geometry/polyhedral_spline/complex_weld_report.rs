use super::complex_weld_scope::PolyhedralBoundaryWeldScope;
use crate::topology::{
    BRep, BrepEdgeCurveCompletionDeclineReason, BrepRepairActionArtifact,
    BrepScopedExportQualityReport, BrepSewingClassificationBucket,
};

#[derive(Debug, Clone)]
pub struct PolyhedralWeldedSurfaceComplex {
    pub brep: BRep,
    pub report: PolyhedralBoundaryWeldReport,
    pub repair_actions: BrepRepairActionArtifact,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PolyhedralBoundaryWeldReport {
    pub candidate_pairs: usize,
    pub welded_pairs: usize,
    pub declined_pairs: usize,
    pub remaining_boundary_edges: usize,
    pub max_weld_gap: f64,
    pub trim_repairs: usize,
    pub edge_curve_repairs: usize,
    pub edge_curve_declines: usize,
    pub edge_curve_decline_reasons: Vec<(BrepEdgeCurveCompletionDeclineReason, usize)>,
    pub classification_buckets: Vec<(BrepSewingClassificationBucket, usize)>,
    pub degenerate_edges: usize,
    pub collapsed_surface_trims: usize,
    pub zero_length_vertex_edges: usize,
    pub ambiguous_degenerate_edges: usize,
    pub topology_errors: usize,
    pub topology_warnings: usize,
    pub export_ready: bool,
    pub expected_open_sheet_boundaries: bool,
    pub export_blockers: usize,
    pub missing_edge_curves: usize,
    pub trim_blockers: usize,
    pub scope: PolyhedralBoundaryWeldScope,
    pub pre_repair_scoped_quality: BrepScopedExportQualityReport,
    pub post_repair_scoped_quality: BrepScopedExportQualityReport,
    pub timing: PolyhedralBoundaryWeldTiming,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PolyhedralBoundaryWeldTiming {
    pub assemble_ms: u128,
    pub candidate_match_ms: u128,
    pub repair_pipeline_ms: u128,
    pub export_quality_ms: u128,
    pub repair_artifact_ms: u128,
    pub total_ms: u128,
}
