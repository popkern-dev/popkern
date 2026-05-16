//! Scoped read-only BRep export-quality diagnostics.

use crate::topology::{
    BRep, BrepExportQualityOptions, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle,
    VertexHandle,
};

/// Topology handles to inspect for scoped export-quality diagnostics.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct BrepExportQualityScope {
    /// Edges whose edge-local export geometry should be inspected.
    pub edges: Vec<EdgeHandle>,
    /// Vertices adjacent to the scoped topology.
    pub vertices: Vec<VertexHandle>,
    /// Half-edges adjacent to the scoped topology.
    pub half_edges: Vec<HalfEdgeHandle>,
    /// Loops adjacent to the scoped topology.
    pub loops: Vec<LoopHandle>,
    /// Faces adjacent to the scoped topology.
    pub faces: Vec<FaceHandle>,
}

/// Edge-local export-quality report for a scoped topology set.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BrepScopedExportQualityReport {
    /// Number of valid unique scoped edges inspected.
    pub edges: usize,
    /// Number of valid unique scoped vertices.
    pub vertices: usize,
    /// Number of valid unique scoped half-edges.
    pub half_edges: usize,
    /// Number of valid unique scoped loops.
    pub loops: usize,
    /// Number of valid unique scoped faces.
    pub faces: usize,
    /// Number of scoped edges missing 3D curve geometry.
    pub missing_edge_curves: usize,
    /// Number of scoped edges whose pcurve deviation exceeds the requested tolerance.
    pub pcurve_blockers: usize,
    /// Maximum measured scoped edge-to-pcurve deviation.
    pub max_pcurve_deviation: f64,
    /// Number of scoped edge handles no longer present in the BRep.
    pub unavailable_edges: usize,
}

/// Analyze edge-local export quality only inside a caller-supplied scope.
///
/// Whole-model export quality remains authoritative for final readiness. This
/// scoped report is intentionally limited to checks whose ownership is local to
/// the supplied handles, so callers can replace global scans incrementally.
pub fn analyze_scoped_brep_export_quality(
    brep: &BRep,
    scope: &BrepExportQualityScope,
    options: BrepExportQualityOptions,
) -> BrepScopedExportQualityReport {
    let edges = unique_handles(&scope.edges);
    let mut report = BrepScopedExportQualityReport {
        vertices: valid_handle_count(&scope.vertices, |handle| brep.vertices.contains_key(handle)),
        half_edges: valid_handle_count(&scope.half_edges, |handle| {
            brep.half_edges.contains_key(handle)
        }),
        loops: valid_handle_count(&scope.loops, |handle| brep.loops.contains_key(handle)),
        faces: valid_handle_count(&scope.faces, |handle| brep.faces.contains_key(handle)),
        ..Default::default()
    };

    for edge in edges {
        let Some(edge_data) = brep.edges.get(edge) else {
            report.unavailable_edges += 1;
            continue;
        };
        report.edges += 1;
        if edge_data.curve.is_none() {
            report.missing_edge_curves += 1;
        }
        let Some(deviation) = brep.measure_edge_pcurve_deviation(edge) else {
            continue;
        };
        report.max_pcurve_deviation = report.max_pcurve_deviation.max(deviation.max_deviation);
        if let Some(tolerance) = options.max_pcurve_deviation
            && deviation.max_deviation > tolerance
        {
            report.pcurve_blockers += 1;
        }
    }
    report
}

fn valid_handle_count<T: Copy + Eq>(handles: &[T], is_valid: impl Fn(T) -> bool) -> usize {
    unique_handles(handles)
        .into_iter()
        .filter(|handle| is_valid(*handle))
        .count()
}

fn unique_handles<T: Copy + Eq>(handles: &[T]) -> Vec<T> {
    let mut unique = Vec::new();
    for &handle in handles {
        if !unique.contains(&handle) {
            unique.push(handle);
        }
    }
    unique
}
