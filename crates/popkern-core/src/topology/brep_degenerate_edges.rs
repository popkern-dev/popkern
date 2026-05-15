//! Read-only degenerate edge diagnostics for BReps.

use crate::geometry::{Point3, SurfaceEval};
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle, edge_segment};

/// Options controlling degenerate edge diagnostics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BrepDegenerateEdgeOptions {
    /// Local geometric tolerance for endpoint and sampled span checks.
    pub tolerance: f64,
    /// Number of samples used to measure pcurve/surface image span.
    pub samples: usize,
}

/// Aggregate degenerate edge diagnostic report.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepDegenerateEdgeReport {
    /// Number of edges inspected.
    pub inspected_edges: usize,
    /// Number of edges whose endpoints are coincident within tolerance.
    pub degenerate_edges: usize,
    /// Number of degenerate edges classified as valid collapsed surface trims.
    pub collapsed_surface_trims: usize,
    /// Number of degenerate edges with no surface-trim evidence.
    pub zero_length_vertex_edges: usize,
    /// Number of degenerate edges with conflicting or non-collapsed trim evidence.
    pub ambiguous_edges: usize,
    /// Maximum endpoint gap among degenerate edges.
    pub max_endpoint_gap: f64,
    /// Maximum sampled surface-image span among degenerate edges.
    pub max_surface_span: f64,
    /// Per-edge findings.
    pub findings: Vec<BrepDegenerateEdgeFinding>,
}

/// Per-edge degenerate diagnostic finding.
#[derive(Debug, Clone, PartialEq)]
pub struct BrepDegenerateEdgeFinding {
    /// Edge being classified.
    pub edge: EdgeHandle,
    /// Classification kind.
    pub kind: BrepDegenerateEdgeKind,
    /// Distance between edge endpoint vertices.
    pub endpoint_gap: f64,
    /// Maximum sampled span of owning pcurve/surface images.
    pub surface_span: f64,
    /// Maximum UV endpoint span among owning pcurves.
    pub pcurve_span: f64,
    /// Number of pcurve/surface sources used.
    pub source_trims: usize,
}

/// Machine-readable degenerate edge classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepDegenerateEdgeKind {
    /// Edge has coincident 3D endpoints and a collapsed pcurve/surface image.
    CollapsedSurfaceTrim,
    /// Edge has coincident endpoints and no surface-trim evidence.
    ZeroLengthVertexEdge,
    /// Edge has coincident endpoints but non-collapsed or conflicting trim evidence.
    AmbiguousCollapsedTrim,
    /// Edge topology was incomplete during classification.
    MissingTopology,
}

struct SurfaceTrimSource {
    half_edge: HalfEdgeHandle,
}

/// Classify degenerate edges without mutating the BRep.
pub fn analyze_brep_degenerate_edges(
    brep: &BRep,
    options: BrepDegenerateEdgeOptions,
) -> BrepDegenerateEdgeReport {
    let mut report = BrepDegenerateEdgeReport {
        inspected_edges: brep.edges.len(),
        degenerate_edges: 0,
        collapsed_surface_trims: 0,
        zero_length_vertex_edges: 0,
        ambiguous_edges: 0,
        max_endpoint_gap: 0.0,
        max_surface_span: 0.0,
        findings: Vec::new(),
    };

    for edge in brep.edges.keys() {
        let finding = classify_edge(brep, edge, options);
        if let Some(finding) = finding {
            report.degenerate_edges += 1;
            report.max_endpoint_gap = report.max_endpoint_gap.max(finding.endpoint_gap);
            report.max_surface_span = report.max_surface_span.max(finding.surface_span);
            match finding.kind {
                BrepDegenerateEdgeKind::CollapsedSurfaceTrim => {
                    report.collapsed_surface_trims += 1;
                }
                BrepDegenerateEdgeKind::ZeroLengthVertexEdge => {
                    report.zero_length_vertex_edges += 1;
                }
                BrepDegenerateEdgeKind::AmbiguousCollapsedTrim
                | BrepDegenerateEdgeKind::MissingTopology => {
                    report.ambiguous_edges += 1;
                }
            }
            report.findings.push(finding);
        }
    }
    report
}

fn classify_edge(
    brep: &BRep,
    edge: EdgeHandle,
    options: BrepDegenerateEdgeOptions,
) -> Option<BrepDegenerateEdgeFinding> {
    let (start, end) = match edge_segment(brep, edge) {
        Some(segment) => segment,
        None => {
            return Some(BrepDegenerateEdgeFinding {
                edge,
                kind: BrepDegenerateEdgeKind::MissingTopology,
                endpoint_gap: f64::INFINITY,
                surface_span: f64::INFINITY,
                pcurve_span: f64::INFINITY,
                source_trims: 0,
            });
        }
    };
    let endpoint_gap = start.distance_to(&end);
    if endpoint_gap > options.tolerance {
        return None;
    }

    let sources = surface_trim_sources(brep, edge);
    let mut max_surface_span = 0.0_f64;
    let mut max_pcurve_span = 0.0_f64;
    for source in &sources {
        max_surface_span = max_surface_span
            .max(surface_trim_span(brep, source, options.samples).unwrap_or(f64::INFINITY));
        max_pcurve_span =
            max_pcurve_span.max(pcurve_endpoint_span(brep, source).unwrap_or(f64::INFINITY));
    }

    let kind = if sources.is_empty() {
        BrepDegenerateEdgeKind::ZeroLengthVertexEdge
    } else if max_surface_span <= options.tolerance {
        BrepDegenerateEdgeKind::CollapsedSurfaceTrim
    } else {
        BrepDegenerateEdgeKind::AmbiguousCollapsedTrim
    };

    Some(BrepDegenerateEdgeFinding {
        edge,
        kind,
        endpoint_gap,
        surface_span: max_surface_span,
        pcurve_span: max_pcurve_span,
        source_trims: sources.len(),
    })
}

fn surface_trim_sources(brep: &BRep, edge: EdgeHandle) -> Vec<SurfaceTrimSource> {
    brep.half_edges
        .iter()
        .filter_map(|(half_edge, data)| {
            if data.edge != edge || data.pcurve.is_none() {
                return None;
            }
            let face = brep.faces.get(brep.loops.get(data.loop_)?.face)?;
            face.surface
                .is_some()
                .then_some(SurfaceTrimSource { half_edge })
        })
        .collect()
}

fn surface_trim_span(brep: &BRep, source: &SurfaceTrimSource, samples: usize) -> Option<f64> {
    let points = sample_surface_trim(brep, source, samples)?;
    let first = *points.first()?;
    Some(
        points
            .iter()
            .map(|point| point.distance_to(&first))
            .fold(0.0_f64, f64::max),
    )
}

fn sample_surface_trim(
    brep: &BRep,
    source: &SurfaceTrimSource,
    samples: usize,
) -> Option<Vec<Point3>> {
    let half_edge = brep.half_edges.get(source.half_edge)?;
    let pcurve = half_edge.pcurve.as_ref()?;
    let face = brep.faces.get(brep.loops.get(half_edge.loop_)?.face)?;
    let surface = &face.surface.as_ref()?.geometry;
    let samples = samples.max(2);
    let mut points = Vec::with_capacity(samples);
    for i in 0..samples {
        let t = pcurve.t_start + (pcurve.t_end - pcurve.t_start) * i as f64 / (samples - 1) as f64;
        let uv = pcurve.point_at(t);
        points.push(surface.point_at(uv.x, uv.y));
    }
    Some(points)
}

fn pcurve_endpoint_span(brep: &BRep, source: &SurfaceTrimSource) -> Option<f64> {
    let pcurve = brep.half_edges.get(source.half_edge)?.pcurve.as_ref()?;
    let start = pcurve.point_at(pcurve.t_start);
    let end = pcurve.point_at(pcurve.t_end);
    Some(start.distance_to(&end))
}

#[cfg(test)]
#[path = "brep_degenerate_edges_tests.rs"]
mod tests;
