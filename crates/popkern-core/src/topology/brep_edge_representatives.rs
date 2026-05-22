//! Edge representative and curve-storage normalization for topology rewrites.

use super::BRep;
use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use crate::topology::EdgeHandle;

impl BRep {
    /// Restore the invariant that each edge's representative half-edge points
    /// back to that same edge, and remove edge records with no half-edges.
    pub fn normalize_edge_representatives(&mut self) {
        let mut representative_by_edge = std::collections::BTreeMap::new();
        for (half_edge_handle, half_edge) in &self.half_edges {
            representative_by_edge
                .entry(half_edge.edge)
                .or_insert(half_edge_handle);
        }

        let edge_handles = self.edges.keys().collect::<Vec<_>>();
        for edge_handle in edge_handles {
            let Some(representative) = representative_by_edge.get(&edge_handle).copied() else {
                self.edges.remove(edge_handle);
                self.boundary_edges.remove(&edge_handle);
                continue;
            };
            let primary_is_valid = self
                .edges
                .get(edge_handle)
                .and_then(|edge| self.half_edges.get(edge.half_edge))
                .is_some_and(|half_edge| half_edge.edge == edge_handle);
            if !primary_is_valid && let Some(edge) = self.edges.get_mut(edge_handle) {
                edge.half_edge = representative;
            }
        }
    }

    /// Canonicalize edge curves that are geometrically straight over their
    /// stored trim range to [`CurveGeometry::Line`].
    ///
    /// Construction stages can temporarily represent straight trim edges as
    /// degree-1 or fitted NURBS while repairing topology. This pass restores
    /// the semantic storage invariant without changing edge endpoints,
    /// parameter ranges, or half-edge ownership.
    pub fn normalize_semantic_line_edge_curves(&mut self, tolerance: f64) -> usize {
        let mut normalized = 0;
        for (_, edge) in self.edges.iter_mut() {
            let Some(curve) = edge.curve.as_mut() else {
                continue;
            };
            if matches!(curve.geometry, CurveGeometry::Line { .. }) {
                continue;
            }
            if let Some(line) = curve
                .geometry
                .canonical_line_in_range((curve.t_start, curve.t_end), tolerance)
            {
                curve.geometry = line;
                normalized += 1;
            }
        }
        normalized
    }

    /// Canonicalize geometrically straight edge curves only when existing
    /// pcurves still reproduce the candidate line within the same tolerance.
    ///
    /// Exact-boundary and imported trim edges can carry fitted 3-D curves that
    /// are visually close to straight but are semantically coupled to authored
    /// pcurves. Blind line canonicalization breaks those triples. This guarded
    /// variant preserves the shared edge/pcurve/surface invariant.
    pub fn normalize_semantic_line_edge_curves_preserving_pcurves(
        &mut self,
        tolerance: f64,
    ) -> usize {
        let candidates = self
            .edges
            .iter()
            .filter_map(|(edge_handle, edge)| {
                let curve = edge.curve.as_ref()?;
                if matches!(curve.geometry, CurveGeometry::Line { .. }) {
                    return None;
                }
                let line = curve
                    .geometry
                    .canonical_line_in_range((curve.t_start, curve.t_end), tolerance)?;
                Some((edge_handle, line))
            })
            .collect::<Vec<(EdgeHandle, CurveGeometry)>>();

        let mut normalized = 0;
        for (edge_handle, line) in candidates {
            let Some(previous) = self
                .edges
                .get(edge_handle)
                .and_then(|edge| edge.curve.as_ref())
                .map(|curve| curve.geometry.clone())
            else {
                continue;
            };
            if let Some(curve) = self
                .edges
                .get_mut(edge_handle)
                .and_then(|edge| edge.curve.as_mut())
            {
                curve.geometry = line;
            }
            let preserves_pcurves = self
                .measure_edge_pcurve_deviation(edge_handle)
                .map(|report| report.max_deviation <= tolerance)
                .unwrap_or(true);
            if preserves_pcurves {
                normalized += 1;
            } else if let Some(curve) = self
                .edges
                .get_mut(edge_handle)
                .and_then(|edge| edge.curve.as_mut())
            {
                curve.geometry = previous;
            }
        }
        normalized
    }

    /// Repair straight edge curves whose stored endpoints no longer match the
    /// edge's topological vertices.
    ///
    /// Topology rewrites can retarget vertices while leaving an analytic line
    /// representative on the old support span. For straight edges, the
    /// commercial-grade invariant is unambiguous: the bounded edge curve must
    /// terminate at the two vertices owned by the edge's twin half-edges.
    pub fn repair_line_edge_curve_endpoint_coherence(&mut self, tolerance: f64) -> usize {
        let repairs = self
            .edges
            .iter()
            .filter_map(|(edge_handle, edge)| {
                let curve = edge.curve.as_ref()?;
                if !matches!(curve.geometry, CurveGeometry::Line { .. }) {
                    return None;
                }
                let (start, end) = edge_vertices(self, edge.half_edge)?;
                let endpoint_gap = edge_endpoint_gap(curve, start, end);
                (endpoint_gap > tolerance).then_some((edge_handle, start, end))
            })
            .collect::<Vec<_>>();

        let mut repaired = 0;
        for (edge_handle, start, end) in repairs {
            let Some(curve) = CurveOnEdge::line_segment(start, end, tolerance) else {
                continue;
            };
            if let Some(edge) = self.edges.get_mut(edge_handle) {
                edge.curve = Some(curve);
                repaired += 1;
            }
        }
        repaired
    }
}

fn edge_vertices(
    brep: &BRep,
    half_edge: crate::topology::HalfEdgeHandle,
) -> Option<(Point3, Point3)> {
    let half_edge = brep.half_edges.get(half_edge)?;
    let twin = brep.half_edges.get(half_edge.twin)?;
    Some((
        brep.vertices.get(half_edge.origin)?.position,
        brep.vertices.get(twin.origin)?.position,
    ))
}

fn edge_endpoint_gap(curve: &CurveOnEdge, start: Point3, end: Point3) -> f64 {
    let curve_start = curve.geometry.point_at(curve.t_start);
    let curve_end = curve.geometry.point_at(curve.t_end);
    start
        .distance_to(&curve_start)
        .max(end.distance_to(&curve_end))
        .min(
            start
                .distance_to(&curve_end)
                .max(end.distance_to(&curve_start)),
        )
}
