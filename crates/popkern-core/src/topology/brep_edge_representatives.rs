//! Edge representative and curve-storage normalization for topology rewrites.

use std::collections::BTreeSet;

use super::BRep;
use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use crate::geometry::fitting::fit_curve;
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

    /// Repair edge curves whose stored endpoints no longer match the edge's
    /// topological vertices.
    ///
    /// For straight edges, the line representative is rebuilt exactly between
    /// the owning vertices. For clamped NURBS representatives whose trim range
    /// is the full curve domain, only the first and last control points are
    /// moved, preserving the interior shape and weights while restoring the
    /// endpoint/topology invariant.
    pub fn repair_edge_curve_endpoint_coherence(&mut self, tolerance: f64) -> usize {
        let repairs = self
            .half_edges
            .iter()
            .filter_map(|(_, half_edge)| {
                let curve = self.edges.get(half_edge.edge)?.curve.as_ref()?;
                let start = self.vertices.get(half_edge.origin)?.position;
                let end = self
                    .half_edges
                    .get(half_edge.next)
                    .and_then(|next| self.vertices.get(next.origin))
                    .map(|vertex| vertex.position)?;
                let endpoint_gap = edge_endpoint_gap(curve, start, end);
                (endpoint_gap > tolerance).then_some((half_edge.edge, start, end))
            })
            .collect::<Vec<_>>();

        let mut repaired = 0;
        for (edge_handle, start, end) in repairs {
            let Some(edge) = self.edges.get_mut(edge_handle) else {
                continue;
            };
            let Some(curve) = edge.curve.as_mut() else {
                continue;
            };
            if repair_curve_endpoint_coherence(curve, start, end, tolerance) {
                repaired += 1;
            }
        }
        let loop_closures = self
            .loops
            .iter()
            .filter_map(|(_, loop_data)| {
                let first = self.half_edges.get(loop_data.half_edge)?;
                let last = self.half_edges.get(first.prev)?;
                let start = self.vertices.get(last.origin)?.position;
                let end = self.vertices.get(first.origin)?.position;
                let curve = self.edges.get(last.edge)?.curve.as_ref()?;
                let endpoint_gap = edge_endpoint_gap(curve, start, end);
                (endpoint_gap > tolerance).then_some((last.edge, start, end))
            })
            .collect::<Vec<_>>();
        for (edge_handle, start, end) in loop_closures {
            let Some(edge) = self.edges.get_mut(edge_handle) else {
                continue;
            };
            let Some(curve) = edge.curve.as_mut() else {
                continue;
            };
            if repair_curve_endpoint_coherence(curve, start, end, tolerance) {
                repaired += 1;
            }
        }
        repaired
    }

    /// Give each radial half-edge pair its own edge record when previous
    /// topology surgery left more than one twin pair sharing an edge handle.
    ///
    /// A BRep edge owns exactly one undirected boundary between two vertices.
    /// When several unrelated half-edge pairs share the same edge handle, no
    /// single curve can satisfy every loop-closure endpoint contract. This
    /// pass preserves the first pair on the original edge and clones the edge
    /// representative for each additional pair.
    pub fn separate_nonmanifold_edge_uses(&mut self) -> usize {
        let edge_handles = self.edges.keys().collect::<Vec<_>>();
        let mut separated = 0;
        for edge_handle in edge_handles {
            let half_edges = self
                .half_edges
                .iter()
                .filter_map(|(handle, data)| (data.edge == edge_handle).then_some(handle))
                .collect::<Vec<_>>();
            if half_edges.len() <= 2 {
                continue;
            }
            let Some(template) = self.edges.get(edge_handle).cloned() else {
                continue;
            };
            let mut seen = BTreeSet::new();
            let mut kept_first_pair = false;
            for half_edge in half_edges {
                if !seen.insert(half_edge) {
                    continue;
                }
                let Some(twin) = self.half_edges.get(half_edge).map(|he| he.twin) else {
                    continue;
                };
                seen.insert(twin);
                if kept_first_pair {
                    let mut edge = template.clone();
                    edge.half_edge = half_edge;
                    let new_edge = self.insert_edge_tracked(edge);
                    if let Some(he) = self.half_edges.get_mut(half_edge) {
                        he.edge = new_edge;
                    }
                    if let Some(twin_he) = self.half_edges.get_mut(twin)
                        && twin_he.edge == edge_handle
                    {
                        twin_he.edge = new_edge;
                    }
                    separated += 1;
                } else {
                    if let Some(edge) = self.edges.get_mut(edge_handle) {
                        edge.half_edge = half_edge;
                    }
                    kept_first_pair = true;
                }
            }
        }
        separated
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

fn repair_curve_endpoint_coherence(
    curve: &mut CurveOnEdge,
    start: Point3,
    end: Point3,
    tolerance: f64,
) -> bool {
    let (curve_start, curve_end) = coherent_endpoint_targets(curve, start, end);
    let t_start = curve.t_start;
    let t_end = curve.t_end;
    match &mut curve.geometry {
        CurveGeometry::Line { .. } => {
            let Some(line) = CurveOnEdge::line_segment(curve_start, curve_end, tolerance) else {
                return false;
            };
            *curve = line;
            true
        }
        CurveGeometry::Nurbs(nurbs) => {
            let original = nurbs.clone();
            let (domain_start, domain_end) = nurbs.parameter_range();
            if (curve.t_start - domain_start).abs() <= tolerance
                && (curve.t_end - domain_end).abs() <= tolerance
                && nurbs.control_points.len() >= 2
            {
                let last = nurbs.control_points.len() - 1;
                let previous_first = nurbs.control_points[0];
                let previous_last = nurbs.control_points[last];
                nurbs.control_points[0] = curve_start;
                nurbs.control_points[last] = curve_end;
                if nurbs_endpoint_gap(nurbs, t_start, t_end, start, end) <= tolerance {
                    return true;
                }
                nurbs.control_points[0] = previous_first;
                nurbs.control_points[last] = previous_last;
            }
            let Some(refit) = refit_nurbs_with_endpoint_targets(
                &original,
                t_start,
                t_end,
                curve_start,
                curve_end,
            ) else {
                return repair_curve_as_line_segment(curve, curve_start, curve_end, tolerance);
            };
            let previous = std::mem::replace(nurbs, refit);
            if nurbs_endpoint_gap(nurbs, t_start, t_end, start, end) <= tolerance {
                true
            } else {
                *nurbs = previous;
                repair_curve_as_line_segment(curve, curve_start, curve_end, tolerance)
            }
        }
        _ => false,
    }
}

fn repair_curve_as_line_segment(
    curve: &mut CurveOnEdge,
    curve_start: Point3,
    curve_end: Point3,
    tolerance: f64,
) -> bool {
    let Some(line) = CurveOnEdge::line_segment(curve_start, curve_end, tolerance) else {
        return false;
    };
    *curve = line;
    true
}

fn nurbs_endpoint_gap(
    nurbs: &crate::geometry::NurbsCurve,
    t_start: f64,
    t_end: f64,
    start: Point3,
    end: Point3,
) -> f64 {
    start
        .distance_to(&nurbs.point_at(t_start))
        .max(end.distance_to(&nurbs.point_at(t_end)))
        .min(
            start
                .distance_to(&nurbs.point_at(t_end))
                .max(end.distance_to(&nurbs.point_at(t_start))),
        )
}

fn refit_nurbs_with_endpoint_targets(
    original: &crate::geometry::NurbsCurve,
    t_start: f64,
    t_end: f64,
    curve_start: Point3,
    curve_end: Point3,
) -> Option<crate::geometry::NurbsCurve> {
    let samples = original.control_points.len().max(8);
    let mut points = Vec::with_capacity(samples);
    let mut params = Vec::with_capacity(samples);
    for index in 0..samples {
        let s = index as f64 / (samples - 1) as f64;
        let t = t_start + (t_end - t_start) * s;
        let point = if index == 0 {
            curve_start
        } else if index + 1 == samples {
            curve_end
        } else {
            original.point_at(t)
        };
        points.push(point);
        params.push(t);
    }
    let degree = original.degree.min(samples.saturating_sub(1)).max(1);
    let control_points = original.control_points.len().max(degree + 1);
    fit_curve(&points, &params, degree, control_points).ok()
}

fn coherent_endpoint_targets(curve: &CurveOnEdge, start: Point3, end: Point3) -> (Point3, Point3) {
    let curve_start = curve.geometry.point_at(curve.t_start);
    let curve_end = curve.geometry.point_at(curve.t_end);
    let natural_gap = start
        .distance_to(&curve_start)
        .max(end.distance_to(&curve_end));
    let swapped_gap = start
        .distance_to(&curve_end)
        .max(end.distance_to(&curve_start));
    if natural_gap <= swapped_gap {
        (start, end)
    } else {
        (end, start)
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
