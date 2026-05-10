//! Boundary-edge split — extracted from `split.rs` to keep
//! that file under the 500-line structural cap.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve::CurveOnEdge;
use crate::geometry::pcurve::Pcurve;
use crate::topology::*;

use super::split::SplitEdgeResult;

impl BRep {
    /// Boundary-edge split: only the surviving half-edge side
    /// exists (twin is null). Creates vertex V, splits he_ab
    /// into he_ab (A→V) + he_vb (V→B), and marks both sub-edges
    /// as boundary.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn split_boundary_edge(
        &mut self,
        edge: EdgeHandle,
        t: f64,
        he_ab: HalfEdgeHandle,
        vertex_a: VertexHandle,
        he_ab_next: HalfEdgeHandle,
        he_ab_loop: LoopHandle,
        pcurve_ab: &Option<Pcurve>,
        curve_on_edge: &Option<CurveOnEdge>,
    ) -> Result<SplitEdgeResult, CoreError> {
        let vertex_b = self.half_edges[he_ab_next].origin;
        let pos_a = self.vertices[vertex_a].position;
        let pos_b = self.vertices[vertex_b].position;

        let split_pos = match curve_on_edge {
            Some(curve) => {
                use crate::geometry::curve::CurveEval;
                let p_start = curve.geometry.point_at(curve.t_start);
                let aligned = pos_a.distance_to(&p_start)
                    <= pos_a.distance_to(&curve.geometry.point_at(curve.t_end));
                let at = if aligned { t } else { 1.0 - t };
                let t_param = curve.t_start + at * (curve.t_end - curve.t_start);
                curve.geometry.point_at(t_param)
            }
            None => Point3::new(
                pos_a.x + t * (pos_b.x - pos_a.x),
                pos_a.y + t * (pos_b.y - pos_a.y),
                pos_a.z + t * (pos_b.z - pos_a.z),
            ),
        };

        let vertex_v = self.vertices.insert(Vertex {
            position: split_pos,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });

        let (pcurve_left, pcurve_right) = super::split::split_pcurve_at(pcurve_ab, t);

        let edge_b = self.edges.insert(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });

        let he_vb = self.half_edges.insert(HalfEdge {
            origin: vertex_v,
            twin: HalfEdgeHandle::default(),
            next: he_ab_next,
            prev: he_ab,
            edge: edge_b,
            loop_: he_ab_loop,
            pcurve: pcurve_right,
        });

        self.half_edges[he_ab].next = he_vb;
        self.half_edges[he_ab_next].prev = he_vb;
        self.half_edges[he_ab].pcurve = pcurve_left;

        self.edges[edge_b].half_edge = he_vb;
        self.vertices[vertex_v].outgoing = he_vb;

        self.boundary_edges.insert(edge_b);

        if let Some(curve) = curve_on_edge {
            use crate::geometry::curve::{CurveEval, CurveGeometry};
            let p_start = curve.geometry.point_at(curve.t_start);
            let aligned = pos_a.distance_to(&p_start)
                <= pos_a.distance_to(&curve.geometry.point_at(curve.t_end));
            let at = if aligned { t } else { 1.0 - t };
            let t_param = curve.t_start + at * (curve.t_end - curve.t_start);
            let nurbs = curve
                .geometry
                .to_nurbs_in_range((curve.t_start, curve.t_end));
            if let Ok((left, right)) = nurbs.split_at(t_param) {
                // Route the two halves to the correct edges — see
                // the matching block in `split_edge` for the
                // rationale. The curve may have been authored in
                // the opposite direction to the half-edge, in
                // which case `left` actually covers the segment
                // near `pos_b` and should go on the new edge, not
                // the original. Pick the assignment that minimises
                // the endpoint mismatch.
                let left_endpoints = {
                    let (t0, t1) = left.parameter_range();
                    (left.point_at(t0), left.point_at(t1))
                };
                let right_endpoints = {
                    let (t0, t1) = right.parameter_range();
                    (right.point_at(t0), right.point_at(t1))
                };
                let dist2 = |a: &Point3, b: &Point3| -> f64 {
                    (a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)
                };
                let cost_natural = dist2(&left_endpoints.0, &pos_a)
                    + dist2(&left_endpoints.1, &split_pos)
                    + dist2(&right_endpoints.0, &split_pos)
                    + dist2(&right_endpoints.1, &pos_b);
                let cost_swapped = dist2(&right_endpoints.0, &pos_a)
                    + dist2(&right_endpoints.1, &split_pos)
                    + dist2(&left_endpoints.0, &split_pos)
                    + dist2(&left_endpoints.1, &pos_b);
                let (old_half, new_half) = if cost_natural <= cost_swapped {
                    (left, right)
                } else {
                    (right, left)
                };

                let lr = old_half.parameter_range();
                let old_geometry = CurveGeometry::Nurbs(old_half);
                let old_geometry = old_geometry
                    .canonical_line_in_range(lr, self.tolerance.linear().max(1.0e-9))
                    .unwrap_or(old_geometry);
                self.edges[edge].curve = Some(CurveOnEdge {
                    geometry: old_geometry,
                    t_start: lr.0,
                    t_end: lr.1,
                });
                let rr = new_half.parameter_range();
                let new_geometry = CurveGeometry::Nurbs(new_half);
                let new_geometry = new_geometry
                    .canonical_line_in_range(rr, self.tolerance.linear().max(1.0e-9))
                    .unwrap_or(new_geometry);
                self.edges[edge_b].curve = Some(CurveOnEdge {
                    geometry: new_geometry,
                    t_start: rr.0,
                    t_end: rr.1,
                });
            }
        }

        Ok(SplitEdgeResult {
            vertex: vertex_v,
            edge_a: edge,
            edge_b,
        })
    }
}
