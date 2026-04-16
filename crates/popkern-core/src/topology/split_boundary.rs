//! Boundary-edge split — extracted from `split.rs` to keep
//! that file under the 500-line structural cap.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve::CurveOnEdge;
use crate::geometry::nurbs_curve::NurbsCurve;
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
        pcurve_ab: &Option<NurbsCurve>,
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
                let lr = left.parameter_range();
                self.edges[edge].curve = Some(CurveOnEdge {
                    geometry: CurveGeometry::Nurbs(left),
                    t_start: lr.0,
                    t_end: lr.1,
                });
                let rr = right.parameter_range();
                self.edges[edge_b].curve = Some(CurveOnEdge {
                    geometry: CurveGeometry::Nurbs(right),
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
