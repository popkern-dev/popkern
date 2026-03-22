//! Face and edge splitting operations for surface trimming.

use crate::error::CoreError;
use crate::geometry::curve::CurveEval;
use crate::geometry::Point3;
use crate::topology::*;

/// Result of splitting an edge at a parameter.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SplitEdgeResult {
    /// The newly created vertex at the split point.
    pub vertex: VertexHandle,
    /// The original edge (start -> new vertex).
    pub edge_a: EdgeHandle,
    /// The new edge (new vertex -> end).
    pub edge_b: EdgeHandle,
}

impl BRep {
    /// Split an edge at parameter `t` (0-1), inserting a new vertex.
    ///
    /// The original edge is shortened to cover [start, split_point] and a new
    /// edge is created for [split_point, end]. All half-edge chains and loop
    /// references are kept consistent.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if the edge handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if `t` is not in (0, 1).
    pub fn split_edge(
        &mut self,
        edge: EdgeHandle,
        t: f64,
    ) -> Result<SplitEdgeResult, CoreError> {
        // Validate parameter range.
        if t <= 0.0 || t >= 1.0 {
            return Err(CoreError::PreconditionFailed(format!(
                "split parameter t={t} must be in (0, 1)"
            )));
        }

        // Look up the edge.
        let edge_data = self
            .edges
            .get(edge)
            .ok_or_else(|| CoreError::InvalidHandle("edge handle is invalid".into()))?;

        let he_ab = edge_data.half_edge;
        let curve_on_edge = edge_data.curve.clone();

        // Get both half-edges.
        let he_ab_data = self
            .half_edges
            .get(he_ab)
            .ok_or_else(|| CoreError::InvalidHandle("half_edge on edge is invalid".into()))?;

        let he_ba = he_ab_data.twin;
        let vertex_a = he_ab_data.origin;
        let he_ab_next = he_ab_data.next;
        let he_ab_loop = he_ab_data.loop_;

        let he_ba_data = self
            .half_edges
            .get(he_ba)
            .ok_or_else(|| CoreError::InvalidHandle("twin half_edge is invalid".into()))?;

        let vertex_b = he_ba_data.origin;
        let he_ba_prev = he_ba_data.prev;
        let he_ba_loop = he_ba_data.loop_;

        // Compute the position of the new vertex.
        let pos_a = self.vertices[vertex_a].position;
        let pos_b = self.vertices[vertex_b].position;

        let split_pos = match &curve_on_edge {
            Some(curve) => {
                let t_param = curve.t_start + t * (curve.t_end - curve.t_start);
                curve.geometry.point_at(t_param)
            }
            None => Point3::new(
                pos_a.x + t * (pos_b.x - pos_a.x),
                pos_a.y + t * (pos_b.y - pos_a.y),
                pos_a.z + t * (pos_b.z - pos_a.z),
            ),
        };

        // Create the new vertex V.
        let vertex_v = self.vertices.insert(Vertex {
            position: split_pos,
            outgoing: Default::default(),
        });

        // Create the new edge (edge_b).
        let edge_b = self.edges.insert(Edge {
            half_edge: Default::default(),
            curve: None,
        });

        // Create new half-edges: he_vb (V->B) and he_bv (B->V).
        let he_vb = self.half_edges.insert(HalfEdge {
            origin: vertex_v,
            twin: Default::default(),
            next: he_ab_next,
            prev: he_ab,
            edge: edge_b,
            loop_: he_ab_loop,
        });

        let he_bv = self.half_edges.insert(HalfEdge {
            origin: vertex_b,
            twin: he_vb,
            next: he_ba,
            prev: he_ba_prev,
            edge: edge_b,
            loop_: he_ba_loop,
        });

        // Patch twin of he_vb.
        self.half_edges[he_vb].twin = he_bv;

        // Rewire the chain on he_ab's side:
        // he_ab -> he_vb -> he_ab_next  (was: he_ab -> he_ab_next)
        self.half_edges[he_ab].next = he_vb;
        self.half_edges[he_ab_next].prev = he_vb;

        // Rewire the chain on he_ba's side:
        // he_ba_prev -> he_bv -> he_ba  (was: he_ba_prev -> he_ba)
        self.half_edges[he_ba_prev].next = he_bv;
        self.half_edges[he_ba].prev = he_bv;

        // Update he_ba origin: was B, now V (since the original edge A->B
        // becomes A->V, its twin goes V->A).
        self.half_edges[he_ba].origin = vertex_v;

        // Update edge_b to point to he_vb.
        self.edges[edge_b].half_edge = he_vb;

        // Update vertex outgoing pointers.
        self.vertices[vertex_v].outgoing = he_vb;

        // If B's outgoing was he_ba (which now originates at V), update it to he_bv.
        if self.vertices[vertex_b].outgoing == he_ba {
            self.vertices[vertex_b].outgoing = he_bv;
        }

        // Split curve geometry if present.
        if let Some(curve) = curve_on_edge {
            let t_param = curve.t_start + t * (curve.t_end - curve.t_start);

            // Try to split NURBS; if the curve is analytic, convert to NURBS first.
            let nurbs = curve.geometry.to_nurbs_in_range((curve.t_start, curve.t_end));
            match nurbs.split_at(t_param) {
                Ok((left, right)) => {
                    use crate::geometry::curve::{CurveGeometry, CurveOnEdge};

                    let left_range = left.parameter_range();
                    self.edges[edge].curve = Some(CurveOnEdge {
                        geometry: CurveGeometry::Nurbs(left),
                        t_start: left_range.0,
                        t_end: left_range.1,
                    });

                    let right_range = right.parameter_range();
                    self.edges[edge_b].curve = Some(CurveOnEdge {
                        geometry: CurveGeometry::Nurbs(right),
                        t_start: right_range.0,
                        t_end: right_range.1,
                    });
                }
                Err(_) => {
                    // If split fails (degenerate), keep original on edge_a,
                    // leave edge_b without a curve.
                }
            }
        }

        Ok(SplitEdgeResult {
            vertex: vertex_v,
            edge_a: edge,
            edge_b,
        })
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn split_edge_simple() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();

        let result = brep.split_edge(mev1.edge, 0.5).unwrap();

        assert_eq!(brep.vertices.len(), 3);
        assert_eq!(brep.edges.len(), 2);
        assert_eq!(brep.half_edges.len(), 4);
        let v = &brep.vertices[result.vertex];
        assert!((v.position.x - 0.5).abs() < 1e-10);
    }

    #[test]
    fn split_edge_on_cube() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let v_before = brep.vertices.len();
        let e_before = brep.edges.len();
        let edge = brep.edges.keys().next().unwrap();

        brep.split_edge(edge, 0.5).unwrap();

        assert_eq!(brep.vertices.len(), v_before + 1);
        assert_eq!(brep.edges.len(), e_before + 1);
    }

    #[test]
    fn split_edge_invalid_handle() {
        let mut brep = BRep::new(Tolerance::default());
        assert!(brep.split_edge(Default::default(), 0.5).is_err());
    }

    #[test]
    fn split_edge_half_edge_chain_valid() {
        // After split, walking the loop should still form a cycle
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let _mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();

        brep.split_edge(mev1.edge, 0.5).unwrap();

        // Walk the loop -- should return to start
        let start = brep.loops[mvfs.loop_].half_edge;
        let mut he = start;
        let mut count = 0;
        loop {
            he = brep.half_edges[he].next;
            count += 1;
            if he == start || count > 100 {
                break;
            }
        }
        assert_eq!(he, start, "Loop should be closed after split_edge");
    }

    #[test]
    fn split_edge_t_out_of_range() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();

        assert!(brep.split_edge(mev1.edge, 0.0).is_err());
        assert!(brep.split_edge(mev1.edge, 1.0).is_err());
        assert!(brep.split_edge(mev1.edge, -0.1).is_err());
        assert!(brep.split_edge(mev1.edge, 1.1).is_err());
    }
}
