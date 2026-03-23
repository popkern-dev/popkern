//! Face and edge splitting operations for surface trimming.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::tolerance::Tolerance;
use crate::topology::*;

/// Result of splitting a face along a trim curve.
#[derive(Debug, Clone)]
pub struct SplitFaceResult {
    /// One side of the split (the original face).
    pub face_a: FaceHandle,
    /// The other side of the split (the new face from mef).
    pub face_b: FaceHandle,
    /// The newly created edges along the trim curve.
    pub new_edges: Vec<EdgeHandle>,
}

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
    pub fn split_edge(&mut self, edge: EdgeHandle, t: f64) -> Result<SplitEdgeResult, CoreError> {
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
            pcurve: None,
        });

        let he_bv = self.half_edges.insert(HalfEdge {
            origin: vertex_b,
            twin: he_vb,
            next: he_ba,
            prev: he_ba_prev,
            edge: edge_b,
            loop_: he_ba_loop,
            pcurve: None,
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
            let nurbs = curve
                .geometry
                .to_nurbs_in_range((curve.t_start, curve.t_end));
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

    /// Split a face along a trim curve.
    ///
    /// The trim curve endpoints must lie at existing vertices on the face's
    /// boundary (created by prior `split_edge` calls if needed). This method
    /// finds the entry and exit vertices, calls `mef` to split the topology,
    /// then attaches the trim curve geometry to the new edge.
    ///
    /// Both resulting faces retain the original face's surface geometry.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if the face handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if the trim curve endpoints cannot
    ///   be matched to vertices on the face's outer loop within tolerance.
    pub fn split_face(
        &mut self,
        face: FaceHandle,
        _trim_curve_uv: &NurbsCurve,
        trim_curve_3d: &NurbsCurve,
        tolerance: &Tolerance,
    ) -> Result<SplitFaceResult, CoreError> {
        // Validate face handle.
        let face_data = self
            .faces
            .get(face)
            .ok_or_else(|| CoreError::InvalidHandle("face handle is invalid".into()))?;

        let outer_loop = face_data.outer_loop;

        // Get the trim curve start and end points in 3D.
        let (t_start, t_end) = trim_curve_3d.parameter_range();
        let curve_start = trim_curve_3d.point_at(t_start);
        let curve_end = trim_curve_3d.point_at(t_end);

        // Walk the outer loop to find half-edges whose origin vertices
        // match the curve start and end points.
        let loop_start_he = self
            .loops
            .get(outer_loop)
            .ok_or_else(|| CoreError::InvalidHandle("outer loop handle is invalid".into()))?
            .half_edge;

        let mut he_entry: Option<HalfEdgeHandle> = None;
        let mut he_exit: Option<HalfEdgeHandle> = None;

        let mut cursor = loop_start_he;
        let max_iters = self.half_edges.len() + 1;
        let mut iters = 0;
        loop {
            let he_data = self.half_edges.get(cursor).ok_or_else(|| {
                CoreError::InvalidHandle("half_edge in loop walk is invalid".into())
            })?;
            let vertex_pos = self.vertices[he_data.origin].position;

            if he_entry.is_none() && tolerance.points_coincident(&vertex_pos, &curve_start) {
                he_entry = Some(cursor);
            } else if he_exit.is_none() && tolerance.points_coincident(&vertex_pos, &curve_end) {
                he_exit = Some(cursor);
            }

            cursor = he_data.next;
            iters += 1;
            if cursor == loop_start_he || iters >= max_iters {
                break;
            }
        }

        // If we found exit but not entry on first pass, or entry but not exit,
        // do a second pass checking the other direction assignment. The curve
        // endpoints may appear in either order around the loop.
        if he_entry.is_some() && he_exit.is_none() {
            // Try matching curve_end to remaining vertices
            let mut cursor2 = loop_start_he;
            let mut iters2 = 0;
            loop {
                let he_data = &self.half_edges[cursor2];
                let vertex_pos = self.vertices[he_data.origin].position;
                if cursor2 != he_entry.unwrap_or_default()
                    && tolerance.points_coincident(&vertex_pos, &curve_end)
                {
                    he_exit = Some(cursor2);
                    break;
                }
                cursor2 = he_data.next;
                iters2 += 1;
                if cursor2 == loop_start_he || iters2 >= max_iters {
                    break;
                }
            }
        } else if he_exit.is_some() && he_entry.is_none() {
            // Try matching curve_start to remaining vertices
            let mut cursor2 = loop_start_he;
            let mut iters2 = 0;
            loop {
                let he_data = &self.half_edges[cursor2];
                let vertex_pos = self.vertices[he_data.origin].position;
                if cursor2 != he_exit.unwrap_or_default()
                    && tolerance.points_coincident(&vertex_pos, &curve_start)
                {
                    he_entry = Some(cursor2);
                    break;
                }
                cursor2 = he_data.next;
                iters2 += 1;
                if cursor2 == loop_start_he || iters2 >= max_iters {
                    break;
                }
            }
        }

        let he_entry = he_entry.ok_or_else(|| {
            CoreError::PreconditionFailed(
                "trim curve start point does not match any vertex on the face boundary".into(),
            )
        })?;
        let he_exit = he_exit.ok_or_else(|| {
            CoreError::PreconditionFailed(
                "trim curve end point does not match any vertex on the face boundary".into(),
            )
        })?;

        // Capture the original surface before mef (which doesn't copy it).
        let original_surface = self.faces[face].surface.clone();

        // Call mef to split the face topology.
        let mef_result = self.mef(he_entry, he_exit)?;

        // Attach the trim curve geometry to the new edge.
        self.edges[mef_result.edge].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(trim_curve_3d.clone()),
            t_start: 0.0,
            t_end: 1.0,
        });

        // Both faces should keep the original surface geometry.
        // The original face (containing he_entry's side) already has it.
        // The new face needs a copy.
        if let Some(ref surface) = original_surface {
            self.faces[face].surface = Some(surface.clone());
            self.faces[mef_result.face].surface = Some(surface.clone());
        }

        // Determine which face is face_a (original) and face_b (new).
        // After mef, the old loop stays on the old face, and the new loop
        // is on the new face. The original face handle is still `face`.
        Ok(SplitFaceResult {
            face_a: face,
            face_b: mef_result.face,
            new_edges: vec![mef_result.edge],
        })
    }
}

#[cfg(test)]
#[path = "split_tests.rs"]
mod tests;
