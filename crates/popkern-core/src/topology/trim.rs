//! Inner loop (hole) creation for face trimming.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::tolerance::Tolerance;
use crate::topology::*;

/// Result of adding an inner trim loop to a face.
#[derive(Debug, Clone)]
pub struct AddInnerLoopResult {
    /// The new inner loop.
    pub loop_: LoopHandle,
    /// Vertices created along the curve.
    pub vertices: Vec<VertexHandle>,
    /// Edges created between consecutive vertices.
    pub edges: Vec<EdgeHandle>,
}

impl BRep {
    /// Add an inner loop (hole) to a face along a closed curve.
    ///
    /// Samples the curve at `num_segments` evenly-spaced points, creates
    /// vertices, edges, and half-edges forming a closed inner loop.
    /// Each half-edge gets a pcurve from the UV curve.
    ///
    /// The curve must be closed (or nearly closed within tolerance).
    /// The twin half-edges are self-twinned and belong to no loop; they
    /// will be stitched during the boolean merge stage.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if the face handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if `num_segments < 3` or the curve
    ///   is not closed within tolerance.
    pub fn add_inner_loop(
        &mut self,
        face: FaceHandle,
        curve_3d: &NurbsCurve,
        curve_uv: &NurbsCurve,
        num_segments: usize,
        tolerance: &Tolerance,
    ) -> Result<AddInnerLoopResult, CoreError> {
        // Validate face handle.
        if !self.faces.contains_key(face) {
            return Err(CoreError::InvalidHandle("face handle is invalid".into()));
        }

        // Validate segment count.
        if num_segments < 3 {
            return Err(CoreError::PreconditionFailed(format!(
                "num_segments must be >= 3, got {num_segments}"
            )));
        }

        // Validate curve is closed.
        let (t_start, t_end) = curve_3d.parameter_range();
        let start_pt = curve_3d.point_at(t_start);
        let end_pt = curve_3d.point_at(t_end);
        if !tolerance.points_coincident(&start_pt, &end_pt) {
            return Err(CoreError::PreconditionFailed(format!(
                "curve_3d is not closed: gap = {}",
                start_pt.distance_to(&end_pt)
            )));
        }

        // Sample parameters evenly (excluding the duplicate endpoint).
        let params: Vec<f64> = (0..num_segments)
            .map(|i| t_start + (t_end - t_start) * (i as f64) / (num_segments as f64))
            .collect();

        // Sample 3D positions and create vertices.
        let mut vertices = Vec::with_capacity(num_segments);
        for &t in &params {
            let pos = curve_3d.point_at(t);
            let vh = self.vertices.insert(Vertex {
                position: pos,
                outgoing: Default::default(),
                name: None,
                tolerance: None,
            });
            vertices.push(vh);
        }

        // Sample UV positions for pcurves.
        let uv_positions: Vec<Point3> = params.iter().map(|&t| curve_uv.point_at(t)).collect();

        // Create placeholder loop handle (needed for half-edge construction).
        let inner_loop = self.loops.insert(Loop {
            half_edge: Default::default(),
            face,
        });

        // Create edges and half-edges.
        let mut edges = Vec::with_capacity(num_segments);
        let mut inner_hes = Vec::with_capacity(num_segments);
        let mut twin_hes = Vec::with_capacity(num_segments);

        for i in 0..num_segments {
            let j = (i + 1) % num_segments;

            // Build pcurve for this segment (degree-1 NURBS in UV space).
            let pcurve = NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![uv_positions[i], uv_positions[j]],
                weights: vec![1.0, 1.0],
            };

            // Build 3D curve segment for the edge.
            let pos_i = self.vertices[vertices[i]].position;
            let pos_j = self.vertices[vertices[j]].position;
            let edge_curve = NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![pos_i, pos_j],
                weights: vec![1.0, 1.0],
            };

            // Create the edge.
            let edge_h = self.insert_edge_tracked(Edge {
                half_edge: Default::default(),
                curve: Some(CurveOnEdge {
                    geometry: CurveGeometry::Nurbs(edge_curve),
                    t_start: 0.0,
                    t_end: 1.0,
                }),
                name: None,
                tolerance: None,
                canonical: Default::default(),
                role: Default::default(),
            });
            edges.push(edge_h);

            // Inner half-edge: goes from vertex[i] to vertex[j] on this face.
            let inner_he = self.half_edges.insert(HalfEdge {
                origin: vertices[i],
                twin: Default::default(), // patched below
                next: Default::default(), // patched below
                prev: Default::default(), // patched below
                edge: edge_h,
                loop_: inner_loop,
                pcurve: Some(pcurve.into()),
            });
            inner_hes.push(inner_he);

            // Twin half-edge: self-twinned, belongs to no loop (will be stitched later).
            // Origin is vertex[j] (opposite direction).
            let twin_he = self.half_edges.insert(HalfEdge {
                origin: vertices[j],
                twin: inner_he,
                next: Default::default(), // self-twinned twins have no chain
                prev: Default::default(),
                edge: edge_h,
                loop_: inner_loop, // needs some valid handle
                pcurve: None,
            });
            twin_hes.push(twin_he);

            // Patch inner half-edge twin.
            self.half_edges[inner_he].twin = twin_he;

            // Update edge to point to its inner half-edge.
            self.edges[edge_h].half_edge = inner_he;

            // Update vertex outgoing pointer.
            self.vertices[vertices[i]].outgoing = inner_he;
        }

        // Wire up next/prev chains for inner half-edges (closed loop).
        for i in 0..num_segments {
            let next_idx = (i + 1) % num_segments;
            let prev_idx = (i + num_segments - 1) % num_segments;
            self.half_edges[inner_hes[i]].next = inner_hes[next_idx];
            self.half_edges[inner_hes[i]].prev = inner_hes[prev_idx];
        }

        // Wire up next/prev chains for twin half-edges (reverse direction).
        // Twins go in the opposite winding: twin[i] origin = vertex[i+1],
        // so the chain is twin[n-1] -> twin[n-2] -> ... -> twin[0] -> twin[n-1].
        for i in 0..num_segments {
            let next_idx = (i + num_segments - 1) % num_segments;
            let prev_idx = (i + 1) % num_segments;
            self.half_edges[twin_hes[i]].next = twin_hes[next_idx];
            self.half_edges[twin_hes[i]].prev = twin_hes[prev_idx];
        }

        // Finalize the inner loop: point it at the first inner half-edge.
        self.loops[inner_loop].half_edge = inner_hes[0];

        // Add the inner loop to the face.
        self.faces[face].inner_loops.push(inner_loop);

        Ok(AddInnerLoopResult {
            loop_: inner_loop,
            vertices,
            edges,
        })
    }
}

#[cfg(test)]
#[path = "trim_tests.rs"]
mod tests;
