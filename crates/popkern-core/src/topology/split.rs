//! Face and edge splitting operations for surface trimming.

#[path = "split_edge_at_point.rs"]
mod split_edge_at_point;
#[path = "split_position.rs"]
mod split_position;
#[path = "split_surface_derived.rs"]
mod split_surface_derived;

use crate::error::CoreError;
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
        let pcurve_ab = he_ab_data.pcurve.clone();

        // Boundary edges have a null twin — handle them with a
        // simplified split that only touches the surviving side.
        let is_boundary = self.boundary_edges.contains(&edge);
        if is_boundary {
            return self.split_boundary_edge(
                edge,
                t,
                he_ab,
                vertex_a,
                he_ab_next,
                he_ab_loop,
                &pcurve_ab,
                &curve_on_edge,
            );
        }

        let he_ba_data = self
            .half_edges
            .get(he_ba)
            .ok_or_else(|| CoreError::InvalidHandle("twin half_edge is invalid".into()))?;

        let vertex_b = he_ba_data.origin;
        let he_ba_prev = he_ba_data.prev;
        let he_ba_loop = he_ba_data.loop_;
        let pcurve_ba = he_ba_data.pcurve.clone();

        // Compute the position of the new vertex.
        let pos_a = self.vertices[vertex_a].position;
        let pos_b = self.vertices[vertex_b].position;

        // The caller passes `t` in the half-edge's direction:
        // `t = 0` is at `vertex_a`, `t = 1` is at `vertex_b`. The
        // edge's curve geometry, however, may be authored in the
        // *opposite* direction (a common pattern in boolean-engine
        // outputs). So we detect the curve direction by checking
        // which curve endpoint is closer to `vertex_a`, and flip
        // `t` if needed before sampling the curve. Without this,
        // `split_pos` lands at the wrong end of the edge.
        let (split_pos, edge_a_curve_replacement, edge_b_curve) =
            split_position::split_position_and_propagation(
                self,
                &curve_on_edge,
                he_ab,
                he_ba,
                pos_a,
                pos_b,
                t,
            );

        // Create the new vertex V.
        let vertex_v = self.vertices.insert(Vertex {
            position: split_pos,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });

        // Split pcurves for he_ab (A->B direction) and he_ba (B->A direction).
        // For he_ab: split at parameter t → left goes to he_ab [A->V], right to he_vb [V->B].
        // For he_ba: split at (1-t) → left goes to he_bv [B->V], right to he_ba [V->A].
        let surface_of = |he: HalfEdgeHandle| -> Option<crate::geometry::surface::SurfaceGeometry> {
            let loop_ = self.half_edges.get(he)?.loop_;
            let face = self.loops.get(loop_)?.face;
            self.faces
                .get(face)?
                .surface
                .as_ref()
                .map(|s| s.geometry.clone())
        };
        let (surface_ab, surface_ba) = (surface_of(he_ab), surface_of(he_ba));
        let linear = self.tolerance.linear();
        let (pcurve_ab_left, pcurve_ab_right) = super::split_helpers::split_pcurve_at_point(
            &pcurve_ab,
            t,
            surface_ab.as_ref(),
            &split_pos,
            linear,
        );
        let (pcurve_ba_left, pcurve_ba_right) = super::split_helpers::split_pcurve_at_point(
            &pcurve_ba,
            1.0 - t,
            surface_ba.as_ref(),
            &split_pos,
            linear,
        );

        // Create the new edge (edge_b), inheriting the V→B portion
        // of the parent's curve when present (β-1-A).
        let edge_b = self.insert_edge_tracked(Edge {
            half_edge: Default::default(),
            curve: edge_b_curve,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });

        // Create new half-edges: he_vb (V->B) and he_bv (B->V).
        let he_vb = self.half_edges.insert(HalfEdge {
            origin: vertex_v,
            twin: Default::default(),
            next: he_ab_next,
            prev: he_ab,
            edge: edge_b,
            loop_: he_ab_loop,
            pcurve: pcurve_ab_right,
        });

        let he_bv = self.half_edges.insert(HalfEdge {
            origin: vertex_b,
            twin: he_vb,
            next: he_ba,
            prev: he_ba_prev,
            edge: edge_b,
            loop_: he_ba_loop,
            pcurve: pcurve_ba_left,
        });

        // Patch twin of he_vb.
        self.half_edges[he_vb].twin = he_bv;

        if he_ab_next == he_ba {
            // SEAM TURNAROUND at B (he_ab and its twin are adjacent in
            // one loop — a cone apex seam, or a bigon seam loop). The
            // generic four rewires alias here (he_ab_next == he_ba and
            // he_ba_prev == he_ab) and orphan he_vb from the cycle
            // (measured: the cone-apex seam split left mef's loop walk
            // spinning). Desired chain:
            // he_ab (A→V) → he_vb (V→B) → he_bv (B→V) → he_ba (V→A).
            self.half_edges[he_ab].next = he_vb;
            self.half_edges[he_vb].prev = he_ab;
            self.half_edges[he_vb].next = he_bv;
            self.half_edges[he_bv].prev = he_vb;
            self.half_edges[he_bv].next = he_ba;
            self.half_edges[he_ba].prev = he_bv;
        } else {
            // Rewire the chain on he_ab's side:
            // he_ab -> he_vb -> he_ab_next  (was: he_ab -> he_ab_next)
            self.half_edges[he_ab].next = he_vb;
            self.half_edges[he_ab_next].prev = he_vb;

            // Rewire the chain on he_ba's side:
            // he_ba_prev -> he_bv -> he_ba  (was: he_ba_prev -> he_ba)
            self.half_edges[he_ba_prev].next = he_bv;
            self.half_edges[he_ba].prev = he_bv;
        }

        // Update he_ba origin: was B, now V (since the original edge A->B
        // becomes A->V, its twin goes V->A).
        self.half_edges[he_ba].origin = vertex_v;

        // Assign split pcurves to the original half-edges.
        self.half_edges[he_ab].pcurve = pcurve_ab_left;
        self.half_edges[he_ba].pcurve = pcurve_ba_right;

        // Update edge_b to point to he_vb.
        self.edges[edge_b].half_edge = he_vb;

        // β-1-A: shorten the original edge's curve to its A→V
        // portion (a no-op when the parent had no curve).
        self.edges[edge].curve = edge_a_curve_replacement;

        // Update vertex outgoing pointers.
        self.vertices[vertex_v].outgoing = he_vb;

        // If B's outgoing was he_ba (which now originates at V), update it to he_bv.
        if self.vertices[vertex_b].outgoing == he_ba {
            self.vertices[vertex_b].outgoing = he_bv;
        }

        // Split curve geometry if present (routing lives in
        // `split_position::route_split_curves` — direction decision and
        // half-to-edge assignment must agree with `split_pos` above).
        split_position::route_split_curves(
            self,
            &curve_on_edge,
            he_ab,
            he_ba,
            &pcurve_ab,
            &pcurve_ba,
            edge,
            edge_b,
            pos_a,
            pos_b,
            split_pos,
            t,
        );

        Ok(SplitEdgeResult {
            vertex: vertex_v,
            edge_a: edge,
            edge_b,
        })
    }

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
        trim_curve_uv: &NurbsCurve,
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

        // NEAR-VERTEX rescue for endpoints the strict pass missed —
        // see `nearest_loop_vertex_within_band` for the trust model.
        if he_entry.is_none() {
            he_entry = self.nearest_loop_vertex_within_band(
                loop_start_he,
                &curve_start,
                he_exit,
                tolerance,
            );
        }
        if he_exit.is_none() {
            he_exit = self.nearest_loop_vertex_within_band(
                loop_start_he,
                &curve_end,
                he_entry,
                tolerance,
            );
        }

        // INNER-LOOP cut (hole-rim band split, PK-185/PK-186): both
        // endpoints can live on one HOLE loop — interior-rail spur
        // tips at reflex hole corners; see `match_inner_loop_cut`.
        // Only fires when the outer loop failed to supply both
        // endpoints — previously a guaranteed precondition error.
        if (he_entry.is_none() || he_exit.is_none())
            && let Some((entry, exit)) =
                self.match_inner_loop_cut(face, &curve_start, &curve_end, tolerance)
        {
            he_entry = Some(entry);
            he_exit = Some(exit);
        }

        // Structured decline evidence: how far the endpoint actually
        // sits from the loop (a reader can classify the miss — fit
        // slop, amplified graze slide, or a genuinely absent vertex —
        // without re-running under a debugger).
        let nearest_miss = |target: &crate::geometry::Point3| -> f64 {
            let mut best = f64::INFINITY;
            let mut cursor = loop_start_he;
            for _ in 0..self.half_edges.len() + 1 {
                let Some(he_data) = self.half_edges.get(cursor) else {
                    break;
                };
                if let Some(v) = self.vertices.get(he_data.origin) {
                    best = best.min(v.position.distance_to(target));
                }
                cursor = he_data.next;
                if cursor == loop_start_he {
                    break;
                }
            }
            best
        };
        let he_entry = he_entry.ok_or_else(|| {
            CoreError::PreconditionFailed(format!(
                "trim curve start point ({:.6},{:.6},{:.6}) does not match any vertex on the \
                 face boundary (nearest loop vertex {:.3e} away)",
                curve_start.x,
                curve_start.y,
                curve_start.z,
                nearest_miss(&curve_start)
            ))
        })?;
        let he_exit = he_exit.ok_or_else(|| {
            CoreError::PreconditionFailed(format!(
                "trim curve end point ({:.6},{:.6},{:.6}) does not match any vertex on the \
                 face boundary (nearest loop vertex {:.3e} away)",
                curve_end.x,
                curve_end.y,
                curve_end.z,
                nearest_miss(&curve_end)
            ))
        })?;

        // Capture the original surface before mef (which doesn't copy it).
        let original_surface = self.faces[face].surface.clone();

        let entry_at_start = tolerance.points_coincident(
            &self.vertices[self.half_edges[he_entry].origin].position,
            &curve_start,
        );
        if crate::trace::trace_enabled("face-split") {
            eprintln!("    [split-face] entry_at_start={entry_at_start}");
        }

        // Call mef to split the face topology.
        let mef_result = self.mef(he_entry, he_exit)?;

        // Attach the trim curve geometry to the new edge, preserving analytic
        // line semantics when the trim is a straight span.
        let trim_geometry = CurveGeometry::Nurbs(trim_curve_3d.clone());
        let trim_geometry = trim_geometry
            .canonical_line_in_range((t_start, t_end), tolerance.linear().max(1.0e-9))
            .unwrap_or(trim_geometry);
        self.edges[mef_result.edge].curve = Some(CurveOnEdge {
            geometry: trim_geometry,
            t_start,
            t_end,
        });

        // Assign pcurves to the new half-edges by each half-edge's OWN
        // traversal direction. mef's half_edge runs origin(he_entry) →
        // origin(he_exit) — usually the trim's forward direction, but
        // the inner-loop matcher may return the pair SWAPPED relative
        // to the trim's parameterization (the winding rule owns
        // entry/exit order — PK-185/PK-186). A forward pcurve on a
        // backwards half-edge makes the loop's traversal disagree with
        // its geometry (measured, PK-428: both keyed-bore cap trims
        // took the swapped branch, every wall×cap pairing in the merge
        // graph then claimed the SAME direction and minted duplicate
        // edges — 8 self-twinned faces, non-manifold serve).
        let (forward_he, backward_he) = if entry_at_start {
            (mef_result.half_edge, mef_result.half_edge_twin)
        } else {
            (mef_result.half_edge_twin, mef_result.half_edge)
        };
        self.half_edges[forward_he].pcurve = Some(trim_curve_uv.clone().into());
        self.half_edges[backward_he].pcurve = Some(reverse_nurbs_curve(trim_curve_uv).into());

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

pub(crate) use super::split_helpers::{reverse_nurbs_curve, split_pcurve_at};

#[cfg(test)]
#[path = "split_tests.rs"]
mod tests;
