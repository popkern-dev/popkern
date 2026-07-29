//! Loop boundary sampling: vertex positions plus interior samples
//! along edge curves, in half-edge traversal order.
//!
//! Split from `brep.rs` for the structural line ratchet.

use crate::topology::{BRep, FaceHandle, LoopHandle};

impl BRep {
    /// Collect boundary sample positions from a face's outer loop:
    /// every vertex position PLUS interior samples along each edge
    /// CURVE. Vertex-only sampling is blind to faces bounded by
    /// closed or long curved edges — a cylinder cap's rim carries
    /// exactly ONE vertex, so a parameter box derived from vertices
    /// collapses to a point at the seam and intersection discovery
    /// misses the entire face interior (measured: any pocket cut
    /// into a swept cylinder's cap silently no-opped). Curve
    /// samples span the true face extent.
    ///
    /// Edges without curve geometry contribute straight-line
    /// interpolants between their endpoints.
    pub fn face_outer_loop_sample_positions(
        &self,
        face_handle: FaceHandle,
        samples_per_edge: usize,
    ) -> Vec<crate::geometry::Point3> {
        let Some(face) = self.faces.get(face_handle) else {
            return Vec::new();
        };
        self.loop_sample_positions(face.outer_loop, samples_per_edge)
    }

    /// Boundary sample positions of ANY loop (outer or inner), in
    /// half-edge traversal order: every vertex position plus interior
    /// samples along each edge curve (straight interpolants for
    /// curveless edges). The inner-loop counterpart of
    /// [`Self::face_outer_loop_sample_positions`] — hole boundaries
    /// carry traversal semantics of their own (they wind opposite the
    /// outer loop) and consumers matching flow against them need the
    /// samples in stored order.
    pub fn loop_sample_positions(
        &self,
        loop_handle: LoopHandle,
        samples_per_edge: usize,
    ) -> Vec<crate::geometry::Point3> {
        use crate::geometry::curve::CurveEval;
        let Some(loop_) = self.loops.get(loop_handle) else {
            return Vec::new();
        };
        let start_he = loop_.half_edge;
        if !self.half_edges.contains_key(start_he) {
            return Vec::new();
        }
        let mut positions = Vec::new();
        let mut current = start_he;
        let limit = self.half_edges.len() + 1;
        for _ in 0..limit {
            let he = &self.half_edges[current];
            let origin = self.vertices[he.origin].position;
            positions.push(origin);
            if samples_per_edge > 0 {
                // Staleness guard: only trust a curve whose endpoints
                // agree with the edge's vertices — fixtures (and some
                // legacy passes) move vertices without re-authoring
                // curves, and sampling a stale curve inflates the
                // position cloud with points from the old geometry
                // (measured: a vertex-scaled inner cube reported its
                // pre-scale extent and containment booleans collapsed).
                let dest = self.vertices[self.half_edges[he.twin].origin].position;
                let fresh_curve = self
                    .edges
                    .get(he.edge)
                    .and_then(|e| e.curve.as_ref())
                    .filter(|curve| {
                        let ceval: &dyn CurveEval = &curve.geometry;
                        let a = ceval.point_at(curve.t_start);
                        let b = ceval.point_at(curve.t_end);
                        (a.distance_to(&origin) < 1e-6 && b.distance_to(&dest) < 1e-6)
                            || (a.distance_to(&dest) < 1e-6 && b.distance_to(&origin) < 1e-6)
                    });
                if let Some(curve) = fresh_curve {
                    let ceval: &dyn CurveEval = &curve.geometry;
                    // Sample in TRAVERSAL order: the stored curve may
                    // run opposite the half-edge (the staleness guard
                    // above accepts either match), and emitting its
                    // samples in raw parameter order zigzags the
                    // polygon — consumers matching flow direction
                    // against the loop (boundary-substitute
                    // resolution) read alternating segment directions
                    // and their votes cancel.
                    let reversed = ceval.point_at(curve.t_start).distance_to(&origin) >= 1e-6;
                    let (s0, s1) = if reversed {
                        (curve.t_end, curve.t_start)
                    } else {
                        (curve.t_start, curve.t_end)
                    };
                    for i in 1..=samples_per_edge {
                        let t = s0 + (s1 - s0) * i as f64 / (samples_per_edge + 1) as f64;
                        positions.push(ceval.point_at(t));
                    }
                } else {
                    for i in 1..=samples_per_edge {
                        let f = i as f64 / (samples_per_edge + 1) as f64;
                        positions.push(origin + (dest - origin) * f);
                    }
                }
            }
            current = he.next;
            if current == start_he {
                break;
            }
        }
        positions
    }
}
