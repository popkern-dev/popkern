//! Geometric deviation measurement for BRep edges and vertices.
//!
//! Provides methods to measure actual geometric deviations between
//! a BRep's topological entities and their underlying geometry.

use super::brep::BRep;
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;
use crate::topology::{EdgeHandle, VertexHandle};

/// Number of sample points along an edge curve for deviation measurement.
const DEFAULT_EDGE_SAMPLES: usize = 20;

impl BRep {
    /// Maximum deviation between an edge's 3D curve and the surface
    /// evaluations via pcurves on adjacent faces.
    ///
    /// For each adjacent half-edge that has both a pcurve and a face
    /// surface, samples the edge curve at `N` parameter values, evaluates
    /// the pcurve at each to obtain `(u, v)`, then evaluates the surface
    /// at `(u, v)` and compares with the 3D curve point.
    ///
    /// Returns `0.0` if the edge has no curve or no adjacent pcurves.
    pub fn measure_edge_deviation(&self, eh: EdgeHandle) -> f64 {
        let Some(edge) = self.edges.get(eh) else {
            return 0.0;
        };
        let Some(ref curve_on_edge) = edge.curve else {
            return 0.0;
        };

        let n = DEFAULT_EDGE_SAMPLES;
        let mut max_dev = 0.0_f64;

        // Examine both half-edges of this edge.
        let he_handle = edge.half_edge;
        let he_handles = if let Some(he) = self.half_edges.get(he_handle) {
            vec![he_handle, he.twin]
        } else {
            return 0.0;
        };

        for heh in he_handles {
            let Some(he) = self.half_edges.get(heh) else {
                continue;
            };
            let Some(ref pcurve) = he.pcurve else {
                continue;
            };

            // Get the face's surface via loop -> face.
            let Some(loop_) = self.loops.get(he.loop_) else {
                continue;
            };
            let Some(face) = self.faces.get(loop_.face) else {
                continue;
            };
            let Some(ref surface_on_face) = face.surface else {
                continue;
            };

            // Determine half-edge direction relative to edge curve.
            // The canonical half-edge (edge.half_edge) runs t_start -> t_end.
            // The twin runs t_end -> t_start.
            let forward = heh == he_handle;
            let (t_lo, t_hi) = if forward {
                (curve_on_edge.t_start, curve_on_edge.t_end)
            } else {
                (curve_on_edge.t_end, curve_on_edge.t_start)
            };

            // Sample along the edge curve, mapping pcurve fraction
            // to the correct edge curve direction.
            let pc_range = pcurve.parameter_range();
            for i in 0..=n {
                let frac = i as f64 / n as f64;
                let t = t_lo + (t_hi - t_lo) * frac;
                let p3d = curve_on_edge.geometry.point_at(t);

                // Evaluate pcurve at corresponding parameter.
                let pc_t = pc_range.0 + (pc_range.1 - pc_range.0) * frac;
                let uv = pcurve.point_at(pc_t);
                let p_surface = surface_on_face.geometry.point_at(uv.x, uv.y);

                let dev = p3d.distance_to(&p_surface);
                max_dev = max_dev.max(dev);
            }
        }

        max_dev
    }

    /// Maximum distance from a vertex position to incident edge curves
    /// evaluated at the vertex parameter.
    ///
    /// Walks the vertex fan (outgoing half-edges) and for each incident
    /// edge that has a curve, evaluates the curve at `t_start` or `t_end`
    /// depending on which end corresponds to this vertex.
    ///
    /// Returns `0.0` if no incident edges have curves.
    pub fn measure_vertex_deviation(&self, vh: VertexHandle) -> f64 {
        let Some(vertex) = self.vertices.get(vh) else {
            return 0.0;
        };
        let pos = vertex.position;

        // Walk the vertex fan: start at outgoing, follow twin.next.
        let start = vertex.outgoing;
        if !self.half_edges.contains_key(start) {
            return 0.0;
        }

        let mut max_dev = 0.0_f64;
        let mut cursor = start;
        let bound = self.half_edges.len();
        let mut steps = 0;

        loop {
            steps += 1;
            if steps > bound {
                break; // Safety: avoid infinite loop on corrupted data.
            }

            let Some(he) = self.half_edges.get(cursor) else {
                break;
            };
            let edge_handle = he.edge;

            if let Some(edge) = self.edges.get(edge_handle)
                && let Some(ref curve_on_edge) = edge.curve
            {
                // Determine which end of the curve corresponds to this vertex.
                // The edge's canonical half-edge origin is at t_start.
                // Its twin's origin is at t_end.
                let t = if let Some(canonical_he) = self.half_edges.get(edge.half_edge) {
                    if canonical_he.origin == vh {
                        curve_on_edge.t_start
                    } else {
                        curve_on_edge.t_end
                    }
                } else {
                    // Fallback: can't determine parameter, skip.
                    cursor = self.advance_vertex_fan(cursor, start);
                    if cursor == start {
                        break;
                    }
                    continue;
                };

                let curve_pt = curve_on_edge.geometry.point_at(t);
                let dev = pos.distance_to(&curve_pt);
                max_dev = max_dev.max(dev);
            }

            cursor = self.advance_vertex_fan(cursor, start);
            if cursor == start {
                break;
            }
        }

        max_dev
    }

    /// Advance the vertex fan cursor: follow twin.next.
    ///
    /// Returns `start` if the fan is complete or a stale handle is
    /// encountered, signalling the caller to stop.
    fn advance_vertex_fan(
        &self,
        cursor: crate::topology::HalfEdgeHandle,
        start: crate::topology::HalfEdgeHandle,
    ) -> crate::topology::HalfEdgeHandle {
        let Some(he) = self.half_edges.get(cursor) else {
            return start;
        };
        let Some(twin_data) = self.half_edges.get(he.twin) else {
            return start;
        };
        twin_data.next
    }
}

#[cfg(test)]
#[path = "brep_measurement_tests.rs"]
mod brep_measurement_tests;
