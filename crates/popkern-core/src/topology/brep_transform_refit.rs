//! The imprecision a scale-up carries past its edges' bands, recorded as
//! their tolerances.
//!
//! A uniform scale maps every chart value unchanged while edges and
//! surfaces grow, so whatever the input's construction left between an
//! edge's curve and its faces grows with them: a NURBS torus whose edges
//! sat 9e-7 off its surfaces at unit size sat 9e-4 off at a thousand,
//! past the 1e-5 band of the unscaled tolerance, and every boolean on
//! the scaled body failed validation (generated corpus, scale 1000).
//! Re-fitting the pcurves cannot help (the edge curve itself is off the
//! surface by that much). The scaled body states it: an edge whose
//! pcurve image strays past its band takes the measured deviation as
//! its tolerance, as an import records a fit's accuracy (pcurve
//! consistency validation: fit accuracy belongs in an honest edge
//! tolerance), and the model's imprecision budget grows with the body.
//! An edge's curve ending off its vertex grows the same way: a box split
//! by a NURBS dome kept curve ends 1.2e-9 off their vertices at unit
//! size, 1.2e-6 off at a thousand, past the vertices' 1e-6, and the
//! scaled operand failed closure before any boolean ran (generated
//! corpus W at scale 1000: SW-E 772, WS-I 898, YW-E 1380). The vertex
//! takes the miss as its tolerance, as an import records a file's
//! precision there (PK-639).

use crate::geometry::curve::CurveEval;
use crate::topology::{BRep, EdgeHandle, VertexHandle};

/// The band in units of an edge's tolerance, as pcurve consistency
/// validation reads it.
const BAND_FACTOR: f64 = 10.0;

impl BRep {
    /// Give every edge whose pcurve image strays past its band (its
    /// tolerance, else its face's, else the model's, times ten) the
    /// worst measured deviation as its tolerance, and scale the model's
    /// imprecision budget with the body (the imprecision it holds grew
    /// by `scale`; the budget fixed, a valid body scaled to an invalid
    /// one).
    pub(super) fn record_scaled_imprecision(&mut self, scale: f64) {
        self.max_imprecision *= scale;
        let mut widen: Vec<(EdgeHandle, f64)> = Vec::new();
        for (half_edge, data) in &self.half_edges {
            if data.pcurve.is_none() {
                continue;
            }
            let Some(edge) = self.edges.get(data.edge) else {
                continue;
            };
            let Some(report) = self.measure_half_edge_pcurve_deviation(half_edge) else {
                continue;
            };
            let face_tolerance = self.faces.get(report.face).and_then(|f| f.tolerance);
            let own = edge
                .tolerance
                .or(face_tolerance)
                .unwrap_or_else(|| self.tolerance.linear())
                .max(1.0e-9);
            if report.max_deviation.is_finite() && report.max_deviation > own * BAND_FACTOR {
                widen.push((data.edge, report.max_deviation));
            }
        }
        for (edge, deviation) in widen {
            if let Some(e) = self.edges.get_mut(edge) {
                e.tolerance = Some(e.tolerance.unwrap_or(0.0).max(deviation));
            }
        }
        self.record_scaled_end_misses();
    }

    /// Give every vertex an edge's curve ends farther from than its
    /// tolerance that distance as its tolerance (each end read against
    /// the nearer of the edge's two vertices' pairing, as closure
    /// validation reads a curve stored either way round), and its edges
    /// at least as much.
    fn record_scaled_end_misses(&mut self) {
        let mut widen: Vec<(VertexHandle, f64)> = Vec::new();
        for (edge, data) in &self.edges {
            let (Some(c), Some((a, b))) = (data.curve.as_ref(), self.edge_vertices(edge)) else {
                continue;
            };
            let (p0, p1) = (c.geometry.point_at(c.t_start), c.geometry.point_at(c.t_end));
            let (pa, pb) = (self.vertices[a].position, self.vertices[b].position);
            let natural = (pa.distance_to(&p0), pb.distance_to(&p1));
            let swapped = (pa.distance_to(&p1), pb.distance_to(&p0));
            let (da, db) = if natural.0.max(natural.1) <= swapped.0.max(swapped.1) {
                natural
            } else {
                swapped
            };
            for (v, miss) in [(a, da), (b, db)] {
                if miss.is_finite() && miss > self.vertex_tolerance(v) {
                    widen.push((v, miss));
                }
            }
        }
        if widen.is_empty() {
            return;
        }
        for &(v, miss) in &widen {
            if let Some(vertex) = self.vertices.get_mut(v) {
                vertex.tolerance = Some(vertex.tolerance.unwrap_or(0.0).max(miss));
            }
        }
        // A vertex's tolerance stays within each of its edges' (the
        // precision hierarchy): an edge at a widened vertex widens to it.
        let edges: Vec<(EdgeHandle, f64)> = self
            .edges
            .keys()
            .filter_map(|e| {
                let (a, b) = self.edge_vertices(e)?;
                if !widen.iter().any(|&(v, _)| v == a || v == b) {
                    return None;
                }
                let need = self.vertex_tolerance(a).max(self.vertex_tolerance(b));
                (need > self.edge_tolerance(e)).then_some((e, need))
            })
            .collect();
        for (e, need) in edges {
            if let Some(edge) = self.edges.get_mut(e) {
                edge.tolerance = Some(need);
            }
        }
    }
}

#[cfg(test)]
#[path = "brep_transform_refit_tests.rs"]
mod tests;
