//! Pcurve consistency: every half-edge's pcurve, mapped through its
//! face's surface, runs along its edge — same ends, same sense, same
//! trace.
//!
//! Check 9 measures the edge curve's distance to each adjacent surface
//! only, so a pcurve that disagrees with its edge passed validation
//! while every trim consumer (containment polygons, loop windings,
//! tessellation boundaries) read the pcurve: a moved countersink's
//! sub-arc kept its parent arc's whole pcurve (4.9 off, PK-820), and a
//! bored tube's rim pcurve ran against its edge (58 off, the hole on
//! its cap wound like the outer loop, PK-819).

use crate::error::ValidationError;
use crate::topology::BRep;

/// The band in units of the edge's tolerance. This check targets
/// pcurves on the wrong curve, sense or span (measured 0.03..58);
/// fitted pcurves on free-form and toroidal charts sit 1.4e-5..9e-5
/// off their edges at the default 1e-6 tolerance — fit accuracy that
/// belongs in an honest edge tolerance, not in this check.
const PCURVE_BAND_FACTOR: f64 = 100.0;

impl BRep {
    /// [`ValidationError::PcurveDeviation`] for every half-edge whose
    /// pcurve strays from its edge beyond the edge's band (its own
    /// tolerance, else its face's, else the model's linear tolerance,
    /// times [`PCURVE_BAND_FACTOR`]).
    pub(crate) fn validate_pcurve_consistency_all(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();
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
            let band = edge
                .tolerance
                .or(face_tolerance)
                .unwrap_or_else(|| self.tolerance.linear())
                .max(1.0e-9)
                * PCURVE_BAND_FACTOR;
            // A closed edge's sense belongs to the canonical-forward
            // audit (PK-478: twin claims must oppose); here its pcurve
            // must trace the edge and both its ends sit on the edge's
            // one vertex.
            let closed = self.edge_is_closed(data.edge);
            let deviation = if closed {
                let end_gap = self.pcurve_end_gap_to_origin(half_edge).unwrap_or(0.0);
                report.max_trace_deviation.max(end_gap)
            } else {
                report.max_deviation
            };
            if deviation.is_finite() && deviation > band {
                errors.push(ValidationError::PcurveDeviation {
                    half_edge,
                    face: report.face,
                    max_deviation: deviation,
                });
            }
        }
        errors
    }
}

impl BRep {
    /// Whether `edge` starts and ends at one vertex.
    fn edge_is_closed(&self, edge: crate::topology::EdgeHandle) -> bool {
        let Some(he) = self
            .edges
            .get(edge)
            .and_then(|e| self.half_edges.get(e.half_edge))
        else {
            return false;
        };
        self.half_edges
            .get(he.twin)
            .is_some_and(|tw| tw.origin == he.origin)
    }

    /// Worst distance from the surface images of `half_edge`'s pcurve
    /// ends to its origin vertex (a closed edge's only vertex).
    fn pcurve_end_gap_to_origin(&self, half_edge: crate::topology::HalfEdgeHandle) -> Option<f64> {
        use crate::geometry::surface::SurfaceEval;
        let he = self.half_edges.get(half_edge)?;
        let pcurve = he.pcurve.as_ref()?;
        let face = self.loops.get(he.loop_)?.face;
        let surface = &self.faces.get(face)?.surface.as_ref()?.geometry;
        let origin = self.vertices.get(he.origin)?.position;
        let image = |t: f64| {
            let q = pcurve.point_at(t);
            surface.point_at(q.x, q.y).distance_to(&origin)
        };
        Some(image(pcurve.t_start).max(image(pcurve.t_end)))
    }
}

#[cfg(test)]
#[path = "pcurve_consistency_tests.rs"]
mod tests;
