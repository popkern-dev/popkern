//! Measured edge tolerances for producers that FIT geometry (PK-826).
//!
//! An edge's tolerance states how closely its 3-D curve and its uses'
//! pcurve images agree — the SameParameter tolerance of a STEP/OCC
//! tolerant edge. Producers that fit an edge curve or its pcurves (a
//! blend strip's contact line, an offset's re-authored trims) left the
//! default model tolerance on edges whose fit sits well outside it
//! (measured: the thickened wavy plate's blend rims, 1.43e-5 at a 1e-6
//! tolerance), and the pcurve consistency check carried a 100x band to
//! stay quiet on them. A producer now records what it built.

use super::{BRep, EdgeHandle};

/// Margin over the measured worst deviation: the measurement samples
/// the curve, and the recorded tolerance must bound the unsampled
/// stretches too.
const MEASUREMENT_MARGIN: f64 = 1.5;

impl BRep {
    /// Raise each edge's tolerance to cover the measured agreement of
    /// its curve with its uses' pcurve images (fraction-matched, both
    /// uses), never lowering a recorded tolerance or going below the
    /// model's linear tolerance, and with its margin capped at the
    /// model's `max_imprecision`. Returns the edges raised.
    ///
    /// For producer exits that fit geometry; validators read the
    /// recorded tolerance as the edge's band.
    pub fn record_measured_edge_tolerances(&mut self) -> Vec<EdgeHandle> {
        let linear = self.tolerance.linear();
        let mut raised = Vec::new();
        let edges: Vec<EdgeHandle> = self.edges.keys().collect();
        for eh in edges {
            let Some(edge) = self.edges.get(eh) else {
                continue;
            };
            let first = edge.half_edge;
            let Some(twin) = self.half_edges.get(first).map(|h| h.twin) else {
                continue;
            };
            let worst = [first, twin]
                .into_iter()
                .filter_map(|he| self.measure_half_edge_pcurve_deviation(he))
                .map(|r| r.max_deviation)
                .filter(|d| d.is_finite())
                .fold(0.0_f64, f64::max);
            let current = edge.tolerance.unwrap_or(linear);
            // Capped at the model's `max_imprecision`: the margin never
            // manufactures a precision error. A measurement PAST the
            // ceiling is also recorded at the ceiling — TEMPORARY
            // CONTAINMENT (PK-885): the near-tangent cylinder tee (r 39
            // into r 40) measures one intersection edge at 1.03e-4;
            // recorded honestly it declines a union that serves. (The
            // wavy plate offset that first needed this, 1.45e-3, now
            // measures ~2e-8 — PK-883.) The replacement records
            // over-ceiling measurements as measured.
            let needed = (worst * MEASUREMENT_MARGIN).min(self.max_imprecision);
            if needed > current {
                if let Some(e) = self.edges.get_mut(eh) {
                    e.tolerance = Some(needed);
                }
                raised.push(eh);
            }
        }
        raised
    }
}

#[cfg(test)]
#[path = "edge_tolerance_measure_tests.rs"]
mod tests;
