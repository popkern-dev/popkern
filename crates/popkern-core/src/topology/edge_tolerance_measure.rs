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
    /// model's `max_imprecision` (a measurement past the ceiling is
    /// recorded as measured). Returns the edges raised.
    ///
    /// For producer exits that fit geometry; validators read the
    /// recorded tolerance as the edge's band.
    pub fn record_measured_edge_tolerances(&mut self) -> Vec<EdgeHandle> {
        self.record_measured_edge_tolerances_except(&std::collections::BTreeSet::new())
    }

    /// [`Self::record_measured_edge_tolerances`] over every edge but
    /// `unchanged`: edges a producer carried over without re-fitting
    /// (their curve copied, their pcurves copied or re-derived by the
    /// measured projection fit, which holds a tenth of the linear
    /// tolerance). Measuring them is the whole-body cost a local
    /// operation must not pay — a boolean's 64th slot cut on a disc
    /// spent 249 ms of its 265 here (PK-826 regression).
    pub fn record_measured_edge_tolerances_except(
        &mut self,
        unchanged: &std::collections::BTreeSet<EdgeHandle>,
    ) -> Vec<EdgeHandle> {
        self.record_measured_edge_tolerances_where(|e| !unchanged.contains(&e))
    }

    /// [`Self::record_measured_edge_tolerances`] over the edges `measure`
    /// admits (a dense set's test where the unchanged edges are most of
    /// a body: a `BTreeSet` lookup per edge was a third of this pass on
    /// a 771-edge body, PK-815).
    pub fn record_measured_edge_tolerances_where(
        &mut self,
        measure: impl Fn(EdgeHandle) -> bool,
    ) -> Vec<EdgeHandle> {
        let linear = self.tolerance.linear();
        let mut raised = Vec::new();
        let edges: Vec<EdgeHandle> = self.edges.keys().filter(|&e| measure(e)).collect();
        for eh in edges {
            let Some(edge) = self.edges.get(eh) else {
                continue;
            };
            let first = edge.half_edge;
            let Some(twin) = self.half_edges.get(first).map(|h| h.twin) else {
                continue;
            };
            let current = edge.tolerance.unwrap_or(linear);
            // Screen: the trace-free bound on identical samples; an edge
            // whose bound cannot raise it skips the trace search.
            let bound = [first, twin]
                .into_iter()
                .filter_map(|he| self.half_edge_pcurve_deviation_bound(he))
                .map(|d| if d.is_finite() { d } else { f64::INFINITY })
                .fold(0.0_f64, f64::max);
            if bound * MEASUREMENT_MARGIN <= current {
                continue;
            }
            let worst = [first, twin]
                .into_iter()
                .filter_map(|he| self.measure_half_edge_pcurve_deviation(he))
                .map(|r| r.max_deviation)
                .filter(|d| d.is_finite())
                .fold(0.0_f64, f64::max);
            // The margin never manufactures a precision error: it is
            // capped at the model's `max_imprecision`. A measurement past
            // the ceiling is recorded as measured, and precision
            // validation reports it (the wavy plate offset's 1.45e-3
            // rails and a near-tangent tee's 1.03e-4 pcurves were
            // recorded at the ceiling until PK-883/PK-885 authored them
            // to the linear tolerance).
            let needed = if worst > self.max_imprecision {
                worst
            } else {
                (worst * MEASUREMENT_MARGIN).min(self.max_imprecision)
            };
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
