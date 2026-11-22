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

use crate::topology::{BRep, EdgeHandle};

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
    }
}
