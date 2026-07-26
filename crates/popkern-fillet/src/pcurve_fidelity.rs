//! Same-fraction pcurve fidelity — a *validated* instrument.
//!
//! # Why this module exists at all
//!
//! Three separate ad-hoc probes of this quantity were written during
//! the PnS campaign and two of them were wrong, one badly enough to
//! produce a retracted root cause. The failures were not subtle once
//! seen, but they were invisible without a known-good reference:
//!
//! - **Orientation.** A half-edge pcurve runs along its *loop*, which
//!   for one of the two half-edges opposes the edge curve's own
//!   direction. Sampling both forward reads a perfectly good reversed
//!   pcurve as a deviation of `~0.9 * chord`. Every "defect" found
//!   that way showed a dev/chord ratio of exactly 0.90.
//! - **Scope.** Measured across every face of a mid-pipeline BRep,
//!   the answer was 4163 on a unit-scale model, because transient
//!   faces carry incomplete pcurves that this question is meaningless
//!   for.
//!
//! So the instrument is here, in library code, with a validation test
//! that pins both a known-good case and a known-bad one. An
//! instrument that has only ever been run on the thing it is
//! investigating is not evidence.
//!
//! # What it measures
//!
//! For each half-edge of the requested faces: the distance between
//! the pcurve's image on its face surface and the edge's 3D curve, at
//! matched fractions, with orientation resolved from the endpoints.
//! That is the agreement validation and tessellation consume — the
//! same quantity `cap_outer_acceptance_measure` gates on.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Samples per half-edge, inset from both ends.
const SAMPLES: usize = 21;
/// Fractional inset. Endpoints are excluded because a pcurve is
/// pinned there by construction and would flatter the measurement.
const END_MARGIN: f64 = 0.05;

/// Worst same-fraction deviation found, and where.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PcurveFidelity {
    /// Worst deviation in model units.
    pub worst: f64,
    /// The edge carrying it.
    pub worst_edge: EdgeHandle,
    /// The face whose pcurve deviated.
    pub worst_face: FaceHandle,
    /// Fraction along the edge where it occurred.
    pub worst_fraction: f64,
    /// Half-edges actually measured.
    pub measured: usize,
}

/// Measure pcurve fidelity over `faces`.
///
/// Returns `None` when nothing was measurable — never a zero, which
/// would read as a perfect score and is precisely how an unmeasured
/// thing comes to look like a passing one.
pub fn measure_pcurve_fidelity(brep: &BRep, faces: &[FaceHandle]) -> Option<PcurveFidelity> {
    let mut best: Option<PcurveFidelity> = None;
    let mut measured = 0usize;
    for &face in faces {
        let Some(face_data) = brep.faces.get(face) else {
            continue;
        };
        let Some(surface) = face_data.surface.as_ref() else {
            continue;
        };
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut he = start;
        for _ in 0..64 {
            let Some(hd) = brep.half_edges.get(he) else {
                break;
            };
            if let Some(pc) = hd.pcurve.as_ref()
                && let Some(curve) = brep.edges.get(hd.edge).and_then(|e| e.curve.as_ref())
            {
                let img_at = |t: f64| {
                    let uv = pc.point_at(t);
                    surface.geometry.point_at(uv.x, uv.y)
                };
                let c_at = |t: f64| (&curve.geometry as &dyn CurveEval).point_at(t);
                let (c_lo, c_hi) = (c_at(curve.t_start), c_at(curve.t_end));
                // Resolve orientation from the endpoints before
                // measuring anything — see the module docs.
                let forward = c_lo.distance_to(&img_at(pc.t_start))
                    + c_hi.distance_to(&img_at(pc.t_end))
                    <= c_lo.distance_to(&img_at(pc.t_end)) + c_hi.distance_to(&img_at(pc.t_start));
                measured += 1;
                for k in 0..SAMPLES {
                    let f = END_MARGIN + (1.0 - 2.0 * END_MARGIN) * k as f64 / (SAMPLES - 1) as f64;
                    let s = if forward { f } else { 1.0 - f };
                    let img = img_at(pc.t_start + s * (pc.t_end - pc.t_start));
                    let cp = c_at(curve.t_start + f * (curve.t_end - curve.t_start));
                    let d = img.distance_to(&cp);
                    if best.is_none_or(|b| d > b.worst) {
                        best = Some(PcurveFidelity {
                            worst: d,
                            worst_edge: hd.edge,
                            worst_face: face,
                            worst_fraction: f,
                            measured,
                        });
                    }
                }
            }
            he = hd.next;
            if he == start {
                break;
            }
        }
    }
    best.map(|b| PcurveFidelity { measured, ..b })
}

#[cfg(test)]
#[path = "pcurve_fidelity_tests.rs"]
mod tests;
