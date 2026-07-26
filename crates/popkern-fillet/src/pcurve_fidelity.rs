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
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

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
    /// Worst deviation measured as NEAREST POINT on the edge curve
    /// rather than at matched fractions.
    ///
    /// The difference between this and [`Self::worst`] is diagnostic:
    /// a pcurve that traces the right curve at a different *rate* has
    /// a small nearest-point error and a large same-fraction one. That
    /// is a parameterisation mismatch, not a wrong pcurve, and the two
    /// call for completely different repairs — so they are reported
    /// separately rather than collapsed into one number.
    pub worst_nearest_point: f64,
}

/// Rim/spoke partition of a cap's boundary fidelity.
///
/// Each half-edge on a cap face's outer loop is classified by the face
/// on the far side of its edge: another cap face makes it a **spoke**
/// (an internal cap–cap seam), anything else makes it **rim** (the
/// boundary the cap shares with its supports). The two families answer
/// different questions — rim error is the cap's border construction
/// against the geometry it was asked to interpolate; spoke error is
/// the cap's internal seam wiring — so they are reported separately
/// rather than collapsed into one number that cannot attribute.
#[derive(Debug, Clone, Copy)]
pub struct CapBoundaryFidelitySplit {
    /// Fidelity over rim half-edges (cap-vs-support seams).
    /// `None` when no rim half-edge was measurable.
    pub rim: Option<PcurveFidelity>,
    /// Fidelity over spoke half-edges (cap-vs-cap seams).
    /// `None` when no spoke half-edge was measurable.
    pub spoke: Option<PcurveFidelity>,
}

/// Measure the rim/spoke split over `cap_faces`.
///
/// The partition is exhaustive over the measurable half-edges of the
/// cap faces' outer loops: every half-edge `measure_pcurve_fidelity`
/// would count lands in exactly one of the two families.
pub fn measure_cap_boundary_fidelity_split(
    brep: &BRep,
    cap_faces: &[FaceHandle],
) -> CapBoundaryFidelitySplit {
    let cap_set: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();
    let neighbour_is_cap = |he: HalfEdgeHandle| -> bool {
        brep.half_edges
            .get(he)
            .and_then(|hd| brep.half_edges.get(hd.twin))
            .and_then(|twin| brep.loops.get(twin.loop_))
            .map(|l| cap_set.contains(&l.face))
            .unwrap_or(false)
    };
    CapBoundaryFidelitySplit {
        rim: measure_pcurve_fidelity_filtered(brep, cap_faces, &|he| !neighbour_is_cap(he)),
        spoke: measure_pcurve_fidelity_filtered(brep, cap_faces, &neighbour_is_cap),
    }
}

/// Measure pcurve fidelity over `faces`.
///
/// Returns `None` when nothing was measurable — never a zero, which
/// would read as a perfect score and is precisely how an unmeasured
/// thing comes to look like a passing one.
pub fn measure_pcurve_fidelity(brep: &BRep, faces: &[FaceHandle]) -> Option<PcurveFidelity> {
    measure_pcurve_fidelity_filtered(brep, faces, &|_| true)
}

/// Shared measuring core: identical to [`measure_pcurve_fidelity`]
/// but restricted to half-edges accepted by `keep`. The filter exists
/// so partitioned measurements (e.g. rim/spoke) reuse this exact loop
/// instead of re-implementing it — three ad-hoc copies of this
/// measurement were written during the PnS campaign and two were
/// wrong.
fn measure_pcurve_fidelity_filtered(
    brep: &BRep,
    faces: &[FaceHandle],
    keep: &dyn Fn(HalfEdgeHandle) -> bool,
) -> Option<PcurveFidelity> {
    let mut best: Option<PcurveFidelity> = None;
    let mut measured = 0usize;
    let mut worst_nearest = 0.0_f64;
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
            if keep(he)
                && let Some(pc) = hd.pcurve.as_ref()
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
                // Dense reference samples for the nearest-point arm.
                let reference: Vec<_> = (0..=64)
                    .map(|r| {
                        c_at(curve.t_start + f64::from(r) / 64.0 * (curve.t_end - curve.t_start))
                    })
                    .collect();
                for k in 0..SAMPLES {
                    let f = END_MARGIN + (1.0 - 2.0 * END_MARGIN) * k as f64 / (SAMPLES - 1) as f64;
                    let s = if forward { f } else { 1.0 - f };
                    let img = img_at(pc.t_start + s * (pc.t_end - pc.t_start));
                    let cp = c_at(curve.t_start + f * (curve.t_end - curve.t_start));
                    let d = img.distance_to(&cp);
                    let nearest = reference
                        .iter()
                        .map(|r| img.distance_to(r))
                        .fold(f64::MAX, f64::min);
                    worst_nearest = worst_nearest.max(nearest);
                    if best.is_none_or(|b| d > b.worst) {
                        best = Some(PcurveFidelity {
                            worst: d,
                            worst_edge: hd.edge,
                            worst_face: face,
                            worst_fraction: f,
                            measured,
                            worst_nearest_point: nearest,
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
    best.map(|b| PcurveFidelity {
        measured,
        worst_nearest_point: worst_nearest,
        ..b
    })
}

#[cfg(test)]
#[path = "pcurve_fidelity_tests.rs"]
mod tests;
