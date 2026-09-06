//! Check 9: Curve-surface consistency — edge curves must lie on the
//! surfaces of their adjacent faces.
//!
//! For every edge that carries a curve, each adjacent face that
//! carries a surface must contain the curve: the curve is sampled at
//! its endpoints and interior stations and each sample's distance to
//! the surface must stay inside the edge's tolerance band. An edge
//! curve off its face surface is a lie every downstream consumer
//! inherits (trims, tessellation, STEP export), and it is exactly the
//! class external validators (OCC BRepCheck) reject while this check
//! was a stub (PK-165).
//!
//! Tolerance is explicit and local: the edge's own tolerance when
//! set, else the adjacent face's, else the model tolerance — widened
//! by the kernel's shared 10× repair/weld band. Samples that cannot
//! be projected onto the surface (no closest-point solve for the
//! surface kind) are skipped: an unverifiable pairing is not treated
//! as a failure, and the projection-capable kinds (all analytic
//! surfaces and NURBS) cover the production corpus.
//!
//! Two entry points (PK-235): [`BRep::validate_curve_surface`] is the
//! early-exit form for boolean gates, and
//! [`BRep::validate_curve_surface_all`] reports EVERY failing edge
//! with its per-edge worst deviation. `validate_all` uses the
//! collect-all form — the early-exit form once fed it, which meant a
//! model with a dozen broken edges reported one error carrying the
//! FIRST deviation found, not the worst (measured on the PK-231
//! saddle fixture: reported 4.685e-4 while the true worst edge sat
//! at ~1.3e-1).

use crate::error::ValidationError;

/// The `check9` trace topic: name the failing (edge side, face)
/// pairing and sample when the deviation check fires — the error
/// type carries only the edge, and locating the face cost a probe
/// during PK-427's chamfer diagnosis.
fn popkern_trace_check9() -> bool {
    crate::trace::trace_enabled("check9")
}
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::uv_projection::project_point_to_uv_multisheet;
use crate::topology::{BRep, Edge};

/// Curve stations checked per (edge, face) pairing: both endpoints
/// plus three interior points — enough to catch off-surface chords
/// and wrong-branch arcs without measurable validate() cost.
const SAMPLE_FRACTIONS: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];

impl BRep {
    /// Check 9: Curve-surface consistency (early-exit form).
    ///
    /// Verifies that edge curves lie on the surfaces of their
    /// adjacent faces, within the edge's explicit tolerance band.
    /// Pairings without geometry on either side pass; pairings whose
    /// surface kind offers no closest-point projection are skipped.
    ///
    /// Returns the FIRST violation found — use
    /// [`BRep::validate_curve_surface_all`] when the full failing
    /// set (or the true worst deviation) matters.
    pub fn validate_curve_surface(&self) -> Result<(), ValidationError> {
        for (edge_handle, edge) in &self.edges {
            if let Some(deviation) = self.curve_surface_edge_violation(edge, true) {
                return Err(ValidationError::CurveSurfaceDeviation {
                    edge: edge_handle,
                    max_deviation: deviation,
                });
            }
        }
        Ok(())
    }

    /// Check 9: Curve-surface consistency (collect-all form).
    ///
    /// Same measurement as [`BRep::validate_curve_surface`], but
    /// visits every edge and reports one
    /// [`ValidationError::CurveSurfaceDeviation`] per failing edge,
    /// carrying that edge's worst deviation across both adjacent
    /// faces and all samples. This is the form `validate_all` uses:
    /// acceptance gates that read `max_deviation` get the true
    /// per-edge worst, and the number of errors is the number of
    /// broken edges.
    pub fn validate_curve_surface_all(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        for (edge_handle, edge) in &self.edges {
            if let Some(deviation) = self.curve_surface_edge_violation(edge, false) {
                errors.push(ValidationError::CurveSurfaceDeviation {
                    edge: edge_handle,
                    max_deviation: deviation,
                });
            }
        }
        errors
    }

    /// Measure one edge against both adjacent face surfaces.
    ///
    /// Returns `Some(worst_over_band_deviation)` when any sample
    /// exceeds the edge's tolerance band after the refinement retry,
    /// `None` when the edge passes (or is unverifiable). With
    /// `stop_at_first` the scan returns at the first confirmed
    /// violation (early-exit form); otherwise it keeps scanning and
    /// returns the per-edge maximum.
    fn curve_surface_edge_violation(&self, edge: &Edge, stop_at_first: bool) -> Option<f64> {
        let curve = edge.curve.as_ref()?;
        let mut sides = Vec::with_capacity(2);
        if let Some(half_edge) = self.half_edges.get(edge.half_edge) {
            sides.push(edge.half_edge);
            if self.half_edges.contains_key(half_edge.twin) {
                sides.push(half_edge.twin);
            }
        }
        let mut worst: Option<f64> = None;
        for side in sides {
            let Some(face) = self
                .half_edges
                .get(side)
                .and_then(|data| self.loops.get(data.loop_))
                .map(|loop_data| loop_data.face)
            else {
                continue;
            };
            let Some(face_data) = self.faces.get(face) else {
                continue;
            };
            let Some(surface) = face_data.surface.as_ref() else {
                continue;
            };
            let tolerance = edge
                .tolerance
                .or(face_data.tolerance)
                .unwrap_or_else(|| self.tolerance.linear())
                .max(1.0e-9)
                * 10.0;
            // Witness-based sampling: the check verifies the
            // sample's distance TO THE SURFACE, and the cold
            // global projection is only a heuristic minimizer —
            // on extended or periodic charts it can clamp to the
            // wrong chart region and misreport an on-surface
            // sample as deviating (measured 2.07e-5 on a boss
            // corner-adjacent edge whose true residual is below
            // the band). Each sample accepts the best of the
            // cold projection and a continuation seed from the
            // previous sample's witness; a real defect has no
            // nearby witness either way.
            let mut prev_uv: Option<(f64, f64)> = None;
            for fraction in SAMPLE_FRACTIONS {
                let t = curve.t_start + (curve.t_end - curve.t_start) * fraction;
                let point = curve.geometry.point_at(t);
                let cold = project_point_to_uv_multisheet(&surface.geometry, &point);
                let seeded = prev_uv.and_then(|(pu, pv)| {
                    crate::geometry::uv_projection::project_point_to_uv_near(
                        &surface.geometry,
                        &point,
                        pu,
                        pv,
                    )
                });
                let best = [cold, seeded]
                    .into_iter()
                    .flatten()
                    .map(|(u, v)| (point.distance_to(&surface.geometry.point_at(u, v)), (u, v)))
                    .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                let Some((deviation, uv)) = best else {
                    // No closest-point solve for this surface
                    // kind: the pairing is unverifiable here, not
                    // failed.
                    break;
                };
                prev_uv = Some(uv);
                if deviation > tolerance {
                    // Refinement retry before failing: the
                    // 5-fraction continuation can stall in a
                    // local minimum a finer walk escapes
                    // (measured 2.07e-5 reported for an edge a
                    // 16-step walk verifies at < 1e-6). A real
                    // defect has no witness at any step size.
                    let mut fine_uv: Option<(f64, f64)> = None;
                    let mut fine_dev = deviation;
                    for k in 0..=16 {
                        let ft = curve.t_start + (t - curve.t_start) * (k as f64) / 16.0;
                        let fp = curve.geometry.point_at(ft);
                        let cold_f = project_point_to_uv_multisheet(&surface.geometry, &fp);
                        let seeded_f = fine_uv.and_then(|(pu, pv)| {
                            crate::geometry::uv_projection::project_point_to_uv_near(
                                &surface.geometry,
                                &fp,
                                pu,
                                pv,
                            )
                        });
                        let best_f = [cold_f, seeded_f]
                            .into_iter()
                            .flatten()
                            .map(|(u, v)| {
                                (fp.distance_to(&surface.geometry.point_at(u, v)), (u, v))
                            })
                            .min_by(|a, b| {
                                a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)
                            });
                        if let Some((d, uv_f)) = best_f {
                            fine_uv = Some(uv_f);
                            if k == 16 {
                                fine_dev = d;
                            }
                        }
                    }
                    if fine_dev > tolerance {
                        if popkern_trace_check9() {
                            eprintln!(
                                "[check9] side {side:?} face {face:?} ({}) t {t:.3} point ({:.6},{:.6},{:.6}) dev {fine_dev:.3e} trange ({:.3},{:.3})",
                                surface.geometry.kind_name(),
                                point.x,
                                point.y,
                                point.z,
                                curve.t_start,
                                curve.t_end
                            );
                        }
                        if stop_at_first {
                            return Some(fine_dev);
                        }
                        worst = Some(worst.map_or(fine_dev, |w: f64| w.max(fine_dev)));
                    }
                    prev_uv = fine_uv;
                }
            }
        }
        worst
    }
}

#[cfg(test)]
#[path = "curve_surface_tests.rs"]
mod tests;
