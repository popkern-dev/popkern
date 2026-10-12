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

    /// Check 9 over a chosen edge set (collect-all form).
    ///
    /// Same measurement and error shape as
    /// [`BRep::validate_curve_surface_all`], restricted to `edges`
    /// (handles no longer in the BRep are skipped). For postcondition
    /// gates that own a known territory — a blend corner's face star —
    /// and must not pay, or be failed by, a whole-body sweep.
    pub fn validate_curve_surface_edges(
        &self,
        edges: &[crate::topology::EdgeHandle],
    ) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        for &edge_handle in edges {
            let Some(edge) = self.edges.get(edge_handle) else {
                continue;
            };
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
            // the band). Each sample accepts the best of the cold
            // projection, a continuation seed from the previous
            // sample's witness, and a seed from THIS side's stored
            // PCURVE — the authoritative uv chart of the edge on
            // this face, immune to the wrong-basin clamp (measured,
            // PK-472: a 5.8e-4 graze-junction edge on the m6 nut's
            // mouth cone read 7.9e-5 deviating under cold +
            // continuation while its true residual was 1.2e-8; the
            // otherwise-valid boolean declined). The pcurve's own
            // orientation relative to the edge curve is unknown
            // here, so both fractions seed. A real defect has no
            // nearby witness under any seed.
            let pcurve = self.half_edges.get(side).and_then(|d| d.pcurve.as_ref());
            let mut prev_uv: Option<(f64, f64)> = None;
            for fraction in SAMPLE_FRACTIONS {
                let t = curve.t_start + (curve.t_end - curve.t_start) * fraction;
                let point = curve.geometry.point_at(t);
                let seeded = prev_uv.and_then(|(pu, pv)| {
                    crate::geometry::uv_projection::project_point_to_uv_near(
                        &surface.geometry,
                        &point,
                        pu,
                        pv,
                    )
                });
                let pcurve_seeded = pcurve.iter().flat_map(|pc| {
                    [fraction, 1.0 - fraction].into_iter().filter_map(|f| {
                        let q = pc.point_at(pc.t_start + (pc.t_end - pc.t_start) * f);
                        crate::geometry::uv_projection::project_point_to_uv_near(
                            &surface.geometry,
                            &point,
                            q.x,
                            q.y,
                        )
                    })
                });
                let witness = |uv: (f64, f64)| {
                    (
                        point.distance_to(&surface.geometry.point_at(uv.0, uv.1)),
                        uv,
                    )
                };
                let nearer = |a: &(f64, (f64, f64)), b: &(f64, (f64, f64))| {
                    a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)
                };
                // The seeded witnesses first: a sample one of them
                // places within the band passes whatever the cold
                // search finds, and the cold multisheet search was the
                // cost of every boolean's and fastener's check 9 over
                // a tapped housing's earlier threads (PK-744).
                let warm = seeded
                    .into_iter()
                    .chain(pcurve_seeded)
                    .map(witness)
                    .min_by(nearer);
                let best = match warm {
                    Some(w) if w.0 <= tolerance => Some(w),
                    _ => project_point_to_uv_multisheet(&surface.geometry, &point)
                        .map(witness)
                        .into_iter()
                        .chain(warm)
                        .min_by(nearer),
                };
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
                    // Knot-seeded multi-start rescue (PK-472): on a
                    // multi-turn NURBS chart the cold multisheet AND
                    // continuation witnesses can BOTH converge into
                    // an adjacent turn's basin (measured: a 5.8e-4
                    // graze-junction stub on the m6 mouth chamfer
                    // read 7.9e-5 deviating while the knot-span
                    // seeds verify 1.2e-8) — the same rescue the
                    // PK-463 closest-point campaign added to the
                    // kernel's projection.
                    if fine_dev > tolerance
                        && let crate::geometry::surface::SurfaceGeometry::Nurbs(nurbs) =
                            &surface.geometry
                        && let Some(seeds) = nurbs.projection_seeds()
                        && let Some((u, v)) = crate::geometry::surface::closest_point_from_seeds(
                            nurbs,
                            &point,
                            &seeds,
                            (false, false),
                            6,
                        )
                    {
                        let d = point.distance_to(&surface.geometry.point_at(u, v));
                        if d < fine_dev {
                            fine_dev = d;
                        }
                    }
                    if fine_dev > tolerance {
                        if popkern_trace_check9() {
                            let ends = {
                                let p0 = curve.geometry.point_at(curve.t_start);
                                let p1 = curve.geometry.point_at(curve.t_end);
                                format!(
                                    "({:.6},{:.6},{:.6})->({:.6},{:.6},{:.6}) chord {:.3e}",
                                    p0.x,
                                    p0.y,
                                    p0.z,
                                    p1.x,
                                    p1.y,
                                    p1.z,
                                    p0.distance_to(&p1)
                                )
                            };
                            eprintln!(
                                "[check9] side {side:?} face {face:?} ({}) t {t:.3} point ({:.6},{:.6},{:.6}) dev {fine_dev:.3e} trange ({:.3},{:.3}) ends {ends}",
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
