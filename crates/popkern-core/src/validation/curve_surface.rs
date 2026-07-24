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

use crate::error::ValidationError;
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::uv_projection::project_point_to_uv;
use crate::topology::BRep;

/// Curve stations checked per (edge, face) pairing: both endpoints
/// plus three interior points — enough to catch off-surface chords
/// and wrong-branch arcs without measurable validate() cost.
const SAMPLE_FRACTIONS: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];

impl BRep {
    /// Check 9: Curve-surface consistency.
    ///
    /// Verifies that edge curves lie on the surfaces of their
    /// adjacent faces, within the edge's explicit tolerance band.
    /// Pairings without geometry on either side pass; pairings whose
    /// surface kind offers no closest-point projection are skipped.
    pub fn validate_curve_surface(&self) -> Result<(), ValidationError> {
        for (edge_handle, edge) in &self.edges {
            let Some(curve) = edge.curve.as_ref() else {
                continue;
            };
            let mut sides = Vec::with_capacity(2);
            if let Some(half_edge) = self.half_edges.get(edge.half_edge) {
                sides.push(edge.half_edge);
                if self.half_edges.contains_key(half_edge.twin) {
                    sides.push(half_edge.twin);
                }
            }
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
                    let cold = project_point_to_uv(&surface.geometry, &point);
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
                            let cold_f = project_point_to_uv(&surface.geometry, &fp);
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
                            return Err(ValidationError::CurveSurfaceDeviation {
                                edge: edge_handle,
                                max_deviation: fine_dev,
                            });
                        }
                        prev_uv = fine_uv;
                    }
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    fn make_cube() -> BRep {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev3 = brep
            .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
        brep.extrude(mef.face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
            .unwrap();
        brep
    }

    #[test]
    fn topology_only_cube_passes_curve_surface() {
        let brep = make_cube();
        assert!(brep.validate_curve_surface().is_ok());
    }

    #[test]
    fn geometric_cube_passes_curve_surface() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        assert!(brep.validate_curve_surface().is_ok());
    }

    #[test]
    fn off_surface_edge_curve_is_rejected() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        // Replace one edge's curve with a chord displaced off every
        // adjacent face plane.
        let edge = brep.edges.keys().next().unwrap();
        let half_edge = brep.edges[edge].half_edge;
        let start = brep.vertices[brep.half_edges[half_edge].origin].position;
        let twin = brep.half_edges[half_edge].twin;
        let end = brep.vertices[brep.half_edges[twin].origin].position;
        let lifted = Point3::new(
            0.5 * (start.x + end.x) + 0.3,
            0.5 * (start.y + end.y) + 0.3,
            0.5 * (start.z + end.z) + 0.3,
        );
        brep.edges[edge].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(crate::geometry::NurbsCurve {
                degree: 2,
                knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                control_points: vec![start, lifted, end],
                weights: vec![1.0, 1.0, 1.0],
            }),
            t_start: 0.0,
            t_end: 1.0,
        });
        let result = brep.validate_curve_surface();
        assert!(
            result.is_err(),
            "an edge curve bulging off its face planes must fail check 9"
        );
    }
}
