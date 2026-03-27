//! Pcurve recomputation for BRep faces.
//!
//! Projects edge curve endpoints (or vertex positions) to surface UV space
//! and creates degree-1 NURBS pcurves for half-edges that lack them.

use super::FaceHandle;
use super::brep::BRep;

impl BRep {
    /// Recompute pcurves for all half-edges on a face by projecting
    /// edge curve endpoints to the face surface's UV space.
    ///
    /// Skips half-edges that already have pcurves. Falls back to vertex
    /// positions for edges that lack curve geometry. Clamps UV coordinates
    /// to the surface's finite parameter range.
    pub fn recompute_face_pcurves(
        &mut self,
        face_handle: FaceHandle,
    ) -> Result<(), crate::error::CoreError> {
        use crate::geometry::curve::CurveEval;
        use crate::geometry::surface::SurfaceEval;
        use crate::geometry::uv_projection::project_point_to_uv;
        use crate::geometry::{NurbsCurve, Point3};

        let face = self
            .faces
            .get(face_handle)
            .ok_or_else(|| crate::error::CoreError::InvalidHandle("face not found".into()))?;
        let surface = match &face.surface {
            Some(sof) => sof.geometry.clone(),
            None => return Ok(()),
        };

        let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();

        // Collect all loop handles for this face.
        let outer_loop = face.outer_loop;
        let inner_loops = face.inner_loops.clone();

        let mut all_loop_handles = vec![outer_loop];
        all_loop_handles.extend(inner_loops);

        // Collect half-edge handles from all loops before mutating.
        let mut he_handles = Vec::new();
        for lh in &all_loop_handles {
            let loop_data = match self.loops.get(*lh) {
                Some(l) => l,
                None => continue,
            };
            let start_he = loop_data.half_edge;
            if !self.half_edges.contains_key(start_he) {
                continue;
            }
            he_handles.push(start_he);
            let mut current = self.half_edges[start_he].next;
            while current != start_he {
                he_handles.push(current);
                current = self.half_edges[current].next;
            }
        }

        // For each half-edge without a pcurve, compute one.
        for heh in he_handles {
            if self.half_edges[heh].pcurve.is_some() {
                continue;
            }

            let he_origin_pos = self.vertices[self.half_edges[heh].origin].position;
            let next_heh = self.half_edges[heh].next;
            let he_dest_pos = self.vertices[self.half_edges[next_heh].origin].position;

            let edge_handle = self.half_edges[heh].edge;
            let (start_3d, end_3d) = if let Some(ref curve_on_edge) = self.edges[edge_handle].curve
            {
                // Determine direction: compare he.origin position with curve start.
                let curve_start_pos = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
                let curve_end_pos = curve_on_edge.geometry.point_at(curve_on_edge.t_end);

                if he_origin_pos.distance_to(&curve_start_pos)
                    <= he_origin_pos.distance_to(&curve_end_pos)
                {
                    (curve_start_pos, curve_end_pos)
                } else {
                    (curve_end_pos, curve_start_pos)
                }
            } else {
                // No edge curve: fall back to vertex positions.
                (he_origin_pos, he_dest_pos)
            };

            // Clamp UV to surface parameter range (if finite).
            let clamp_if_finite = |val: f64, lo: f64, hi: f64| -> f64 {
                if lo.is_finite() && hi.is_finite() {
                    val.clamp(lo, hi)
                } else {
                    val
                }
            };

            // For planar surfaces, degree-1 (endpoint-to-endpoint) is exact.
            // For curved surfaces, sample the edge curve at intermediate
            // points and fit a higher-degree pcurve to capture the UV-space
            // curvature.
            let is_planar = matches!(
                surface,
                crate::geometry::surface::SurfaceGeometry::Plane { .. }
            );
            let edge_curve = &self.edges[edge_handle].curve;

            let pcurve = if is_planar || edge_curve.is_none() {
                // Degree-1: straight line in UV.
                let (u_s, v_s) = match project_point_to_uv(&surface, &start_3d) {
                    Some(uv) => uv,
                    None => continue,
                };
                let (u_e, v_e) = match project_point_to_uv(&surface, &end_3d) {
                    Some(uv) => uv,
                    None => continue,
                };
                NurbsCurve {
                    degree: 1,
                    knots: vec![0.0, 0.0, 1.0, 1.0],
                    control_points: vec![
                        Point3::new(
                            clamp_if_finite(u_s, u_min, u_max),
                            clamp_if_finite(v_s, v_min, v_max),
                            0.0,
                        ),
                        Point3::new(
                            clamp_if_finite(u_e, u_min, u_max),
                            clamp_if_finite(v_e, v_min, v_max),
                            0.0,
                        ),
                    ],
                    weights: vec![1.0, 1.0],
                }
            } else {
                // Multi-point sampling for curved surfaces.
                let Some(curve_on_edge) = edge_curve.as_ref() else {
                    continue;
                };
                let forward = he_origin_pos
                    .distance_to(&curve_on_edge.geometry.point_at(curve_on_edge.t_start))
                    <= he_origin_pos
                        .distance_to(&curve_on_edge.geometry.point_at(curve_on_edge.t_end));
                let (t_lo, t_hi) = if forward {
                    (curve_on_edge.t_start, curve_on_edge.t_end)
                } else {
                    (curve_on_edge.t_end, curve_on_edge.t_start)
                };

                // Dense sampling + interpolation for accurate torus/cylinder
                // pcurves. Use enough samples for the curve to be
                // indistinguishable from exact at STEP export tolerance.
                let n_samples = 64usize;
                let mut uv_points = Vec::with_capacity(n_samples);
                let mut params = Vec::with_capacity(n_samples);
                for i in 0..n_samples {
                    let t = t_lo + (t_hi - t_lo) * (i as f64 / (n_samples - 1).max(1) as f64);
                    let p = curve_on_edge.geometry.point_at(t);
                    if let Some((u, v)) = project_point_to_uv(&surface, &p) {
                        uv_points.push(Point3::new(
                            clamp_if_finite(u, u_min, u_max),
                            clamp_if_finite(v, v_min, v_max),
                            0.0,
                        ));
                        params.push(i as f64 / (n_samples - 1).max(1) as f64);
                    }
                }

                if uv_points.len() >= 4 {
                    // Dense sampling + interpolating fit for accurate
                    // torus/cylinder pcurves. Use n_cp close to n_samples
                    // so the curve passes near all sample points.
                    let degree = 3.min(uv_points.len() - 1);
                    let n_cp = uv_points.len().min(50);
                    match crate::geometry::fitting::fit_curve(&uv_points, &params, degree, n_cp) {
                        Ok(mut fitted) => {
                            for cp in &mut fitted.control_points {
                                cp.x = clamp_if_finite(cp.x, u_min, u_max);
                                cp.y = clamp_if_finite(cp.y, v_min, v_max);
                            }
                            fitted
                        }
                        Err(_) => NurbsCurve {
                            degree: 1,
                            knots: vec![0.0, 0.0, 1.0, 1.0],
                            control_points: vec![uv_points[0], uv_points[uv_points.len() - 1]],
                            weights: vec![1.0, 1.0],
                        },
                    }
                } else if uv_points.len() >= 2 {
                    NurbsCurve {
                        degree: 1,
                        knots: vec![0.0, 0.0, 1.0, 1.0],
                        control_points: vec![uv_points[0], uv_points[uv_points.len() - 1]],
                        weights: vec![1.0, 1.0],
                    }
                } else {
                    continue;
                }
            };

            self.half_edges[heh].pcurve = Some(pcurve);
        }

        Ok(())
    }

    /// Recompute pcurves for ALL faces in the BRep.
    pub fn recompute_all_pcurves(&mut self) -> Result<(), crate::error::CoreError> {
        let faces: Vec<_> = self.faces.keys().collect();
        for fh in faces {
            self.recompute_face_pcurves(fh)?;
        }
        Ok(())
    }

    /// Detect and replace stale pcurves whose endpoints no longer match
    /// the actual edge curve positions.
    ///
    /// A pcurve is stale when the 3D position at pcurve(t=0) projected
    /// to the face surface doesn't match the half-edge origin vertex
    /// position within tolerance.  This happens when the fillet or other
    /// operations move edge positions without updating existing pcurves.
    pub fn fix_stale_pcurves(&mut self) {
        use crate::geometry::curve::CurveEval;
        use crate::geometry::surface::SurfaceEval;
        use crate::geometry::uv_projection::project_point_to_uv;

        let he_keys: Vec<_> = self.half_edges.keys().collect();
        for heh in he_keys {
            let pcurve = match self.half_edges[heh].pcurve.as_ref() {
                Some(pc) => pc.clone(),
                None => continue,
            };

            // Get the face surface for this half-edge.
            let loop_handle = self.half_edges[heh].loop_;
            let face_handle = match self.loops.get(loop_handle) {
                Some(l) => l.face,
                None => continue,
            };
            let surface = match self.faces.get(face_handle).and_then(|f| f.surface.as_ref()) {
                Some(s) => s.geometry.clone(),
                None => continue,
            };

            // Check: evaluate pcurve at start → project that UV to 3D →
            // compare with the actual half-edge origin vertex position.
            let (t0, _t1) = pcurve.parameter_range();
            let uv_start = pcurve.point_at(t0);
            let pos_from_pcurve = surface.point_at(uv_start.x, uv_start.y);
            let actual_pos = self.vertices[self.half_edges[heh].origin].position;

            let error = pos_from_pcurve.distance_to(&actual_pos);
            if error > 0.05 {
                // Pcurve is stale. Clear it so recompute_face_pcurves
                // will regenerate it on the next pass.
                self.half_edges[heh].pcurve = None;
            }
        }

        // Recompute cleared pcurves.
        let _ = self.recompute_all_pcurves();
    }
}

#[cfg(test)]
mod tests {
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn recompute_pcurves_on_cube_is_idempotent() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        // All half-edges already have pcurves from unit_cube.
        brep.recompute_all_pcurves().unwrap();
        for (_, he) in &brep.half_edges {
            assert!(he.pcurve.is_some(), "pcurve should still be present");
        }
    }

    #[test]
    fn recompute_fills_missing_pcurves() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        // Strip all pcurves.
        let he_keys: Vec<_> = brep.half_edges.keys().collect();
        for heh in &he_keys {
            brep.half_edges[*heh].pcurve = None;
        }
        // Recompute should fill them back.
        brep.recompute_all_pcurves().unwrap();
        for (_, he) in &brep.half_edges {
            assert!(he.pcurve.is_some(), "pcurve should be recomputed");
        }
    }
}
