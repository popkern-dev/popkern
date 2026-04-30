//! Pcurve recomputation for BRep faces.
//!
//! Projects edge curve endpoints (or vertex positions) to surface UV space
//! and creates degree-1 NURBS pcurves for half-edges that lack them.

use super::FaceHandle;
use super::brep::BRep;
use crate::geometry::curve::CurveEval;

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
        use crate::geometry::surface::SurfaceEval;
        use crate::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};
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

        let outer_loop = face.outer_loop;
        let inner_loops = face.inner_loops.clone();

        let mut all_loop_handles = vec![outer_loop];
        all_loop_handles.extend(inner_loops);

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

        // For each half-edge without a pcurve (or with a degenerate one
        // on a periodic surface), compute a fresh pcurve.
        let is_periodic = {
            let (pu, pv) = surface.periodicity();
            pu.is_some() || pv.is_some()
        };

        for heh in he_handles {
            if let Some(ref existing) = self.half_edges[heh].pcurve {
                // On periodic surfaces, clear degenerate pcurves (UV span ~0
                // but edge has a real curve) — e.g. self-loop circles on
                // a cylinder. Non-periodic self-loops use 3D trim fallback.
                if !is_periodic {
                    continue;
                }
                let (t0, t1) = existing.parameter_range();
                let s = existing.point_at(t0);
                let e = existing.point_at(t1);
                let uv_span = ((e.x - s.x).powi(2) + (e.y - s.y).powi(2)).sqrt();
                let edge_h = self.half_edges[heh].edge;
                let long_curve = self.edges[edge_h]
                    .curve
                    .as_ref()
                    .is_some_and(|c| (c.t_end - c.t_start).abs() > 1e-6);
                if uv_span > 1e-6 || !long_curve {
                    continue;
                }
                self.half_edges[heh].pcurve = None;
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

                // Quick linearity check: if the UV midpoint lies on the
                // line between start/end UV, the pcurve is a straight
                // line in UV space.  Use degree-1 (exact) instead of the
                // 64-point degree-3 fit, avoiding control-point overshoot.
                let (u_s, v_s) = match project_point_to_uv(&surface, &start_3d) {
                    Some(uv) => uv,
                    None => continue,
                };
                let (u_e, v_e) = match project_point_to_uv(&surface, &end_3d) {
                    Some(uv) => uv,
                    None => continue,
                };
                let t_mid = (t_lo + t_hi) * 0.5;
                let p_mid = curve_on_edge.geometry.point_at(t_mid);
                let uv_linear = if let Some((u_mid, v_mid)) = project_point_to_uv(&surface, &p_mid)
                {
                    let expected_u = (u_s + u_e) * 0.5;
                    let expected_v = (v_s + v_e) * 0.5;
                    let deviation =
                        ((u_mid - expected_u).powi(2) + (v_mid - expected_v).powi(2)).sqrt();
                    deviation < 1e-4
                } else {
                    false
                };

                if uv_linear {
                    // UV path is linear -- exact degree-1 pcurve.
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
                    // UV path is curved — dense sampling.
                    // OCCT-style pcurve computation (per ShapeConstruct_
                    // ProjectCurveOnSurface): sample 3D curve, project
                    // each point to UV using CONTINUATION from the
                    // previous UV (not independent projection), unwrap
                    // seam crossings, then INTERPOLATE (not fit).
                    let n_samples = 23usize; // OCCT default
                    let mut uv_points: Vec<Point3> = Vec::with_capacity(n_samples);
                    let (period_u, period_v) = surface.periodicity();
                    let mut prev_uv: Option<(f64, f64)> = None;

                    for i in 0..n_samples {
                        let frac = i as f64 / (n_samples - 1).max(1) as f64;
                        let t = t_lo + (t_hi - t_lo) * frac;
                        let p = curve_on_edge.geometry.point_at(t);

                        // Continuation: use previous UV as starting
                        // guess for Gauss-Newton. This prevents branch
                        // jumping on multi-valued inverse maps.
                        let uv = if let Some((pu, pv)) = prev_uv {
                            project_point_to_uv_near(&surface, &p, pu, pv)
                        } else {
                            project_point_to_uv(&surface, &p)
                        };

                        if let Some((mut u, mut v)) = uv {
                            // Seam unwrap relative to PREVIOUS point
                            // (OCCT style: half-period threshold).
                            if let (Some(pu), Some((prev_u, _))) = (period_u, prev_uv) {
                                while u - prev_u > pu * 0.5 {
                                    u -= pu;
                                }
                                while u - prev_u < -pu * 0.5 {
                                    u += pu;
                                }
                            }
                            if let (Some(pv), Some((_, prev_v))) = (period_v, prev_uv) {
                                while v - prev_v > pv * 0.5 {
                                    v -= pv;
                                }
                                while v - prev_v < -pv * 0.5 {
                                    v += pv;
                                }
                            }
                            prev_uv = Some((u, v));
                            uv_points.push(Point3::new(u, v, 0.0));
                        }
                    }

                    if uv_points.len() >= 4 {
                        // Interpolation (not approximation): curve
                        // passes through every sample point exactly.
                        // Use degree 3 with n_cp = n_samples. With
                        // continuation-based projection (no branch
                        // jumping), the points are smooth enough that
                        // interpolation doesn't oscillate.
                        let params: Vec<f64> = (0..uv_points.len())
                            .map(|i| i as f64 / (uv_points.len() - 1).max(1) as f64)
                            .collect();
                        let degree = 3.min(uv_points.len() - 1);
                        let n_cp = uv_points.len(); // interpolation
                        match crate::geometry::fitting::fit_curve(&uv_points, &params, degree, n_cp)
                        {
                            Ok(fitted) => fitted,
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
                }
            };

            self.half_edges[heh].pcurve = Some(pcurve.into());
        }

        self.fix_pcurve_endpoints(face_handle);

        Ok(())
    }

    /// Force pcurve endpoints to match vertex UV projections exactly.
    ///
    /// For each half-edge with a pcurve, verify that the pcurve evaluated
    /// at t=0 produces the correct UV for he.origin, and the pcurve at t=1
    /// produces the correct UV for the next half-edge's origin. Adjust
    /// the first/last control points if needed.
    ///
    /// On periodic surfaces (cylinder u, torus u/v), the projected UV may
    /// differ from the pcurve endpoint by a full period. This is correct
    /// for edges that span the seam (e.g., self-loop circles), so the fix
    /// is skipped when the error is close to a multiple of the period.
    fn fix_pcurve_endpoints(&mut self, face_handle: FaceHandle) {
        use crate::geometry::uv_projection::project_point_to_uv;

        let surface = match self.faces.get(face_handle).and_then(|f| f.surface.as_ref()) {
            Some(s) => s.geometry.clone(),
            None => return,
        };
        let (period_u, period_v) = surface.periodicity();

        let face = match self.faces.get(face_handle) {
            Some(f) => f,
            None => return,
        };
        let outer_loop = face.outer_loop;
        let inner_loops = face.inner_loops.clone();
        let mut all_loop_handles = vec![outer_loop];
        all_loop_handles.extend(inner_loops);

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

        let tol = 1e-10;

        for heh in he_handles {
            let pcurve = match self.half_edges[heh].pcurve.as_ref() {
                Some(pc) => pc.clone(),
                None => continue,
            };

            // §6c-D Stage 4.6: analytic pcurves are exact-by-construction;
            // skip endpoint surgery for them.  Closest-point projection on
            // Gregory caps frequently lands on the wrong boundary.
            use crate::geometry::pcurve::PcurveGeometry;
            if !matches!(pcurve.geometry, PcurveGeometry::Nurbs(_)) {
                continue;
            }

            let origin_pos = self.vertices[self.half_edges[heh].origin].position;
            let next_heh = self.half_edges[heh].next;
            let dest_pos = self.vertices[self.half_edges[next_heh].origin].position;

            let exact_uv_start = match project_point_to_uv(&surface, &origin_pos) {
                Some((u, v)) => crate::geometry::Point3::new(u, v, 0.0),
                None => continue,
            };
            let exact_uv_end = match project_point_to_uv(&surface, &dest_pos) {
                Some((u, v)) => crate::geometry::Point3::new(u, v, 0.0),
                None => continue,
            };

            let (t0, t1) = pcurve.parameter_range();
            let eval_start = pcurve.point_at(t0);
            let eval_end = pcurve.point_at(t1);

            // Compute UV error modulo periodicity: if the difference is
            // close to a multiple of the period, the pcurve is correct
            // (it wraps around the seam) and should not be "fixed".
            let periodic_err = |exact: f64, eval: f64, period: Option<f64>| -> f64 {
                let diff = (exact - eval).abs();
                if let Some(p) = period.filter(|&p| p > 0.0) {
                    let remainder = diff % p;
                    return remainder.min(p - remainder);
                }
                diff
            };

            let start_err_u = periodic_err(exact_uv_start.x, eval_start.x, period_u);
            let start_err_v = periodic_err(exact_uv_start.y, eval_start.y, period_v);
            let start_err = (start_err_u.powi(2) + start_err_v.powi(2)).sqrt();

            let end_err_u = periodic_err(exact_uv_end.x, eval_end.x, period_u);
            let end_err_v = periodic_err(exact_uv_end.y, eval_end.y, period_v);
            let end_err = (end_err_u.powi(2) + end_err_v.powi(2)).sqrt();

            if start_err < tol && end_err < tol {
                continue;
            }

            // Pole-aware bypass: at a sphere pole, `project_point_to_uv`
            // returns a `u` from `atan2(0,0)` (effectively undefined),
            // so `start_err`/`end_err` can be huge even when the
            // existing pcurve is geometrically correct.  Verify in
            // 3-D: if `surface.point_at(eval_start)` matches
            // `origin_pos` within tolerance, the pcurve is fine and
            // converting it to NURBS just to overwrite the CP with a
            // noisy `u` would *introduce* error.
            use crate::geometry::surface::SurfaceEval;
            let eval_start_3d = surface.point_at(eval_start.x, eval_start.y);
            let eval_end_3d = surface.point_at(eval_end.x, eval_end.y);
            let pos_tol = 1e-9_f64;
            let start_3d_ok = eval_start_3d.distance_to(&origin_pos) < pos_tol;
            let end_3d_ok = eval_end_3d.distance_to(&dest_pos) < pos_tol;
            if start_3d_ok && end_3d_ok {
                continue;
            }

            // Materialise the pcurve as NURBS for endpoint surgery;
            // analytic variants (Linear) round-trip through `to_nurbs`
            // so the same code path serves both.
            let mut fixed = pcurve.to_nurbs();
            let n_cp = fixed.control_points.len();

            if fixed.degree == 1 && n_cp == 2 {
                // Degree-1: control points ARE the evaluated endpoints.
                if start_err >= tol {
                    fixed.control_points[0] = exact_uv_start;
                }
                if end_err >= tol {
                    fixed.control_points[n_cp - 1] = exact_uv_end;
                }
            } else {
                // Higher degree: shift first/last control points by the delta.
                if start_err >= tol {
                    let dx = exact_uv_start.x - eval_start.x;
                    let dy = exact_uv_start.y - eval_start.y;
                    fixed.control_points[0].x += dx;
                    fixed.control_points[0].y += dy;
                }
                if end_err >= tol {
                    let dx = exact_uv_end.x - eval_end.x;
                    let dy = exact_uv_end.y - eval_end.y;
                    fixed.control_points[n_cp - 1].x += dx;
                    fixed.control_points[n_cp - 1].y += dy;
                }
            }

            self.half_edges[heh].pcurve = Some(fixed.into());
        }
    }

    /// Recompute pcurves for ALL faces in the BRep.
    pub fn recompute_all_pcurves(&mut self) -> Result<(), crate::error::CoreError> {
        let faces: Vec<_> = self.faces.keys().collect();
        for fh in faces {
            self.recompute_face_pcurves(fh)?;
        }
        Ok(())
    }

    /// Like [`fix_stale_pcurves`] but only processes half-edges on
    /// Detect and replace stale pcurves whose endpoints no longer match
    /// the actual edge curve positions.
    ///
    /// A pcurve is stale when the 3D position at pcurve(t=0) projected
    /// to the face surface doesn't match the half-edge origin vertex
    /// position within tolerance.  This happens when the fillet or other
    /// operations move edge positions without updating existing pcurves.
    pub fn fix_stale_pcurves(&mut self) {
        use crate::geometry::pcurve::PcurveGeometry;
        use crate::geometry::surface::SurfaceEval;

        let he_keys: Vec<_> = self.half_edges.keys().collect();
        for heh in he_keys {
            let pcurve = match self.half_edges[heh].pcurve.as_ref() {
                Some(pc) => pc.clone(),
                None => continue,
            };

            // §6c-D Stage 4.6: analytic pcurves are authoritative;
            // run the staleness check only on Nurbs pcurves.
            if !matches!(pcurve.geometry, PcurveGeometry::Nurbs(_)) {
                continue;
            }

            let loop_handle = self.half_edges[heh].loop_;
            let face_handle = match self.loops.get(loop_handle) {
                Some(l) => l.face,
                None => continue,
            };
            let surface = match self.faces.get(face_handle).and_then(|f| f.surface.as_ref()) {
                Some(s) => s.geometry.clone(),
                None => continue,
            };

            let (t0, _t1) = pcurve.parameter_range();
            let uv_start = pcurve.point_at(t0);
            let pos_from_pcurve = surface.point_at(uv_start.x, uv_start.y);
            let actual_pos = self.vertices[self.half_edges[heh].origin].position;

            let error = pos_from_pcurve.distance_to(&actual_pos);
            if error > 0.05 {
                self.half_edges[heh].pcurve = None;
            }
        }

        let _ = self.recompute_all_pcurves();

        let face_keys: Vec<_> = self.faces.keys().collect();
        for fh in face_keys {
            self.fix_pcurve_endpoints(fh);
        }
    }
}

#[cfg(test)]
#[path = "pcurve_recompute_tests.rs"]
mod tests;
