//! Pcurve recomputation for BRep faces.
//!
//! Projects trim vertices to surface UV space and uses edge curves, when
//! present, for interior pcurve shape on non-planar surfaces.

use super::FaceHandle;
use super::brep::BRep;
use crate::geometry::curve::CurveEval;

impl BRep {
    /// Recompute pcurves for all half-edges on a face by projecting trim
    /// vertices to the face surface's UV space.
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
            let edge_direction = if let Some(ref curve_on_edge) = self.edges[edge_handle].curve {
                // Determine direction: compare he.origin position with curve start.
                let curve_start_pos = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
                let curve_end_pos = curve_on_edge.geometry.point_at(curve_on_edge.t_end);

                let closed = curve_start_pos.distance_to(&curve_end_pos) < 1e-9;
                if closed {
                    // Endpoint distances are EQUAL on a closed curve —
                    // the comparison below always picked forward, so
                    // BOTH half-edges of a closed rim were given
                    // same-direction pcurves and a periodic face's
                    // loop wound ±2 periods (untessellatable). Closed
                    // edges use the canonical convention instead: the
                    // edge's own `half_edge` traverses the curve
                    // forward, its twin backward.
                    if self.edges[edge_handle].half_edge == heh {
                        (curve_on_edge.t_start, curve_on_edge.t_end)
                    } else {
                        (curve_on_edge.t_end, curve_on_edge.t_start)
                    }
                } else if he_origin_pos.distance_to(&curve_start_pos)
                    <= he_origin_pos.distance_to(&curve_end_pos)
                {
                    (curve_on_edge.t_start, curve_on_edge.t_end)
                } else {
                    (curve_on_edge.t_end, curve_on_edge.t_start)
                }
            } else {
                (0.0, 1.0)
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

            // The plane fast path is only valid for STRAIGHT edges: a
            // curved edge lying IN a plane (fillet end-runoff arc) has
            // a curved UV image, and a straight pcurve misstates it by
            // the full sagitta (PK-096: 2.4e-2 on runoff arcs — STEP
            // consumers flag the edge and refuse to stitch).
            let edge_is_straight = match edge_curve {
                None => true,
                Some(c) => {
                    let eval: &dyn crate::geometry::curve::CurveEval = &c.geometry;
                    let a = eval.point_at(c.t_start);
                    let b = eval.point_at(c.t_end);
                    let chord = b - a;
                    let len2 = chord.dot(&chord);
                    let mut straight = true;
                    for i in 1..8 {
                        let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 8.0;
                        let p = eval.point_at(t);
                        let ap = p - a;
                        let f = if len2 < 1e-300 {
                            0.0
                        } else {
                            ap.dot(&chord) / len2
                        };
                        let q =
                            Point3::new(a.x + chord.x * f, a.y + chord.y * f, a.z + chord.z * f);
                        if p.distance_to(&q) > 1.0e-9 {
                            straight = false;
                            break;
                        }
                    }
                    straight
                }
            };

            let pcurve = if (is_planar && edge_is_straight) || edge_curve.is_none() {
                // Degree-1: straight line in UV.
                let (u_s, v_s) = match project_point_to_uv(&surface, &he_origin_pos) {
                    Some(uv) => uv,
                    None => continue,
                };
                let (u_e, v_e) = match project_point_to_uv(&surface, &he_dest_pos) {
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
                let (t_lo, t_hi) = edge_direction;

                // Quick linearity check: if the UV midpoint lies on the
                // line between start/end UV, the pcurve is a straight
                // line in UV space.  Use degree-1 (exact) instead of the
                // 64-point degree-3 fit, avoiding control-point overshoot.
                let (u_s, v_s) = match project_point_to_uv(&surface, &he_origin_pos) {
                    Some(uv) => uv,
                    None => continue,
                };
                let (u_e, v_e) = match project_point_to_uv(&surface, &he_dest_pos) {
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
                        let p = if i == 0 {
                            he_origin_pos
                        } else if i + 1 == n_samples {
                            he_dest_pos
                        } else {
                            curve_on_edge.geometry.point_at(t)
                        };

                        // Continuation: use previous UV as starting
                        // guess for Gauss-Newton. This prevents branch
                        // jumping on multi-valued inverse maps. If the
                        // continuation step fails to converge, RESEED
                        // with an independent projection rather than
                        // dropping the sample — dropped samples used to
                        // collapse the pcurve to a straight 2-point
                        // fallback (PK-095: 5e-2 UV error on fillet
                        // strip trims, the dominant tessellation crack).
                        let seed_u = u_s + (u_e - u_s) * frac;
                        let seed_v = v_s + (v_e - v_s) * frac;
                        let uv = if let Some((pu, pv)) = prev_uv {
                            project_point_to_uv_near(&surface, &p, pu, pv)
                                .or_else(|| project_point_to_uv(&surface, &p))
                                .or_else(|| {
                                    // Endpoint-interpolated seed: the
                                    // independent Newton can diverge
                                    // along fillet-sweep trims (PK-095
                                    // — sample loss collapsed these
                                    // pcurves to straight 2-point
                                    // lines, 5e-2 off).
                                    project_point_to_uv_near(&surface, &p, seed_u, seed_v)
                                })
                        } else {
                            project_point_to_uv(&surface, &p)
                                .or_else(|| project_point_to_uv_near(&surface, &p, seed_u, seed_v))
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

    /// Recompute pcurves for ALL faces in the BRep.
    pub fn recompute_all_pcurves(&mut self) -> Result<(), crate::error::CoreError> {
        let faces: Vec<_> = self.faces.keys().collect();
        for fh in faces {
            self.recompute_face_pcurves(fh)?;
        }
        Ok(())
    }

    /// Detect and replace stale pcurves whose endpoints no longer match
    /// the actual edge curve positions or whose sampled trace no longer
    /// matches the authoritative 3D edge curve.
    ///
    /// A pcurve is stale when either endpoint fails to evaluate to the
    /// half-edge's vertices, or when the full measured edge/pcurve deviation
    /// exceeds the BRep imprecision budget. This happens when fillets or other
    /// operations move edge positions or reassign face surfaces without
    /// updating inherited pcurves.
    pub fn fix_stale_pcurves(&mut self) {
        use crate::geometry::surface::SurfaceEval;

        let he_keys: Vec<_> = self.half_edges.keys().collect();
        for heh in he_keys {
            let pcurve = match self.half_edges[heh].pcurve.as_ref() {
                Some(pc) => pc.clone(),
                None => continue,
            };

            let loop_handle = self.half_edges[heh].loop_;
            let face_handle = match self.loops.get(loop_handle) {
                Some(l) => l.face,
                None => continue,
            };
            let surface = match self.faces.get(face_handle).and_then(|f| f.surface.as_ref()) {
                Some(s) => s.geometry.clone(),
                None => continue,
            };

            let (t0, t1) = pcurve.parameter_range();
            let uv_start = pcurve.point_at(t0);
            let uv_end = pcurve.point_at(t1);
            let start_from_pcurve = surface.point_at(uv_start.x, uv_start.y);
            let end_from_pcurve = surface.point_at(uv_end.x, uv_end.y);
            let actual_start = self.vertices[self.half_edges[heh].origin].position;
            let next_heh = self.half_edges[heh].next;
            let actual_end = self.vertices[self.half_edges[next_heh].origin].position;

            let endpoint_error = start_from_pcurve
                .distance_to(&actual_start)
                .max(end_from_pcurve.distance_to(&actual_end));
            let trace_error = self
                .measure_half_edge_pcurve_deviation(heh)
                .map_or(0.0, |report| report.max_deviation);
            let stale_tolerance = self.max_imprecision.max(self.tolerance.linear());
            if endpoint_error > stale_tolerance || trace_error > stale_tolerance {
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
