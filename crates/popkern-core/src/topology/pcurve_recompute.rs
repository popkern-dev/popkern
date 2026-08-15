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
            // Existing pcurves are kept UNLESS degenerate or
            // stale-imaged — see keep_or_clear_existing_pcurve for
            // the verdict (and PK_PCURVE_TRACE for its trace).
            if self.half_edges[heh].pcurve.is_some()
                && self.keep_or_clear_existing_pcurve(heh, &surface, is_periodic)
            {
                continue;
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
            let edge_is_straight = super::pcurve_recompute_analytic::edge_is_straight(edge_curve);

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
                // A CLOSED edge on a u-periodic surface must sweep a
                // full period in UV: both endpoints project to the
                // same u, and walks pinned to equal endpoint seeds
                // collapse the rim to a zero-span pcurve (measured:
                // drafted-cone rims at (2π,v)→(2π,v), the mesher then
                // drops the cone side of both rims — ~640 cracked
                // weld edges on an otherwise exact frustum). Unwrap
                // the destination seed by ±period, direction taken
                // from an early interior sample; the unwrapped u must
                // NOT be clamped back into the base period.
                let closed_edge = he_origin_pos.distance_to(&he_dest_pos) < 1e-9;
                let period_unwrapped =
                    closed_edge && (t_hi - t_lo).abs() > 1e-9 && surface.periodicity().0.is_some();
                let (u_e, v_e) = if period_unwrapped {
                    let period = surface.periodicity().0.unwrap_or(std::f64::consts::TAU);
                    let t_probe = t_lo + (t_hi - t_lo) / 16.0;
                    let p_probe = curve_on_edge.geometry.point_at(t_probe);
                    match project_point_to_uv(&surface, &p_probe) {
                        Some((u_probe, _)) => {
                            let mut du = u_probe - u_s;
                            while du > period / 2.0 {
                                du -= period;
                            }
                            while du < -period / 2.0 {
                                du += period;
                            }
                            (u_s + period * du.signum(), v_e)
                        }
                        None => (u_e, v_e),
                    }
                } else {
                    (u_e, v_e)
                };
                // A degree-1 pcurve is uniform-speed in UV, so the
                // shortcut is only exact when the curve's
                // PARAMETERIZATION maps uniformly onto the UV line —
                // not merely when the image is straight. A single
                // midpoint sample accepted image-straight curves with
                // non-uniform parameterization (measured: cube_extreme
                // cap seam, fraction-matched deviation 1.42e-2 while
                // the image was 3.8e-4 off); the quarter points catch
                // the skew the symmetric midpoint cannot.
                let uv_linear = [0.25, 0.5, 0.75].iter().all(|&fraction| {
                    let t = t_lo + (t_hi - t_lo) * fraction;
                    let p = curve_on_edge.geometry.point_at(t);
                    match project_point_to_uv(&surface, &p) {
                        Some((u, v)) => {
                            let expected_u = u_s + (u_e - u_s) * fraction;
                            // Compare u on the expected branch: raw
                            // projections live in the base period,
                            // while an unwrapped sweep walks out of it.
                            let u = if period_unwrapped {
                                let period =
                                    surface.periodicity().0.unwrap_or(std::f64::consts::TAU);
                                let mut u = u;
                                while u - expected_u > period / 2.0 {
                                    u -= period;
                                }
                                while u - expected_u < -period / 2.0 {
                                    u += period;
                                }
                                u
                            } else {
                                u
                            };
                            let expected_v = v_s + (v_e - v_s) * fraction;
                            let deviation =
                                ((u - expected_u).powi(2) + (v - expected_v).powi(2)).sqrt();
                            deviation < 1e-4
                        }
                        None => false,
                    }
                });

                if uv_linear {
                    // UV path is linear -- exact degree-1 pcurve.
                    // Periodic-unwrapped u endpoints stay unclamped —
                    // clamping would fold the full-period sweep back
                    // onto a zero-span line.
                    let clamp_u = |val: f64| {
                        if period_unwrapped {
                            val
                        } else {
                            clamp_if_finite(val, u_min, u_max)
                        }
                    };
                    NurbsCurve {
                        degree: 1,
                        knots: vec![0.0, 0.0, 1.0, 1.0],
                        control_points: vec![
                            Point3::new(clamp_u(u_s), clamp_if_finite(v_s, v_min, v_max), 0.0),
                            Point3::new(clamp_u(u_e), clamp_if_finite(v_e, v_min, v_max), 0.0),
                        ],
                        weights: vec![1.0, 1.0],
                    }
                } else {
                    // UV path is curved — dense sampling.
                    // OCCT-style pcurve computation (per ShapeConstruct_
                    // ProjectCurveOnSurface): sample the 3D curve,
                    // project by CONTINUATION (endpoint- and midpoint-
                    // anchored walks race and the smaller measured
                    // deviation wins — see pcurve_projection_walk),
                    // unwrap seam crossings, then INTERPOLATE (not fit).
                    let uv_points = super::pcurve_projection_walk::anchored_uv_samples(
                        &surface,
                        curve_on_edge,
                        (t_lo, t_hi),
                        he_origin_pos,
                        he_dest_pos,
                        ((u_s, v_s), (u_e, v_e)),
                        23, // OCCT default sample count
                    );

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

            // Survival rule: a recomputed pcurve only replaces an
            // existing one when it satisfies the fraction-matched trim
            // contract at least as well. Construction passes author
            // dense pcurves with better-than-refit fidelity (PK-087:
            // a blanket recompute degraded surgery-authored 33-point
            // pcurves 8.5e-7 → 1.35e-3, and a refit over a partially
            // diverged projection walk tore pentagon5_blunt fan seam
            // normals to 2.9e-1). Half-edges without an existing
            // pcurve or without a 3-D curve to measure against are
            // written unconditionally, as before.
            // The survival rule below still measures an analytic
            // candidate against any existing pcurve — exactness is
            // verified, not assumed (see pcurve_recompute_analytic).
            let mut candidate: crate::geometry::pcurve::Pcurve = edge_curve
                .as_ref()
                .and_then(|c| {
                    super::pcurve_recompute_analytic::exact_plane_arc_pcurve(
                        &surface,
                        &c.geometry,
                        edge_direction.0,
                        edge_direction.1,
                    )
                })
                .unwrap_or_else(|| pcurve.into());
            // Branch consistency (PK-125): on periodic surfaces the
            // projection's principal branch can sit whole periods away
            // from the face's chart window. The 3-D image is identical
            // there, so the contract comparison below is decided by
            // trig ULPs — and a wrong-branch candidate that "wins" the
            // tie throws every chart-space consumer (boundary-row
            // relocks, trim fits, seam measures) a full period off
            // (measured: pyramid4_mixed support pcurve replaced by its
            // u − 2π copy; the export-repair relock chased the
            // off-window trace into a cap-surface fold). Translate the
            // candidate onto the EXISTING pcurve's branch first — an
            // exact, image-preserving move.
            if let Some(existing) = self.half_edges[heh].pcurve.as_ref() {
                let (period_u, period_v) = surface.periodicity();
                let branch_shift = |delta: f64, period: Option<f64>| -> f64 {
                    match period {
                        Some(period) if period > 0.0 => (delta / period).round() * period,
                        _ => 0.0,
                    }
                };
                let midpoint_of = |pc: &crate::geometry::pcurve::Pcurve| {
                    let (t0, t1) = pc.parameter_range();
                    pc.point_at(0.5 * (t0 + t1))
                };
                let target = midpoint_of(existing);
                let current = midpoint_of(&candidate);
                let shift_u = branch_shift(target.x - current.x, period_u);
                let shift_v = branch_shift(target.y - current.y, period_v);
                if shift_u != 0.0 || shift_v != 0.0 {
                    candidate = candidate.translated_uv(shift_u, shift_v);
                }
            }
            let keep_existing = match (&self.half_edges[heh].pcurve, edge_curve) {
                (Some(existing), Some(curve_on_edge)) => {
                    pcurve_fraction_contract_error(
                        &surface,
                        existing,
                        curve_on_edge,
                        edge_direction,
                    ) <= pcurve_fraction_contract_error(
                        &surface,
                        &candidate,
                        curve_on_edge,
                        edge_direction,
                    )
                }
                _ => false,
            };
            if !keep_existing {
                self.half_edges[heh].pcurve = Some(candidate);
            }
        }

        self.fix_pcurve_endpoints(face_handle);

        Ok(())
    }

    /// Recompute pcurves for ALL faces in the BRep.
    pub fn recompute_all_pcurves(&mut self) -> Result<(), crate::error::CoreError> {
        let faces: Vec<_> = self.faces.keys().collect();
        for fh in faces {
            self.recompute_face_pcurves(fh)?;
            self.close_periodic_loop_winding(fh);
        }
        Ok(())
    }

    /// [`Self::recompute_all_pcurves`], fault-tolerant per face: one
    /// face's projection failure does not abandon the faces after it
    /// — a dress-up surgery that nulled stale pcurves (fix_stale)
    /// relies on this pass to refill them, and an early abort left
    /// the remaining faces meshing against their pre-surgery extents
    /// (measured: a chamfered arm tip's rest face triangulated clear
    /// across the setbacks). Kept SEPARATE from the strict variant:
    /// the direct-apex pipeline is tuned to the abort semantics
    /// (measured: the lenient loop broke eight of its cap gates).
    pub fn recompute_all_pcurves_lenient(&mut self) {
        let faces: Vec<_> = self.faces.keys().collect();
        for fh in faces {
            let _ = self.recompute_face_pcurves(fh);
            self.close_periodic_loop_winding(fh);
        }
    }

    /// Zero the net U-winding of every loop of a face on a
    /// U-periodic chart.
    ///
    /// The closed-edge pcurve direction convention (canonical
    /// half-edge forward, twin backward) makes the two uses of ONE
    /// edge consistent, but says nothing about a LOOP that chains
    /// several closed rims: a recess wall bounded by its bottom
    /// circle and the mouth ellipse got both rims wound the same
    /// way, the ring's net winding came to ±2 periods, and the
    /// boundary-first mesher rightly declined to chart it. The
    /// chart-valid invariant is per loop: after branch-chaining each
    /// pcurve to its predecessor, the net winding must be zero.
    /// Full-period members are flipped (via [`Pcurve::reversed`])
    /// until it is.
    fn close_periodic_loop_winding(&mut self, face_handle: FaceHandle) {
        let Some(u_period) = self
            .faces
            .get(face_handle)
            .and_then(|f| f.surface.as_ref())
            .and_then(|s| s.geometry.periodicity().0)
        else {
            return;
        };

        let Some(face) = self.faces.get(face_handle) else {
            return;
        };
        let mut loop_handles = vec![face.outer_loop];
        loop_handles.extend(face.inner_loops.iter().copied());

        for lh in loop_handles {
            let Some(loop_data) = self.loops.get(lh) else {
                continue;
            };
            let start = loop_data.half_edge;
            if !self.half_edges.contains_key(start) {
                continue;
            }
            let mut hes = Vec::new();
            let mut cursor = start;
            let limit = self.half_edges.len() + 1;
            for _ in 0..limit {
                hes.push(cursor);
                cursor = self.half_edges[cursor].next;
                if cursor == start {
                    break;
                }
            }
            if hes.len() < 2 {
                continue;
            }

            // Iterate: chain branches, measure net winding, flip one
            // reversible full-period member, repeat. Bounded by the
            // loop length.
            for _ in 0..hes.len() {
                let mut prev_end: Option<f64> = None;
                let mut first_start: Option<f64> = None;
                let mut deltas: Vec<f64> = Vec::with_capacity(hes.len());
                let mut complete = true;
                for &heh in &hes {
                    let Some(pc) = self.half_edges[heh].pcurve.as_ref() else {
                        complete = false;
                        break;
                    };
                    let (t0, t1) = pc.parameter_range();
                    let mut su = pc.point_at(t0).x;
                    let mut eu = pc.point_at(t1).x;
                    if let Some(pe) = prev_end {
                        let k = ((pe - su) / u_period).round();
                        if k != 0.0 {
                            if let Some(pc_mut) = self.half_edges[heh].pcurve.as_mut() {
                                pc_mut.shift_u(k * u_period);
                            }
                            su += k * u_period;
                            eu += k * u_period;
                        }
                    } else {
                        first_start = Some(su);
                    }
                    deltas.push(eu - su);
                    prev_end = Some(eu);
                }
                if !complete {
                    break;
                }
                let (Some(fs), Some(pe)) = (first_start, prev_end) else {
                    break;
                };
                let net = pe - fs;
                if net.abs() <= 0.5 * u_period {
                    break;
                }
                let mut flipped = false;
                for (i, &heh) in hes.iter().enumerate() {
                    let d = deltas[i];
                    if d.abs() >= u_period * (1.0 - 1e-6)
                        && (net - 2.0 * d).abs() < net.abs() - 0.25 * u_period
                        && let Some(rev) =
                            self.half_edges[heh].pcurve.as_ref().map(|p| p.reversed())
                    {
                        self.half_edges[heh].pcurve = Some(rev);
                        flipped = true;
                        break;
                    }
                }
                if !flipped {
                    break;
                }
            }
            // PRINCIPAL-CELL NORMALIZATION: branch-chaining anchors
            // the loop to whatever cell its first pcurve landed on,
            // and analytic inversion of a rim whose arc starts past
            // +π puts the whole loop on the 2π..4π cell. That is
            // chart-equivalent for the kernel but OUTSIDE the
            // surface's principal parameter cell for interchange
            // (measured: OCC's BRepCheck rejected a twice-chamfered
            // column's bottom cone strip while the mesh welded
            // clean). Shift the WHOLE loop by whole periods so its
            // mean u lands in [0, period) — a rigid translation that
            // preserves every relative relationship.
            {
                let mut min_u = f64::MAX;
                let mut max_u = f64::MIN;
                for &heh in &hes {
                    if let Some(pc) = self.half_edges[heh].pcurve.as_ref() {
                        let (t0, t1) = pc.parameter_range();
                        for t in [t0, t1] {
                            let u = pc.point_at(t).x;
                            min_u = min_u.min(u);
                            max_u = max_u.max(u);
                        }
                    }
                }
                // Scoped to loops WHOLLY outside [0, period): loops
                // that merely straddle a cell edge (charts anchored
                // at −π, seam-adjacent trims) are tuned behavior and
                // must not move (measured: shifting by the mean
                // broke eight direct-apex cap gates).
                if min_u <= max_u {
                    let eps = 1e-9 * u_period;
                    let k = if min_u >= u_period - eps {
                        (min_u / u_period + eps).floor()
                    } else if max_u <= eps {
                        (max_u / u_period - eps).ceil() - 1.0
                    } else {
                        0.0
                    };
                    if k != 0.0 {
                        for &heh in &hes {
                            if let Some(pc) = self.half_edges[heh].pcurve.as_mut() {
                                pc.shift_u(-k * u_period);
                            }
                        }
                    }
                }
            }
        }
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
    /// Re-window edge 3D curves whose parameter window no longer
    /// matches the edge's vertices. A topology surgery that splices
    /// new vertices into existing edges (a chamfer shortening the
    /// edges adjacent to its strip) can leave the CURVE window
    /// covering the old span; every pcurve rebuilt from that curve
    /// then faithfully reproduces the stale extent (measured: a
    /// chamfered slotted plate's wall pcurves sat 0.7 off their own
    /// vertices after recompute, because recompute projects curve
    /// samples). The repair inverts both vertex positions onto the
    /// stored curve and narrows the window to them; edges whose
    /// vertices don't lie ON the curve (within the stale tolerance)
    /// are left alone.
    pub fn fix_stale_edge_windows(&mut self) {
        use crate::geometry::curve::CurveEval;
        let stale_tolerance = self.max_imprecision.max(self.tolerance.linear());
        let edge_keys: Vec<_> = self.edges.keys().collect();
        for eh in edge_keys {
            let Some(edge) = self.edges.get(eh) else {
                continue;
            };
            let Some(curve) = edge.curve.clone() else {
                continue;
            };
            let he = edge.half_edge;
            let Some(he_data) = self.half_edges.get(he) else {
                continue;
            };
            let v0h = he_data.origin;
            let Some(twin) = self.half_edges.get(he_data.twin) else {
                continue;
            };
            let v1h = twin.origin;
            if v0h == v1h {
                continue; // closed edge — windows are full-period
            }
            let v0 = self.vertices[v0h].position;
            let v1 = self.vertices[v1h].position;
            let p0 = curve.geometry.point_at(curve.t_start);
            let p1 = curve.geometry.point_at(curve.t_end);
            let err = p0.distance_to(&v0).max(p1.distance_to(&v1));
            if err <= stale_tolerance {
                continue;
            }
            // Nearest parameter on the stored curve to a point:
            // coarse scan + ternary refinement.
            let nearest = |p: &crate::geometry::Point3| -> (f64, f64) {
                let n = 128;
                let mut best = (curve.t_start, f64::INFINITY);
                for i in 0..=n {
                    let t = curve.t_start + (curve.t_end - curve.t_start) * i as f64 / n as f64;
                    let d = curve.geometry.point_at(t).distance_to(p);
                    if d < best.1 {
                        best = (t, d);
                    }
                }
                let step = (curve.t_end - curve.t_start).abs() / n as f64;
                let (mut lo, mut hi) = (best.0 - step, best.0 + step);
                for _ in 0..48 {
                    let m1 = lo + (hi - lo) / 3.0;
                    let m2 = hi - (hi - lo) / 3.0;
                    if curve.geometry.point_at(m1).distance_to(p)
                        <= curve.geometry.point_at(m2).distance_to(p)
                    {
                        hi = m2;
                    } else {
                        lo = m1;
                    }
                }
                let t = 0.5 * (lo + hi);
                (t, curve.geometry.point_at(t).distance_to(p))
            };
            let (t0, d0) = nearest(&v0);
            let (t1, d1) = nearest(&v1);
            if d0 <= stale_tolerance
                && d1 <= stale_tolerance
                && (t1 - t0).abs() > 1e-12
                && let Some(edge_mut) = self.edges.get_mut(eh)
                && let Some(c) = edge_mut.curve.as_mut()
            {
                c.t_start = t0;
                c.t_end = t1;
            }
        }
    }
}

#[path = "pcurve_stale_repair.rs"]
mod pcurve_stale_repair;

#[path = "pcurve_region_recompute.rs"]
mod pcurve_region_recompute;
#[path = "pcurve_survival.rs"]
mod pcurve_survival;
use pcurve_survival::pcurve_fraction_contract_error;

#[cfg(test)]
#[path = "pcurve_recompute_tests.rs"]
mod tests;
