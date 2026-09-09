//! Helpers for [`super::pcurve_recompute`].  Extracted to keep
//! that file under the 500-line structural cap.

use super::FaceHandle;
use super::brep::BRep;

impl BRep {
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
    pub(super) fn fix_pcurve_endpoints(&mut self, face_handle: FaceHandle) {
        for (heh, fixed) in self.compute_pcurve_endpoint_fixes(face_handle) {
            self.half_edges[heh].pcurve = Some(fixed);
        }
    }

    /// The read-only compute half of [`BRep::fix_pcurve_endpoints`]:
    /// per half-edge of `face_handle`, the endpoint-snapped
    /// replacement pcurve (when a fix is needed), in loop order. Each
    /// half-edge's fix reads only its own stored pcurve, the face
    /// surface, and vertex positions, so faces (and half-edges) are
    /// independent — callers may compute several faces in parallel
    /// and apply the writes serially (PK-302).
    pub(super) fn compute_pcurve_endpoint_fixes(
        &self,
        face_handle: FaceHandle,
    ) -> Vec<(
        crate::topology::HalfEdgeHandle,
        crate::geometry::pcurve::Pcurve,
    )> {
        use crate::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};

        let mut fixes = Vec::new();
        let surface = match self.faces.get(face_handle).and_then(|f| f.surface.as_ref()) {
            Some(s) => s.geometry.clone(),
            None => return fixes,
        };
        let (period_u, period_v) = surface.periodicity();

        let face = match self.faces.get(face_handle) {
            Some(f) => f,
            None => return fixes,
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

            let origin_pos = self.vertices[self.half_edges[heh].origin].position;
            let next_heh = self.half_edges[heh].next;
            let dest_pos = self.vertices[self.half_edges[next_heh].origin].position;

            let (t0, t1) = pcurve.parameter_range();
            let eval_start = pcurve.point_at(t0);
            let eval_end = pcurve.point_at(t1);

            // CLOSED edges: a closed ring's pcurve may legitimately
            // start anywhere on the ring — the vertex is an arbitrary
            // seam anchor, not a geometric constraint (sweep-authored
            // cylinder caps start 90° from the seam vertex and are
            // correct). Dragging a closed pcurve's end control points
            // to the vertex UV WARPS the ring (measured: a boolean-
            // transplanted cap rim became an off-circle curve whose
            // phantom second-lap crossing broke the SECOND boolean on
            // the body). If the ring closes on itself consistently,
            // there is nothing to anchor — leave it alone.
            if origin_pos.distance_to(&dest_pos) < 1e-9 {
                let uv_closure = ((eval_end.x - eval_start.x).powi(2)
                    + (eval_end.y - eval_start.y).powi(2))
                .sqrt();
                if uv_closure < tol {
                    continue;
                }
            }
            let exact_uv_start =
                match project_point_to_uv_near(&surface, &origin_pos, eval_start.x, eval_start.y)
                    .or_else(|| project_point_to_uv(&surface, &origin_pos))
                {
                    Some((u, v)) => crate::geometry::Point3::new(u, v, 0.0),
                    None => continue,
                };
            let exact_uv_end =
                match project_point_to_uv_near(&surface, &dest_pos, eval_end.x, eval_end.y)
                    .or_else(|| project_point_to_uv(&surface, &dest_pos))
                {
                    Some((u, v)) => crate::geometry::Point3::new(u, v, 0.0),
                    None => continue,
                };

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
            // Model-tolerance-scaled: a fixed 1e-9 flapped against
            // join surfaces whose radius is MEASURED from projections
            // onto fitted offsets (deviation ~1e-9..1e-8) — the bypass
            // failed by a hair at sphere poles and the endpoint
            // surgery dragged the pole control point to the undefined
            // atan2(0,0) u, tearing the tessellation chain ~3e-4 off
            // its neighbours (PK-325 wavy-plate corner caps, 23 open
            // welded mesh edges).
            let pos_tol = self.tolerance.linear().max(1e-9_f64);
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

            fixes.push((heh, fixed.into()));
        }
        fixes
    }
}

impl BRep {
    /// Verdict on an EXISTING pcurve during face recompute: `true` =
    /// keep it (skip recompute), `false` = it was cleared and must be
    /// recomputed. `POPKERN_TRACE=pcurve` prints each verdict.
    ///
    /// Cleared cases, periodic charts only:
    /// - degenerate span (~0 UV extent) under a real edge curve —
    ///   e.g. self-loop circles on a cylinder;
    /// - STALE IMAGE: a healthy SPAN is not a healthy pcurve — after
    ///   a surface swap (draft's cylinder→cone) the old chart's
    ///   pcurve keeps its span but maps off the edge curve, and
    ///   fix_pcurve_endpoints then snaps its ends to the vertex UV
    ///   and collapses the rim to zero span. The check is
    ///   min-distance from the pcurve's midpoint image to SAMPLES of
    ///   the edge curve (not the fraction-matched point: a healthy
    ///   pcurve with different parameter speed images ONTO the curve
    ///   at some other fraction).
    pub(super) fn keep_or_clear_existing_pcurve(
        &mut self,
        heh: crate::topology::HalfEdgeHandle,
        surface: &crate::geometry::surface::SurfaceGeometry,
        is_periodic: bool,
    ) -> bool {
        use crate::geometry::curve::CurveEval;
        use crate::geometry::surface::SurfaceEval;
        let trace = crate::trace::trace_enabled("pcurve");
        let Some(existing) = self.half_edges[heh].pcurve.as_ref() else {
            return false;
        };
        let (t0, t1) = existing.parameter_range();
        let s = existing.point_at(t0);
        let e = existing.point_at(t1);
        if trace {
            eprintln!(
                "  [pcurve] {heh:?} existing uv ({:.3},{:.3})->({:.3},{:.3})",
                s.x, s.y, e.x, e.y
            );
        }
        // Non-periodic self-loops use the 3D trim fallback.
        if !is_periodic {
            return true;
        }
        let uv_span = ((e.x - s.x).powi(2) + (e.y - s.y).powi(2)).sqrt();
        let edge_h = self.half_edges[heh].edge;
        let long_curve = self.edges[edge_h]
            .curve
            .as_ref()
            .is_some_and(|c| (c.t_end - c.t_start).abs() > 1e-6);
        let image_stale = long_curve && uv_span > 1e-6 && {
            let mid_uv = existing.point_at(0.5 * (t0 + t1));
            let on_surface = surface.point_at(mid_uv.x, mid_uv.y);
            self.edges[edge_h].curve.as_ref().is_some_and(|c| {
                let gate = self.tolerance.linear().max(1.0e-9) * 10.0;
                let min_d = (0..=8)
                    .map(|i| {
                        let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 8.0;
                        on_surface.distance_to(&c.geometry.point_at(t))
                    })
                    .fold(f64::INFINITY, f64::min);
                min_d > gate
            })
        };
        if (uv_span > 1e-6 && !image_stale) || !long_curve {
            if trace {
                eprintln!("  [pcurve] {heh:?} KEPT (span {uv_span:.3e} long {long_curve})");
            }
            return true;
        }
        if trace {
            eprintln!(
                "  [pcurve] {heh:?} cleared ({}), recomputing",
                if image_stale {
                    "stale image"
                } else {
                    "degenerate span"
                }
            );
        }
        self.half_edges[heh].pcurve = None;
        false
    }
}

/// One-line skip diagnostic for the per-half-edge recompute loop:
/// every silent `continue` names itself under `POPKERN_TRACE=pcurve`
/// (PK-463: a projection-failure skip served a pcurve-less face and
/// the mesh tore — the skip was invisible).
pub(super) fn trace_skip(heh: crate::topology::HalfEdgeHandle, reason: &str) {
    if crate::trace::trace_enabled("pcurve") {
        eprintln!("  [pcurve] {heh:?} SKIP: {reason}");
    }
}

impl BRep {
    /// Survival rule (PK-087) + branch consistency (PK-125) for one
    /// recomputed pcurve candidate. A recomputed pcurve only replaces
    /// an existing one when it satisfies the fraction-matched trim
    /// contract at least as well — construction passes author dense
    /// pcurves with better-than-refit fidelity (a blanket recompute
    /// degraded surgery-authored 33-point pcurves 8.5e-7 → 1.35e-3).
    /// Half-edges without an existing pcurve or without a 3-D curve
    /// to measure against are written unconditionally. An analytic
    /// exact-plane candidate replaces the projected fit when
    /// available (exactness verified, see pcurve_recompute_analytic).
    ///
    /// Branch consistency: on periodic surfaces the projection's
    /// principal branch can sit whole periods from the face's chart
    /// window; the 3-D image is identical there, so the contract
    /// comparison is decided by trig ULPs — and a wrong-branch
    /// candidate that "wins" the tie throws every chart-space
    /// consumer a full period off (measured: pyramid4_mixed support
    /// pcurve replaced by its u − 2π copy). The candidate is
    /// translated onto the EXISTING pcurve's branch first — an exact,
    /// image-preserving move.
    ///
    /// Returns `true` when the candidate was installed.
    pub(super) fn install_pcurve_with_survival(
        &mut self,
        heh: crate::topology::HalfEdgeHandle,
        surface: &crate::geometry::surface::SurfaceGeometry,
        edge_direction: (f64, f64),
        pcurve: crate::geometry::NurbsCurve,
    ) -> bool {
        use super::pcurve_recompute::pcurve_survival::pcurve_fraction_contract_error;
        let edge_curve = &self.edges[self.half_edges[heh].edge].curve;
        let mut candidate: crate::geometry::pcurve::Pcurve = edge_curve
            .as_ref()
            .and_then(|c| {
                super::pcurve_recompute_analytic::exact_plane_pcurve(
                    surface,
                    &c.geometry,
                    edge_direction.0,
                    edge_direction.1,
                )
            })
            .unwrap_or_else(|| pcurve.into());
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
                pcurve_fraction_contract_error(surface, existing, curve_on_edge, edge_direction)
                    <= pcurve_fraction_contract_error(
                        surface,
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
        !keep_existing
    }
}
