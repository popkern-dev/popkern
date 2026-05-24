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
        use crate::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};

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

            let origin_pos = self.vertices[self.half_edges[heh].origin].position;
            let next_heh = self.half_edges[heh].next;
            let dest_pos = self.vertices[self.half_edges[next_heh].origin].position;

            let (t0, t1) = pcurve.parameter_range();
            let eval_start = pcurve.point_at(t0);
            let eval_end = pcurve.point_at(t1);
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
}
