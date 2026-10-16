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
        self.recompute_face_pcurves_impl(face_handle, false)
    }

    /// [`Self::recompute_face_pcurves`] with the loop-winding repair
    /// postcondition armed (see `pcurve_loop_winding_repair`). Reserved
    /// for from-scratch re-authoring of bodies whose orientation flags
    /// are trusted (shell/offset assembly output) — established
    /// pipelines (boolean repair tails) carry closed-rim direction
    /// conventions their consumers rely on, and the repair must not
    /// re-orient those.
    pub fn recompute_face_pcurves_winding_repaired(
        &mut self,
        face_handle: FaceHandle,
    ) -> Result<(), crate::error::CoreError> {
        self.recompute_face_pcurves_impl(face_handle, true)
    }

    fn recompute_face_pcurves_impl(
        &mut self,
        face_handle: FaceHandle,
        repair_loop_winding: bool,
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

        // Half-edges whose pcurves THIS call writes — the loop-winding
        // repair below only re-orients loops re-authored in full.
        let mut freshly_authored: std::collections::HashSet<super::HalfEdgeHandle> =
            std::collections::HashSet::new();
        for heh in he_handles {
            // Existing pcurves are kept UNLESS degenerate or
            // stale-imaged — see keep_or_clear_existing_pcurve for
            // the verdict (and the pcurve trace topic for its trace).
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
                    None => {
                        trace_skip(heh, "origin projection failed");
                        continue;
                    }
                };
                let (u_e, v_e) = match project_point_to_uv(&surface, &he_dest_pos) {
                    Some(uv) => uv,
                    None => {
                        trace_skip(heh, "dest projection failed");
                        continue;
                    }
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

                // Quick linearity check (quarter points, then the 3D trace):
                // a straight UV image takes an exact degree-1 pcurve.
                let (mut u_s, v_s) = match project_point_to_uv(&surface, &he_origin_pos) {
                    Some(uv) => uv,
                    None => {
                        trace_skip(heh, "origin projection failed");
                        continue;
                    }
                };
                let (raw_u_e, raw_v_e) = match project_point_to_uv(&surface, &he_dest_pos) {
                    Some(uv) => uv,
                    None => {
                        trace_skip(heh, "dest projection failed");
                        continue;
                    }
                };
                // Branch-honest sampling on periodic charts (PK-286): cold
                // projections pick a seam side by ULP and spuriously failed
                // linearity (a 7e-15 upstream change flipped a torus seam
                // pcurve's degree on double_boss_plate), so samples chain to
                // the branch nearest their predecessor in every periodic
                // direction. The chain also unwraps CLOSED edges (drafted-
                // cone rims collapsed to zero-span pcurves, ~640 cracked weld
                // edges): the chained destination lands one period on.
                let (period_u, period_v) = surface.periodicity();
                let align = |value: f64, reference: f64, period: Option<f64>| -> f64 {
                    match period {
                        Some(period) if period > 0.0 => {
                            value - ((value - reference) / period).round() * period
                        }
                        _ => value,
                    }
                };
                let mut chain = (u_s, v_s);
                let mut quarter_samples = Vec::with_capacity(3);
                let mut all_projected = true;
                for fraction in [0.25, 0.5, 0.75] {
                    let t = t_lo + (t_hi - t_lo) * fraction;
                    let p = curve_on_edge.geometry.point_at(t);
                    match project_point_to_uv(&surface, &p) {
                        Some((u, v)) => {
                            chain = (align(u, chain.0, period_u), align(v, chain.1, period_v));
                            quarter_samples.push((fraction, chain.0, chain.1));
                        }
                        None => {
                            all_projected = false;
                            break;
                        }
                    }
                }
                let mut u_e = align(raw_u_e, chain.0, period_u);
                let v_e = align(raw_v_e, chain.1, period_v);
                // Degenerate-row endpoints (chart poles: sphere poles,
                // cone apexes) project with an ARBITRARY u — the whole
                // row maps to one 3D point, so ∂u vanishes there and
                // the bogus u defeats the straight-in-UV detection:
                // meridian edges (the PK-325 corner-cap arcs) fell to
                // the projection walk, whose pole wobble (measured
                // 1.13e-4) tore tessellation chains ~3e-4 off their
                // neighbours (23 open welded mesh edges at the wavy
                // plate's four top corners). Snap the free parameter
                // to the interior chain's u — and EXEMPT snapped
                // endpoints from principal-range clamping, which is
                // what re-broke the first fix (the meridian branch
                // u = −π/2 clamped to 0).
                let mut u_snapped = false;
                if let Some(&(_, first_u, _)) = quarter_samples.first() {
                    use crate::geometry::surface::SurfaceEval;
                    // Row-degeneracy gate: ∂u under 1e-4 of the edge's
                    // interior ∂u (a projected pole lands at ∂u ≈ r·1e-6,
                    // so an exact-zero gate never fires). Not ∂v: ∂u/∂v
                    // has length units on a cylinder, so r < 1e-4 read
                    // as a pole and lost its exact line (PK-860).
                    let du_in = quarter_samples
                        .iter()
                        .map(|&(_, u, v)| surface.partial_u(u, v).length())
                        .fold(0.0_f64, f64::max);
                    let pole_at = |u: f64, v: f64| -> bool {
                        du_in > 0.0 && surface.partial_u(u, v).length() <= 1.0e-4 * du_in
                    };
                    // The free azimuth at a pole is the curve's limit:
                    // the direction it arrives from (leaves toward),
                    // read a hair inside the edge. A quarter sample's
                    // azimuth is not it when the curve swings round
                    // the pole (measured: a dome rim through the pole
                    // snapped its arriving end π off — 0.41 out of r
                    // 10; PK-646).
                    // Analytic poles only (sphere, cone apex), whose
                    // azimuth is well conditioned a hair from the pole;
                    // a NURBS chart's collapsed row is not (measured:
                    // the pyramid3 apex cap's projections there lost
                    // the cap's G1 rim, 0.975 vs a 5e-2 gate).
                    let analytic_pole = matches!(
                        surface,
                        crate::geometry::surface::SurfaceGeometry::Sphere { .. }
                            | crate::geometry::surface::SurfaceGeometry::Cone { .. }
                    );
                    let limit_u = |fraction: f64, reference: f64| -> f64 {
                        let t = t_lo + (t_hi - t_lo) * fraction;
                        match analytic_pole {
                            true => {
                                project_point_to_uv(&surface, &curve_on_edge.geometry.point_at(t))
                                    .map_or(reference, |(u, _)| align(u, reference, period_u))
                            }
                            false => reference,
                        }
                    };
                    if pole_at(u_s, v_s) {
                        u_s = limit_u(1.0e-7, first_u);
                        u_snapped = true;
                    }
                    if pole_at(u_e, v_e) {
                        u_e = limit_u(1.0 - 1.0e-7, chain.0);
                        u_snapped = true;
                    }
                }
                // A degree-1 pcurve is uniform-speed in UV: exact only when
                // the PARAMETERIZATION maps uniformly onto the line (a
                // midpoint test passed cube_extreme's skewed cap seam,
                // 1.42e-2 off; quarter points catch it), and only when
                // the line traces the edge everywhere (PK-863).
                let gate = 10.0
                    * self.edges[edge_handle]
                        .tolerance
                        .unwrap_or(self.tolerance.linear());
                let line_holds = |a, b| {
                    super::pcurve_recompute_analytic::uv_line_traces_edge(
                        &surface,
                        curve_on_edge,
                        (t_lo, t_hi),
                        a,
                        b,
                        gate,
                    )
                };
                let uv_linear = all_projected
                    && quarter_samples.iter().all(|&(fraction, u, v)| {
                        let expected_u = u_s + (u_e - u_s) * fraction;
                        let expected_v = v_s + (v_e - v_s) * fraction;
                        let deviation =
                            ((u - expected_u).powi(2) + (v - expected_v).powi(2)).sqrt();
                        deviation < 1e-4
                    })
                    && line_holds((u_s, v_s), (u_e, v_e));

                if uv_linear {
                    // UV path is linear -- exact degree-1 pcurve.
                    // Endpoints stay unclamped in a periodic direction
                    // whose chain left the principal branch — clamping
                    // would fold an unwrapped sweep back onto a
                    // zero-span line. (`align` shifts by exact period
                    // multiples, so an unshifted value compares
                    // bit-equal to its raw projection.)
                    let u_unwrapped = u_e != raw_u_e || u_snapped;
                    let v_unwrapped = v_e != raw_v_e;
                    let clamp_u = |val: f64| {
                        if u_unwrapped {
                            val
                        } else {
                            clamp_if_finite(val, u_min, u_max)
                        }
                    };
                    let clamp_v = |val: f64| {
                        if v_unwrapped {
                            val
                        } else {
                            clamp_if_finite(val, v_min, v_max)
                        }
                    };
                    NurbsCurve {
                        degree: 1,
                        knots: vec![0.0, 0.0, 1.0, 1.0],
                        control_points: vec![
                            Point3::new(clamp_u(u_s), clamp_v(v_s), 0.0),
                            Point3::new(clamp_u(u_e), clamp_v(v_e), 0.0),
                        ],
                        weights: vec![1.0, 1.0],
                    }
                } else if let Some(p) = super::pcurve_recompute_analytic::exact_plane_placeholder(
                    &surface,
                    curve_on_edge,
                    (t_lo, t_hi),
                    [(u_s, v_s), (u_e, v_e)],
                ) {
                    p
                } else {
                    // Curved UV path, OCCT-style: project by continuation
                    // (anchored walks race: pcurve_projection_walk), unwrap
                    // seams, INTERPOLATE piecewise between the curve's
                    // breakpoints, refine until the image follows (PK-555).
                    let sample = |fractions: &[f64]| {
                        super::pcurve_projection_walk::anchored_uv_samples(
                            &surface,
                            curve_on_edge,
                            (t_lo, t_hi),
                            he_origin_pos,
                            he_dest_pos,
                            ((u_s, v_s), (u_e, v_e)),
                            fractions,
                        )
                    };
                    let fitted = super::pcurve_projection_fit::measured_projection_interpolant(
                        &surface,
                        curve_on_edge,
                        (t_lo, t_hi),
                        self.tolerance.linear(),
                        &sample,
                    );
                    if let Some(fitted) = fitted {
                        fitted
                    } else if let uv_points =
                        sample(&super::pcurve_projection_fit::segment_fractions(
                            &[],
                            &[super::pcurve_projection_fit::PROJECTION_SAMPLES - 1],
                        ))
                        && uv_points.len() >= 2
                    {
                        NurbsCurve {
                            degree: 1,
                            knots: vec![0.0, 0.0, 1.0, 1.0],
                            control_points: vec![uv_points[0], uv_points[uv_points.len() - 1]],
                            weights: vec![1.0, 1.0],
                        }
                    } else {
                        trace_skip(heh, "too few uv samples");
                        continue;
                    }
                }
            };

            // Survival rule + branch consistency — extracted to
            // `install_pcurve_with_survival` (pcurve_recompute_helpers)
            // under the file-size ratchet.
            if self.install_pcurve_with_survival(heh, &surface, edge_direction, pcurve) {
                freshly_authored.insert(heh);
            }
        }

        self.fix_pcurve_endpoints(face_handle);

        // Per-loop chart coherence is part of THIS operation's
        // postcondition (PK-286), not a batch-wrapper courtesy: the
        // independent projections above land each pcurve in the
        // principal branch, so a seam-adjacent recompute could leave
        // one loop member a whole period from its mates. Direct
        // per-face callers (seam welds, cap authors, repair tails)
        // hit exactly that tear — measured: the PK-159 parallel-seam
        // completion meshed 12 cracked weld edges when a seam-vertex
        // projection anchored one rim arc a period below its loop.
        self.close_periodic_loop_winding(face_handle);

        // Loop-winding honesty, when armed: a closed-edge member's
        // direction comes from the per-edge canonical bit, which cannot
        // express the direction each face's CHART needs. Fully
        // re-authored outer loops that contradict their face's
        // orientation flag get their closed members flipped
        // transactionally; partially re-authored loops keep their
        // trusted directions (see pcurve_loop_winding_repair).
        if repair_loop_winding {
            self.fix_loop_pcurve_winding(face_handle, &freshly_authored);
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

    /// [`Self::recompute_all_pcurves`] with the loop-winding repair
    /// armed on every face — for from-scratch re-authoring of bodies
    /// with trusted orientation flags (shell/offset assembly output,
    /// which copies topology without pcurves).
    pub fn recompute_all_pcurves_winding_repaired(
        &mut self,
    ) -> Result<(), crate::error::CoreError> {
        let faces: Vec<_> = self.faces.keys().collect();
        for fh in faces {
            self.recompute_face_pcurves_winding_repaired(fh)?;
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

#[path = "pcurve_loop_chart_coherence.rs"]
mod pcurve_loop_chart_coherence;

#[path = "pcurve_stale_repair.rs"]
mod pcurve_stale_repair;

#[path = "pcurve_region_recompute.rs"]
mod pcurve_region_recompute;
#[path = "pcurve_survival.rs"]
pub(super) mod pcurve_survival;
use super::pcurve_recompute_helpers::trace_skip;

#[cfg(test)]
#[path = "pcurve_recompute_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "pcurve_winding_reauthor_tests.rs"]
mod winding_reauthor_tests;
