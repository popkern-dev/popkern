//! The pair-collecting form of the self-intersection check — the
//! detection authority the validator and the global-SI trim stage
//! share (PK-327). Split from `self_intersection.rs` (PK-329) to
//! honour the 500-line structural cap.

use crate::geometry::intersection::degeneracy::DefaultDegeneracyHandler;
use crate::geometry::intersection::intersect_surfaces_with_handler;
use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::surface::SurfaceEval;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle};
use crate::validation::bvh::Bvh;
use std::collections::HashSet;

use super::self_intersection_helpers::{TraceSurfacePatch, trace_self_intersection_candidate};
use super::self_intersection_same_surface::same_surface_domains_disjoint;

/// A face's pcurve-sampled `(u, v)` trim window.
type UvWindow = ((f64, f64), (f64, f64));

impl BRep {
    /// The pair-collecting form of the self-intersection check — the
    /// SAME adjacency exemption, AABB/BVH pruning, trim-overlap gates,
    /// and pair budget, returning every offending pair instead of the
    /// first (`collect_all = false` stops at the first hit, which is
    /// what the validator needs). This is the detection authority the
    /// global-SI trim stage (PK-327) builds on: trimming exactly what
    /// the validator would reject keeps the serve gate coherent.
    /// `tangential_band_is_clean`: trimming context — one-sided
    /// tolerance-band contact (a join strip riding its neighbour's
    /// skin) has no material crossing and nothing to trim, so the
    /// caller may skip the marching escalation the VALIDATOR keeps
    /// for its coincidence question (measured 36 s per grazing
    /// cylinder×NURBS pair).
    pub fn collect_self_intersecting_pairs_among(
        &self,
        face_keys: &[FaceHandle],
        collect_all: bool,
        tangential_band_is_clean: bool,
    ) -> Vec<(FaceHandle, FaceHandle)> {
        let mut found: Vec<(FaceHandle, FaceHandle)> = Vec::new();
        let face_keys: Vec<FaceHandle> = {
            let mut seen = std::collections::BTreeSet::new();
            face_keys
                .iter()
                .copied()
                .filter(|fh| self.faces.contains_key(*fh) && seen.insert(*fh))
                .collect()
        };
        let n = face_keys.len();
        if n < 2 {
            return found;
        }

        // ------------------------------------------------------------------ //
        // 1. Build adjacency set (by face index pairs, not handles)           //
        // ------------------------------------------------------------------ //
        let mut adjacent: HashSet<(usize, usize)> = HashSet::new();

        for (i, &fk) in face_keys.iter().enumerate() {
            let face = &self.faces[fk];
            // ALL loops: faces welded along an INNER rim are
            // adjacent exactly like outer-loop neighbours; the
            // outer-only walk missed them and ran the SI probe
            // against the pair's own shared weld.
            let loop_handles: Vec<crate::topology::LoopHandle> = std::iter::once(face.outer_loop)
                .chain(face.inner_loops.iter().copied())
                .collect();
            for lh in loop_handles {
                let Some(loop_data) = self.loops.get(lh) else {
                    continue;
                };
                let start_he_handle = loop_data.half_edge;
                // Empty loop — skip
                if !self.half_edges.contains_key(start_he_handle) {
                    continue;
                }
                let mut he_handle = start_he_handle;
                // Walk the outer loop; guard against cycles larger than the mesh.
                let bound = self.half_edges.len().max(1);
                for _ in 0..=bound {
                    let he = match self.half_edges.get(he_handle) {
                        Some(h) => h,
                        None => break,
                    };
                    let twin_handle = he.twin;
                    if let Some(twin) = self.half_edges.get(twin_handle) {
                        let twin_loop_handle = twin.loop_;
                        if let Some(twin_loop) = self.loops.get(twin_loop_handle) {
                            let twin_face = twin_loop.face;
                            if let Some(j) = face_keys.iter().position(|&k| k == twin_face)
                                && i != j
                            {
                                let pair = if i < j { (i, j) } else { (j, i) };
                                adjacent.insert(pair);
                            }
                        }
                    }
                    he_handle = he.next;
                    if he_handle == start_he_handle {
                        break;
                    }
                }
            }
        }

        // ------------------------------------------------------------------ //
        // 2. Compute face AABBs via 5×5 sampling; build BVH                  //
        // ------------------------------------------------------------------ //
        let sample_params = [0.0, 0.25, 0.5, 0.75, 1.0];

        // Per-face windows and boundary positions are cached here and
        // reused by every candidate pair below (PK-329: the pair loop
        // recomputed the pcurve-sampled window and the boundary walk
        // per PAIR — tens of ms of pure recomputation on corpus bodies
        // whose final pair list is empty). The pcurve-sampled window
        // build is itself the collection's dominant fixed cost (~1 ms
        // per face), and each face's window is a pure read of `self`,
        // so it dispatches in parallel with an ordered collect.
        let window_of = |&fk: &FaceHandle| -> Option<UvWindow> {
            self.faces[fk].surface.as_ref().map(|s| {
                // Pcurve-sampled trim window (face_uv_window): the
                // positional window collapses on seam-anchored
                // closed faces and both under-boxed the BVH and
                // tripped the degenerate-range skip below.
                crate::validation::face_uv_window::face_trim_uv_window(self, fk, &s.geometry)
            })
        };
        #[cfg(feature = "parallel")]
        let windows: Vec<Option<UvWindow>> = {
            use rayon::prelude::*;
            face_keys.par_iter().map(window_of).collect()
        };
        #[cfg(not(feature = "parallel"))]
        let windows: Vec<Option<UvWindow>> = face_keys.iter().map(window_of).collect();
        let boundary_positions: Vec<Vec<crate::geometry::Point3>> = face_keys
            .iter()
            .map(|&fk| self.face_outer_loop_positions(fk))
            .collect();

        let aabbs: Vec<Option<Aabb>> = face_keys
            .iter()
            .enumerate()
            .map(|(idx, &fk)| {
                let face = &self.faces[fk];
                face.surface.as_ref().and_then(|s| {
                    let ((u0, u1), (v0, v1)) = windows[idx]?;
                    let mut pts = Vec::with_capacity(25);
                    for &tp in &sample_params {
                        let u = u0 + tp * (u1 - u0);
                        for &tq in &sample_params {
                            let v = v0 + tq * (v1 - v0);
                            pts.push(s.geometry.point_at(u, v));
                        }
                    }
                    Aabb::from_points(&pts)
                })
            })
            .collect();

        // Build BVH items: only include faces that have valid AABBs.
        let bvh_items: Vec<(usize, Aabb)> = aabbs
            .iter()
            .enumerate()
            .filter_map(|(i, aabb)| aabb.as_ref().map(|a| (i, a.clone())))
            .collect();

        if bvh_items.len() < 2 {
            return found;
        }

        let bvh = Bvh::build(&bvh_items);

        // Use 2× linear tolerance for intersection detection: tight enough to catch
        // genuine self-intersections, loose enough for fast intersection testing.
        // This is conservative: may report near-misses as intersections (false positives)
        // but won't miss actual intersections within tolerance (no false negatives).
        let linear_tol = self.tolerance.linear();
        let coarse_tol = Tolerance::default()
            .with_linear(linear_tol * 2.0)
            .with_angular(self.tolerance.angular())
            .with_parametric(self.tolerance.parametric_base());

        // ------------------------------------------------------------------ //
        // 3. Test non-adjacent, AABB-overlapping pairs with intersect_surfaces //
        // ------------------------------------------------------------------ //
        // Validation budget: cap the number of pairs tested to prevent
        // validation from being more expensive than the operation that
        // produced the BRep. 200 pairs covers typical solids; beyond that,
        // diminishing returns on defect detection.
        const MAX_PAIRS_CHECKED: usize = 200;
        let overlapping = bvh.find_overlapping_pairs(linear_tol);
        let mut pairs_checked = 0usize;
        // Gate phase (sequential): cheap exemptions and certificates, in
        // BVH pair order, feeding the budget-capped job list the
        // parallel phase below consumes.
        let mut jobs: Vec<(usize, usize)> = Vec::new();

        for (i, j) in overlapping {
            if pairs_checked >= MAX_PAIRS_CHECKED {
                break;
            }
            let pair = if i < j { (i, j) } else { (j, i) };
            if adjacent.contains(&pair) {
                continue;
            }

            // Canonical pair orientation: the intersection solver is
            // measurably order-sensitive for NURBS×analytic pairs
            // (PK-251: cylinder-first found the strip×cap crossing,
            // NURBS-first returned nothing), so the verdict must not
            // depend on face enumeration order. Analytic surfaces go
            // first. Containment for the solver asymmetry itself —
            // see the free-surface intersection follow-up task.
            let analytic_rank = |face: crate::topology::FaceHandle| match self.faces[face]
                .surface
                .as_ref()
                .map(|s| &s.geometry)
            {
                Some(
                    crate::geometry::surface::SurfaceGeometry::Nurbs(_)
                    | crate::geometry::surface::SurfaceGeometry::Gregory(_),
                ) => 1,
                _ => 0,
            };
            let (i, j) = if analytic_rank(face_keys[j]) < analytic_rank(face_keys[i]) {
                (j, i)
            } else {
                (i, j)
            };

            let surf_i = match self.faces[face_keys[i]].surface.as_ref() {
                Some(s) => s,
                None => continue,
            };
            let surf_j = match self.faces[face_keys[j]].surface.as_ref() {
                Some(s) => s,
                None => continue,
            };

            // Analytic separation certificate (PK-329): exact O(1)
            // lower bound on untrimmed surface distance for analytic
            // class pairs — checked before any window/boundary work.
            // The offset pipeline's dominant verification class
            // (inner/outer skins of one source face: concentric
            // spheres, parallel cylinders a thickness apart) burned the
            // subdivision seeder's budget per pair proving what the
            // closed form states directly (measured 250 ms on the
            // sphere-capped rod's concentric cap pair).
            if let Some(lb) =
                crate::validation::self_intersection_pair_gate::analytic_separation_lower_bound(
                    &surf_i.geometry,
                    &surf_j.geometry,
                )
                && lb > linear_tol * 4.0
            {
                continue;
            }

            // Zero-padding trim windows (no false positives at shared
            // edges), pcurve-sampled so wrapped loops keep their full
            // period (face_uv_window — the positional window collapsed
            // on seam-anchored closed faces and this arm skipped every
            // such pair via the degenerate-range guard below), cached
            // per face above.
            let positions_i = &boundary_positions[i];
            let Some(((u0i, u1i), (v0i, v1i))) = windows[i] else {
                continue;
            };
            let positions_j = &boundary_positions[j];
            let Some(((u0j, u1j), (v0j, v1j))) = windows[j] else {
                continue;
            };

            // Skip faces with degenerate (zero-area) parameter ranges.
            if (u1i - u0i).abs() < f64::EPSILON || (v1i - v0i).abs() < f64::EPSILON {
                continue;
            }
            if (u1j - u0j).abs() < f64::EPSILON || (v1j - v0j).abs() < f64::EPSILON {
                continue;
            }

            // SAME-SURFACE patch exemption (PK-308) — see
            // `self_intersection_same_surface.rs`.
            if same_surface_domains_disjoint(
                &surf_i.geometry,
                &surf_j.geometry,
                positions_i,
                positions_j,
                &self.tolerance,
            ) {
                continue;
            }
            // Certified-separation pre-gate (PK-327): free-form pairs
            // route to the marching solver, which degrades badly on
            // near-parallel non-intersecting sheets (a healthy
            // thickened NURBS face: two skins a thickness apart).
            // Subdivision AABB-gap certificates prove separation
            // quickly; an unproven pair falls through unchanged.
            let free_form =
                crate::validation::self_intersection_trim::surface_uses_iterative_projection(
                    &surf_i.geometry,
                ) || crate::validation::self_intersection_trim::surface_uses_iterative_projection(
                    &surf_j.geometry,
                );
            if free_form
                && crate::validation::self_intersection_pair_gate::pair_clearly_separated(
                    &surf_i.geometry,
                    (u0i, u1i),
                    (v0i, v1i),
                    &surf_j.geometry,
                    (u0j, u1j),
                    (v0j, v1j),
                    linear_tol * 4.0,
                )
            {
                continue;
            }
            // Window-aware coaxial cone certificate (PK-329): coaxial
            // same-angle cones a thickness apart march ~25 ms/pair to
            // prove separation their matching trimmed nappes state in
            // closed form.
            if let Some(lb) =
                crate::validation::self_intersection_pair_gate::coaxial_cone_same_nappe_lower_bound(
                    &surf_i.geometry,
                    (v0i, v1i),
                    &surf_j.geometry,
                    (v0j, v1j),
                )
                && lb > linear_tol * 4.0
            {
                continue;
            }
            pairs_checked += 1;
            trace_self_intersection_candidate(
                self,
                pairs_checked,
                TraceSurfacePatch {
                    face: face_keys[i],
                    u_range: (u0i, u1i),
                    v_range: (v0i, v1i),
                },
                TraceSurfacePatch {
                    face: face_keys[j],
                    u_range: (u0j, u1j),
                    v_range: (v0j, v1j),
                },
            );
            jobs.push((i, j));
        }

        // Expensive phase (parallel, PK-329): the surviving jobs carry the
        // expensive per-pair work — bounded trim-overlap gates, the
        // one-sided-contact sampler, and the marching solver. Each is a
        // pure read of `self`, so jobs dispatch in parallel with an
        // ordered collect; the fold below consumes verdicts in job
        // order, so `found` and the `collect_all = false` first-hit
        // short-circuit are positionally identical to the sequential
        // loop at any thread count. Env-gated traces force serial
        // dispatch so their dumps cannot interleave.
        let evaluate = |&(i, j): &(usize, usize)| -> bool {
            let (Some(surf_i), Some(surf_j)) = (
                self.faces[face_keys[i]].surface.as_ref(),
                self.faces[face_keys[j]].surface.as_ref(),
            ) else {
                return false;
            };
            let (Some(((u0i, u1i), (v0i, v1i))), Some(((u0j, u1j), (v0j, v1j)))) =
                (windows[i], windows[j])
            else {
                return false;
            };
            if let Some(overlaps) = self.bounded_cylinder_cylinder_trim_overlap(
                face_keys[i],
                &surf_i.geometry,
                (u0i, u1i),
                (v0i, v1i),
                face_keys[j],
                &surf_j.geometry,
                (u0j, u1j),
                (v0j, v1j),
                linear_tol,
            ) {
                return overlaps;
            }
            // Cylinder × free-form validation fast path (PK-251): the
            // canonical ordering above guarantees the cylinder is
            // operand i for this class. Implicit-distance sampling
            // localizes the crossing curve where the marching solver's
            // adaptive step collapses on grazing fillet-seam contact
            // (measured 25 minutes for one validate_all without it).
            if let Some(overlaps) = self.bounded_analytic_free_form_trim_overlap(
                face_keys[i],
                &surf_i.geometry,
                face_keys[j],
                &surf_j.geometry,
                (u0j, u1j),
                (v0j, v1j),
                linear_tol,
            ) {
                return overlaps;
            }
            // Trim-scan context: one-sided tolerance-band contact (a
            // join strip riding its neighbour's skin) has no material
            // crossing and nothing to trim — skip the marching
            // escalation the VALIDATOR keeps for its coincidence
            // question (measured 36 s per grazing cylinder×NURBS pair).
            if tangential_band_is_clean
                && crate::validation::self_intersection_pair_gate::one_sided_contact_only(
                    &surf_i.geometry,
                    &surf_j.geometry,
                    (u0j, u1j),
                    (v0j, v1j),
                    linear_tol,
                )
            {
                return false;
            }
            let result = intersect_surfaces_with_handler(
                &surf_i.geometry,
                (u0i, u1i),
                (v0i, v1i),
                &surf_j.geometry,
                (u0j, u1j),
                (v0j, v1j),
                &coarse_tol,
                &DefaultDegeneracyHandler,
            );
            match result {
                Ok(ref r) => r.intersections.iter().any(|intersection| {
                    self.intersection_overlaps_face_trims(
                        face_keys[i],
                        &surf_i.geometry,
                        face_keys[j],
                        &surf_j.geometry,
                        intersection,
                        linear_tol,
                    )
                }),
                _ => false,
            }
        };
        let serial_for_trace = std::env::var_os("POPKERN_SELF_INTERSECTION_TRACE").is_some()
            || std::env::var_os("POPKERN_SI_CANDIDATE_TRACE").is_some();
        #[cfg(feature = "parallel")]
        let verdicts: Vec<bool> = if serial_for_trace {
            jobs.iter().map(evaluate).collect()
        } else {
            use rayon::prelude::*;
            jobs.par_iter().map(evaluate).collect()
        };
        #[cfg(not(feature = "parallel"))]
        let verdicts: Vec<bool> = {
            let _ = serial_for_trace;
            jobs.iter().map(evaluate).collect()
        };

        for (&(i, j), hit) in jobs.iter().zip(verdicts) {
            if !hit {
                continue;
            }
            if std::env::var_os("POPKERN_SELF_INTERSECTION_TRACE").is_some() {
                self.trace_self_intersection_pair(face_keys[i], face_keys[j]);
            }
            found.push((face_keys[i], face_keys[j]));
            if !collect_all {
                return found;
            }
        }

        found
    }
}
