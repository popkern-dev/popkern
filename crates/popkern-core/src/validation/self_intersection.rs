//! Self-intersection validation check (#11).
//!
//! Detects when non-adjacent faces in a solid intersect each other.

use crate::error::ValidationError;
use crate::geometry::intersection::degeneracy::DefaultDegeneracyHandler;
use crate::geometry::intersection::intersect_surfaces_with_handler;
use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::surface::SurfaceEval;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle};
use crate::validation::bvh::Bvh;
pub(super) use self_intersection_helpers::{
    AnalyticFreeFormSample, analytic_free_form_near_samples,
};
use self_intersection_helpers::{
    TraceSurfacePatch, trace_self_intersection_candidate, trace_surface_kind,
};
use std::collections::HashSet;

#[path = "self_intersection_helpers.rs"]
mod self_intersection_helpers;
#[path = "self_intersection_pair_fast_paths.rs"]
mod self_intersection_pair_fast_paths;
#[path = "self_intersection_same_surface.rs"]
mod self_intersection_same_surface;
use self_intersection_same_surface::same_surface_domains_disjoint;

impl BRep {
    /// Check #11: no face intersects another non-adjacent face.
    ///
    /// The check proceeds in three stages:
    ///
    /// 1. Build a set of adjacent face-index pairs by walking shared half-edges.
    /// 2. Sample each face surface on a 5×5 grid to compute an AABB; build a
    ///    BVH from those AABBs for efficient spatial queries.
    /// 3. For each non-adjacent pair whose AABBs overlap (from BVH query), call
    ///    `intersect_surfaces_with_handler` with `DefaultDegeneracyHandler` and
    ///    a coarse tolerance (2× linear) to test for actual surface intersection.
    ///    Any intersections found are reported as a `SelfIntersection` error.
    pub fn validate_self_intersection(&self) -> Result<(), ValidationError> {
        let face_keys: Vec<_> = self.faces.keys().collect();
        self.validate_self_intersection_among(&face_keys)
    }

    /// As [`Self::validate_self_intersection`], restricted to pairs
    /// drawn from `face_keys` (the same adjacency exemption, AABB/BVH
    /// pruning, and pair budget). Operators that mutate a bounded
    /// neighbourhood use this as a commit postcondition: checking the
    /// touched faces alone keeps the gate O(neighbourhood²) instead
    /// of O(model²) — the full-model form was measured at 250 s on a
    /// single N=3 conic corner (NURBS strip pairs dominate).
    /// Unknown or dead handles are ignored.
    pub fn validate_self_intersection_among(
        &self,
        face_keys: &[FaceHandle],
    ) -> Result<(), ValidationError> {
        match self
            .collect_self_intersecting_pairs_among(face_keys, false, false)
            .first()
        {
            Some(&(face_a, face_b)) => Err(ValidationError::SelfIntersection { face_a, face_b }),
            None => Ok(()),
        }
    }

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

        let aabbs: Vec<Option<Aabb>> = face_keys
            .iter()
            .map(|&fk| {
                let face = &self.faces[fk];
                face.surface.as_ref().and_then(|s| {
                    let positions = self.face_outer_loop_positions(fk);
                    let ((u0, u1), (v0, v1)) =
                        crate::geometry::uv_projection::tight_bounded_parameter_range(
                            &positions,
                            &s.geometry,
                        );
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

            // Use exact (zero-padding) ranges for intersection testing to
            // avoid false positives at shared edges where padding from two
            // perpendicular faces would overlap.
            let positions_i = self.face_outer_loop_positions(face_keys[i]);
            let ((u0i, u1i), (v0i, v1i)) =
                crate::geometry::uv_projection::exact_bounded_parameter_range(
                    &positions_i,
                    &surf_i.geometry,
                );
            let positions_j = self.face_outer_loop_positions(face_keys[j]);
            let ((u0j, u1j), (v0j, v1j)) =
                crate::geometry::uv_projection::exact_bounded_parameter_range(
                    &positions_j,
                    &surf_j.geometry,
                );

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
                &positions_i,
                &positions_j,
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
                if overlaps {
                    found.push((face_keys[i], face_keys[j]));
                    if !collect_all {
                        return found;
                    }
                }
                continue;
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
                if overlaps {
                    found.push((face_keys[i], face_keys[j]));
                    if !collect_all {
                        return found;
                    }
                }
                continue;
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
                continue;
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
                Ok(ref r)
                    if r.intersections.iter().any(|intersection| {
                        self.intersection_overlaps_face_trims(
                            face_keys[i],
                            &surf_i.geometry,
                            face_keys[j],
                            &surf_j.geometry,
                            intersection,
                            linear_tol,
                        )
                    }) =>
                {
                    if std::env::var_os("POPKERN_SELF_INTERSECTION_TRACE").is_some() {
                        self.trace_self_intersection_pair(face_keys[i], face_keys[j]);
                    }
                    found.push((face_keys[i], face_keys[j]));
                    if !collect_all {
                        return found;
                    }
                }
                _ => {}
            }
        }

        found
    }

    fn trace_self_intersection_pair(
        &self,
        face_a: crate::topology::FaceHandle,
        face_b: crate::topology::FaceHandle,
    ) {
        eprintln!(
            "[self-intersection pair] {face_a:?}/{} {face_b:?}/{}",
            trace_surface_kind(self, face_a),
            trace_surface_kind(self, face_b),
        );
        for face in [face_a, face_b] {
            let Some(face_data) = self.faces.get(face) else {
                continue;
            };
            let Some(loop_data) = self.loops.get(face_data.outer_loop) else {
                continue;
            };
            let mut entries = Vec::new();
            let mut cursor = loop_data.half_edge;
            for _ in 0..=self.half_edges.len() {
                let Some(he) = self.half_edges.get(cursor) else {
                    break;
                };
                let p = self.vertices[he.origin].position;
                entries.push(format!(
                    "{cursor:?}:e{:?}:v{:?}=({:.4},{:.4},{:.4})",
                    he.edge, he.origin, p.x, p.y, p.z,
                ));
                cursor = he.next;
                if cursor == loop_data.half_edge {
                    break;
                }
            }
            eprintln!(
                "[self-intersection loop] {face:?}: {}",
                entries.join(" -> ")
            );
        }
    }
}

#[cfg(test)]
#[path = "self_intersection_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "self_intersection_order_tests.rs"]
mod order_tests;

#[cfg(test)]
#[path = "self_intersection_grazing_tests.rs"]
mod grazing_tests;
