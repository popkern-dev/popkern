//! Self-intersection validation check (#11).
//!
//! Detects when non-adjacent faces in a solid intersect each other.

use crate::error::ValidationError;
use crate::geometry::intersection::degeneracy::DefaultDegeneracyHandler;
use crate::geometry::intersection::intersect_surfaces_with_handler;
use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::surface::SurfaceEval;
use crate::tolerance::Tolerance;
use crate::topology::BRep;
use crate::validation::bvh::Bvh;
use std::collections::HashSet;

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
        let n = face_keys.len();
        if n < 2 {
            return Ok(());
        }

        // ------------------------------------------------------------------ //
        // 1. Build adjacency set (by face index pairs, not handles)           //
        // ------------------------------------------------------------------ //
        let mut adjacent: HashSet<(usize, usize)> = HashSet::new();

        for (i, &fk) in face_keys.iter().enumerate() {
            let face = &self.faces[fk];
            let start_he_handle = self.loops[face.outer_loop].half_edge;
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
                        crate::geometry::uv_projection::bounded_parameter_range(
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
            return Ok(());
        }

        let bvh = Bvh::build(&bvh_items);

        // Use 2× linear tolerance for intersection detection: tight enough to catch
        // genuine self-intersections, loose enough for fast intersection testing.
        // This is conservative: may report near-misses as intersections (false positives)
        // but won't miss actual intersections within tolerance (no false negatives).
        let linear_tol = self.tolerance.linear;
        let coarse_tol = Tolerance::default()
            .with_linear(linear_tol * 2.0)
            .with_angular(self.tolerance.angular)
            .with_parametric(self.tolerance.parametric);

        // ------------------------------------------------------------------ //
        // 3. Test non-adjacent, AABB-overlapping pairs with intersect_surfaces //
        // ------------------------------------------------------------------ //
        let overlapping = bvh.find_overlapping_pairs(linear_tol);

        for (i, j) in overlapping {
            let pair = if i < j { (i, j) } else { (j, i) };
            if adjacent.contains(&pair) {
                continue;
            }

            let surf_i = match self.faces[face_keys[i]].surface.as_ref() {
                Some(s) => s,
                None => continue,
            };
            let surf_j = match self.faces[face_keys[j]].surface.as_ref() {
                Some(s) => s,
                None => continue,
            };

            let positions_i = self.face_outer_loop_positions(face_keys[i]);
            let ((u0i, u1i), (v0i, v1i)) = crate::geometry::uv_projection::bounded_parameter_range(
                &positions_i,
                &surf_i.geometry,
            );
            let positions_j = self.face_outer_loop_positions(face_keys[j]);
            let ((u0j, u1j), (v0j, v1j)) = crate::geometry::uv_projection::bounded_parameter_range(
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
                Ok(ref r) if !r.intersections.is_empty() => {
                    return Err(ValidationError::SelfIntersection {
                        face_a: face_keys[i],
                        face_b: face_keys[j],
                    });
                }
                _ => {}
            }
        }

        Ok(())
    }
}

#[cfg(test)]
#[path = "self_intersection_tests.rs"]
mod tests;
