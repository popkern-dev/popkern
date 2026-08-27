//! Self-intersection validation check (#11).
//!
//! Detects when non-adjacent faces in a solid intersect each other.

use crate::error::ValidationError;
use crate::topology::{BRep, FaceHandle};
use self_intersection_helpers::trace_surface_kind;
pub(super) use self_intersection_helpers::{
    AnalyticFreeFormSample, analytic_free_form_near_samples,
};

#[path = "self_intersection_helpers.rs"]
mod self_intersection_helpers;
#[path = "self_intersection_pair_fast_paths.rs"]
mod self_intersection_pair_fast_paths;
#[path = "self_intersection_pair_scan.rs"]
mod self_intersection_pair_scan;
#[path = "self_intersection_same_surface.rs"]
mod self_intersection_same_surface;

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
