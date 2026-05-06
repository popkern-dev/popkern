//! Transactional adjacent-face merge.

use crate::error::CoreError;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle};

/// Result of merging two adjacent faces across one shared edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FaceMergeResult {
    /// Face that remains after the merge.
    pub retained_face: FaceHandle,
    /// Face removed by the merge.
    pub removed_face: FaceHandle,
    /// Edge removed by the merge.
    pub removed_edge: EdgeHandle,
    /// Loop that remains after the merge.
    pub retained_loop: LoopHandle,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct MergeCounts {
    vertices: i64,
    edges: i64,
    faces: i64,
    loops: i64,
}

impl BRep {
    /// Merge two adjacent faces by deleting their shared edge.
    ///
    /// This wraps `kef` in a transaction and checks the expected
    /// `dV=0 dE=-1 dF=-1 dL=-1` delta before commit. It only accepts shared
    /// edges whose two half-edges bound distinct faces in the same shell, and
    /// refuses incompatible analytic surfaces instead of silently creating a
    /// geometrically ambiguous face.
    pub fn merge_faces_across_edge(
        &mut self,
        half_edge: HalfEdgeHandle,
    ) -> Result<FaceMergeResult, CoreError> {
        let before = MergeCounts::capture(self);
        let mut candidate = self.clone();
        let result = candidate.apply_merge_faces_across_edge(half_edge)?;
        candidate.verify_face_merge_postconditions(before, result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_merge_faces_across_edge(
        &mut self,
        half_edge: HalfEdgeHandle,
    ) -> Result<FaceMergeResult, CoreError> {
        let data = self
            .half_edges
            .get(half_edge)
            .ok_or_else(|| CoreError::InvalidHandle("merge half-edge is invalid".into()))?;
        let twin = self
            .half_edges
            .get(data.twin)
            .ok_or_else(|| CoreError::InvalidHandle("merge twin half-edge is invalid".into()))?;
        let retained_loop = data.loop_;
        let removed_loop = twin.loop_;
        if retained_loop == removed_loop {
            return Err(CoreError::PreconditionFailed(
                "merge edge must separate two face loops".into(),
            ));
        }
        let retained_face = self.loops[retained_loop].face;
        let removed_face = self.loops[removed_loop].face;
        if retained_face == removed_face {
            return Err(CoreError::PreconditionFailed(
                "merge edge must separate two faces".into(),
            ));
        }
        if self.faces[retained_face].shell != self.faces[removed_face].shell {
            return Err(CoreError::PreconditionFailed(
                "merge faces must belong to the same shell".into(),
            ));
        }
        if !merge_surfaces_are_compatible(
            self.faces[retained_face].surface.as_ref(),
            self.faces[removed_face].surface.as_ref(),
            self.tolerance.linear(),
        ) {
            return Err(CoreError::PreconditionFailed(
                "merge faces require compatible coplanar surfaces".into(),
            ));
        }
        let removed_edge = data.edge;
        self.kef(half_edge)?;
        Ok(FaceMergeResult {
            retained_face,
            removed_face,
            removed_edge,
            retained_loop,
        })
    }

    fn verify_face_merge_postconditions(
        &self,
        before: MergeCounts,
        result: FaceMergeResult,
    ) -> Result<(), CoreError> {
        let actual_delta = MergeCounts::capture(self).delta_from(before);
        let expected_delta = MergeCounts {
            vertices: 0,
            edges: -1,
            faces: -1,
            loops: -1,
        };
        if actual_delta != expected_delta {
            return Err(CoreError::PreconditionFailed(format!(
                "face merge postcondition failed: expected_delta=dV=0 dE=-1 dF=-1 dL=-1 actual_delta=dV={} dE={} dF={} dL={}",
                actual_delta.vertices, actual_delta.edges, actual_delta.faces, actual_delta.loops
            )));
        }
        if self.faces.contains_key(result.removed_face)
            || self.loops.contains_key(result.retained_loop)
                && self.loops[result.retained_loop].face != result.retained_face
        {
            return Err(CoreError::PreconditionFailed(
                "face merge retained/removed ownership postcondition failed".into(),
            ));
        }
        if self.edges.contains_key(result.removed_edge) {
            return Err(CoreError::PreconditionFailed(
                "face merge removed edge still exists".into(),
            ));
        }
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)?;
        Ok(())
    }
}

fn merge_surfaces_are_compatible(
    retained: Option<&SurfaceOnFace>,
    removed: Option<&SurfaceOnFace>,
    tolerance: f64,
) -> bool {
    match (retained, removed) {
        (Some(a), Some(b)) => plane_surfaces_are_coplanar(&a.geometry, &b.geometry, tolerance),
        _ => true,
    }
}

fn plane_surfaces_are_coplanar(
    first: &SurfaceGeometry,
    second: &SurfaceGeometry,
    tolerance: f64,
) -> bool {
    match (first, second) {
        (
            SurfaceGeometry::Plane {
                origin: first_origin,
                normal: first_normal,
                ..
            },
            SurfaceGeometry::Plane {
                origin: second_origin,
                normal: second_normal,
                ..
            },
        ) => {
            let normals_parallel = first_normal.cross(second_normal).length() <= tolerance;
            let origin_delta = *second_origin - *first_origin;
            normals_parallel && origin_delta.dot(first_normal).abs() <= tolerance
        }
        _ => false,
    }
}

impl MergeCounts {
    fn capture(brep: &BRep) -> Self {
        Self {
            vertices: brep.vertices.len() as i64,
            edges: brep.edges.len() as i64,
            faces: brep.faces.len() as i64,
            loops: brep.loops.len() as i64,
        }
    }

    fn delta_from(self, before: Self) -> Self {
        Self {
            vertices: self.vertices - before.vertices,
            edges: self.edges - before.edges,
            faces: self.faces - before.faces,
            loops: self.loops - before.loops,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::tolerance::Tolerance;
    use crate::topology::{BRep, FaceHandle};

    #[test]
    fn merge_faces_across_edge_commits_for_compatible_adjacent_faces() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let face = brep.faces.keys().next().expect("cube has a face");
        let half_edge = brep.loops[brep.faces[face].outer_loop].half_edge;
        let shared = brep.half_edges[half_edge].twin;
        for face in brep.faces.values_mut() {
            face.surface = None;
        }
        let before = (
            brep.vertices.len(),
            brep.edges.len(),
            brep.faces.len(),
            brep.loops.len(),
        );

        let result = brep
            .merge_faces_across_edge(shared)
            .expect("topology-only adjacent faces merge");

        assert_eq!(
            (
                brep.vertices.len() as i64 - before.0 as i64,
                brep.edges.len() as i64 - before.1 as i64,
                brep.faces.len() as i64 - before.2 as i64,
                brep.loops.len() as i64 - before.3 as i64
            ),
            (0, -1, -1, -1)
        );
        assert!(!brep.faces.contains_key(result.removed_face));
        assert!(!brep.edges.contains_key(result.removed_edge));
        brep.validate_half_edge_pairing().unwrap();
        brep.validate_half_edge_cycles().unwrap();
        brep.validate_manifold_edges().unwrap();
    }

    #[test]
    fn merge_faces_across_edge_rolls_back_for_incompatible_cube_planes() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let (first, second) = adjacent_face_pair(&brep);
        let half_edge = shared_half_edge(&brep, first, second);
        let before = (brep.edges.len(), brep.faces.len(), brep.loops.len());

        let err = brep
            .merge_faces_across_edge(half_edge)
            .expect_err("perpendicular cube faces are not coplanar");

        assert!(err.to_string().contains("compatible coplanar"));
        assert_eq!(
            before,
            (brep.edges.len(), brep.faces.len(), brep.loops.len())
        );
    }

    fn adjacent_face_pair(brep: &BRep) -> (FaceHandle, FaceHandle) {
        let faces = brep.faces.keys().collect::<Vec<_>>();
        for (i, first) in faces.iter().enumerate() {
            for second in faces.iter().skip(i + 1) {
                if shared_half_edge_optional(brep, *first, *second).is_some() {
                    return (*first, *second);
                }
            }
        }
        panic!("cube fixture has no adjacent face pair")
    }

    fn shared_half_edge(
        brep: &BRep,
        first: FaceHandle,
        second: FaceHandle,
    ) -> crate::topology::HalfEdgeHandle {
        shared_half_edge_optional(brep, first, second).expect("faces share an edge")
    }

    fn shared_half_edge_optional(
        brep: &BRep,
        first: FaceHandle,
        second: FaceHandle,
    ) -> Option<crate::topology::HalfEdgeHandle> {
        let first_loop = brep.faces[first].outer_loop;
        let mut cursor = brep.loops[first_loop].half_edge;
        loop {
            let twin = brep.half_edges[cursor].twin;
            if brep
                .half_edges
                .get(twin)
                .is_some_and(|data| brep.loops[data.loop_].face == second)
            {
                return Some(cursor);
            }
            cursor = brep.half_edges[cursor].next;
            if cursor == brep.loops[first_loop].half_edge {
                break;
            }
        }
        None
    }
}
