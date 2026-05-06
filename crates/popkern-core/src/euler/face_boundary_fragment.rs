//! Retire a boundary fragment from one face to expose its opposite interval.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle};

/// Result of splitting off and retiring a face-boundary fragment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundaryFragmentRetireResult {
    /// Face that kept the complementary boundary interval.
    pub retained_face: FaceHandle,
    /// Face fragment removed from the topology.
    pub retired_face: FaceHandle,
    /// New edge closing the retained face after the fragment split.
    pub replacement_boundary_edge: EdgeHandle,
    /// Existing edges whose opposite half-edge became boundary.
    pub exposed_edges: Vec<EdgeHandle>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct FragmentCounts {
    vertices: i64,
    edges: i64,
    faces: i64,
    loops: i64,
}

impl BRep {
    /// Split off a contiguous outer-loop interval and retire that fragment.
    ///
    /// The interval is half-open: `first_half_edge` is included and
    /// `after_last_half_edge` is the first half-edge kept by the retained face.
    /// The operator is transactional. It uses `mef` to create a fragment face,
    /// retires that whole fragment, and commits only after the expected
    /// `dV=0 dE=1 dF=0 dL=0` delta and structural validation pass.
    ///
    /// This exposes the twins of the retired interval as boundary edges and
    /// leaves one replacement boundary edge on the retained face. Callers that
    /// need a closed shell must plan subsequent ownership for that replacement
    /// boundary explicitly.
    pub fn retire_face_boundary_fragment(
        &mut self,
        face: FaceHandle,
        first_half_edge: HalfEdgeHandle,
        after_last_half_edge: HalfEdgeHandle,
    ) -> Result<BoundaryFragmentRetireResult, CoreError> {
        let before = FragmentCounts::capture(self);
        let mut candidate = self.clone();
        let result = candidate.apply_retire_face_boundary_fragment(
            face,
            first_half_edge,
            after_last_half_edge,
        )?;
        candidate.verify_boundary_fragment_postconditions(
            before,
            &result,
            first_half_edge,
            after_last_half_edge,
        )?;
        *self = candidate;
        Ok(result)
    }

    fn apply_retire_face_boundary_fragment(
        &mut self,
        face: FaceHandle,
        first_half_edge: HalfEdgeHandle,
        after_last_half_edge: HalfEdgeHandle,
    ) -> Result<BoundaryFragmentRetireResult, CoreError> {
        let interval =
            self.boundary_fragment_interval(face, first_half_edge, after_last_half_edge)?;
        let interval_edges = interval
            .iter()
            .map(|half_edge| self.half_edges[*half_edge].edge)
            .collect::<Vec<_>>();
        let split = self.mef(first_half_edge, after_last_half_edge)?;
        let retained_face = self.loops[self.half_edges[split.half_edge].loop_].face;
        if retained_face != face {
            return Err(CoreError::PreconditionFailed(
                "boundary fragment split did not retain the requested face".into(),
            ));
        }
        self.retire_faces_to_boundary(&[split.face])?;

        let mut exposed_edges = Vec::new();
        for edge in interval_edges {
            if self.boundary_edges.contains(&edge) {
                exposed_edges.push(edge);
            }
        }
        Ok(BoundaryFragmentRetireResult {
            retained_face,
            retired_face: split.face,
            replacement_boundary_edge: split.edge,
            exposed_edges,
        })
    }

    fn boundary_fragment_interval(
        &self,
        face: FaceHandle,
        first_half_edge: HalfEdgeHandle,
        after_last_half_edge: HalfEdgeHandle,
    ) -> Result<Vec<HalfEdgeHandle>, CoreError> {
        let face_data = self.faces.get(face).ok_or_else(|| {
            CoreError::InvalidHandle("boundary fragment face handle is invalid".into())
        })?;
        if !face_data.inner_loops.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "boundary fragment retire does not support faces with inner loops".into(),
            ));
        }
        let loop_handle = face_data.outer_loop;
        let first = self.half_edges.get(first_half_edge).ok_or_else(|| {
            CoreError::InvalidHandle("boundary fragment first half-edge is invalid".into())
        })?;
        let after_last = self.half_edges.get(after_last_half_edge).ok_or_else(|| {
            CoreError::InvalidHandle("boundary fragment after-last half-edge is invalid".into())
        })?;
        if first.loop_ != loop_handle || after_last.loop_ != loop_handle {
            return Err(CoreError::PreconditionFailed(
                "boundary fragment endpoints must be on the requested face outer loop".into(),
            ));
        }
        if first_half_edge == after_last_half_edge {
            return Err(CoreError::PreconditionFailed(
                "boundary fragment interval must not cover the whole loop".into(),
            ));
        }

        let mut cursor = first_half_edge;
        let mut interval = Vec::new();
        for _ in 0..=self.half_edges.len() {
            if cursor == after_last_half_edge {
                if interval.is_empty() {
                    return Err(CoreError::PreconditionFailed(
                        "boundary fragment interval is empty".into(),
                    ));
                }
                return Ok(interval);
            }
            let half_edge = self.half_edges.get(cursor).ok_or_else(|| {
                CoreError::InvalidHandle("boundary fragment interval half-edge is invalid".into())
            })?;
            if half_edge.loop_ != loop_handle {
                return Err(CoreError::PreconditionFailed(
                    "boundary fragment interval left the requested face loop".into(),
                ));
            }
            if self.boundary_edges.contains(&half_edge.edge) {
                return Err(CoreError::PreconditionFailed(
                    "boundary fragment interval already contains a boundary edge".into(),
                ));
            }
            interval.push(cursor);
            cursor = half_edge.next;
        }
        Err(CoreError::PreconditionFailed(
            "boundary fragment interval walk did not reach after-last half-edge".into(),
        ))
    }

    fn verify_boundary_fragment_postconditions(
        &self,
        before: FragmentCounts,
        result: &BoundaryFragmentRetireResult,
        _first_half_edge: HalfEdgeHandle,
        after_last_half_edge: HalfEdgeHandle,
    ) -> Result<(), CoreError> {
        let actual_delta = FragmentCounts::capture(self).delta_from(before);
        let expected_delta = FragmentCounts {
            vertices: 0,
            edges: 1,
            faces: 0,
            loops: 0,
        };
        if actual_delta != expected_delta {
            return Err(CoreError::PreconditionFailed(format!(
                "boundary fragment postcondition failed: expected_delta=dV=0 dE=1 dF=0 dL=0 actual_delta=dV={} dE={} dF={} dL={}",
                actual_delta.vertices, actual_delta.edges, actual_delta.faces, actual_delta.loops
            )));
        }
        if self.faces.contains_key(result.retired_face) {
            return Err(CoreError::PreconditionFailed(
                "boundary fragment retired face still exists after commit".into(),
            ));
        }
        if !self
            .boundary_edges
            .contains(&result.replacement_boundary_edge)
        {
            return Err(CoreError::PreconditionFailed(
                "boundary fragment replacement edge is not boundary".into(),
            ));
        }
        for edge in &result.exposed_edges {
            let exposed_half_edge = self.edges[*edge].half_edge;
            if self.half_edges[exposed_half_edge].twin != HalfEdgeHandle::default() {
                return Err(CoreError::PreconditionFailed(
                    "boundary fragment exposed edge kept a live twin".into(),
                ));
            }
        }
        if self.half_edges[after_last_half_edge].loop_
            != self.faces[result.retained_face].outer_loop
        {
            return Err(CoreError::PreconditionFailed(
                "boundary fragment retained loop does not own after-last half-edge".into(),
            ));
        }
        let retained_loop = self.faces[result.retained_face].outer_loop;
        let replacement_half_edge = self
            .edges
            .get(result.replacement_boundary_edge)
            .ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "boundary fragment replacement edge no longer exists".into(),
                )
            })?
            .half_edge;
        self.verify_boundary_fragment_retained_loop(
            retained_loop,
            &[after_last_half_edge, replacement_half_edge],
        )?;
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)?;
        Ok(())
    }

    fn verify_boundary_fragment_retained_loop(
        &self,
        loop_: LoopHandle,
        required_half_edges: &[HalfEdgeHandle],
    ) -> Result<(), CoreError> {
        let loop_data = self.loops.get(loop_).ok_or_else(|| {
            CoreError::PreconditionFailed(
                "boundary fragment retained face references a retired loop".into(),
            )
        })?;
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut visited = BTreeSet::new();
        for _ in 0..=self.half_edges.len() {
            let half_edge = self.half_edges.get(cursor).ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "boundary fragment retained loop references a missing half-edge".into(),
                )
            })?;
            if half_edge.loop_ != loop_ {
                return Err(CoreError::PreconditionFailed(
                    "boundary fragment retained cycle left the retained loop".into(),
                ));
            }
            if !visited.insert(cursor) {
                return Err(CoreError::PreconditionFailed(
                    "boundary fragment retained cycle repeats before closing".into(),
                ));
            }
            cursor = half_edge.next;
            if cursor == start {
                if required_half_edges
                    .iter()
                    .all(|half_edge| visited.contains(half_edge))
                {
                    return Ok(());
                }
                return Err(CoreError::PreconditionFailed(
                    "boundary fragment retained cycle is missing a required half-edge".into(),
                ));
            }
        }
        Err(CoreError::PreconditionFailed(
            "boundary fragment retained cycle did not close".into(),
        ))
    }
}

impl FragmentCounts {
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
    use crate::topology::BRep;

    #[test]
    fn retire_face_boundary_fragment_exposes_interval_transactionally() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let face = brep.faces.keys().next().expect("cube has a face");
        let loop_ = brep.faces[face].outer_loop;
        let first = brep.loops[loop_].half_edge;
        let after_last = brep.half_edges[brep.half_edges[first].next].next;
        let first_edge = brep.half_edges[first].edge;
        let second_edge = brep.half_edges[brep.half_edges[first].next].edge;
        let before = (
            brep.vertices.len(),
            brep.edges.len(),
            brep.faces.len(),
            brep.loops.len(),
        );

        let result = brep
            .retire_face_boundary_fragment(face, first, after_last)
            .expect("fragment retires");

        assert_eq!(
            (
                brep.vertices.len() as i64 - before.0 as i64,
                brep.edges.len() as i64 - before.1 as i64,
                brep.faces.len() as i64 - before.2 as i64,
                brep.loops.len() as i64 - before.3 as i64
            ),
            (0, 1, 0, 0)
        );
        assert_eq!(result.retained_face, face);
        assert!(!brep.faces.contains_key(result.retired_face));
        assert!(brep.boundary_edges.contains(&first_edge));
        assert!(brep.boundary_edges.contains(&second_edge));
        assert!(
            brep.boundary_edges
                .contains(&result.replacement_boundary_edge)
        );
        brep.validate().unwrap();
    }

    #[test]
    fn retire_face_boundary_fragment_rolls_back_for_empty_interval() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let face = brep.faces.keys().next().expect("cube has a face");
        let first = brep.loops[brep.faces[face].outer_loop].half_edge;
        let before = (brep.edges.len(), brep.faces.len(), brep.loops.len());

        let err = brep
            .retire_face_boundary_fragment(face, first, first)
            .expect_err("whole-loop interval is rejected");

        assert!(err.to_string().contains("whole loop"));
        assert_eq!(
            before,
            (brep.edges.len(), brep.faces.len(), brep.loops.len())
        );
    }
}
