//! Transactional adjacent-face merge followed by retained-loop face splits.

use crate::error::CoreError;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::{FaceMergeResult, MefResult};

/// One retained-loop split requested after an adjacent-face merge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FaceMergeSplitRequest {
    /// Start vertex for the split edge.
    pub start: VertexHandle,
    /// End vertex for the split edge.
    pub end: VertexHandle,
}

/// Result of a transactional merge followed by retained-loop splits.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaceMergeSplitResult {
    /// The adjacent-face merge result.
    pub merge: FaceMergeResult,
    /// The split results, in request order.
    pub splits: Vec<MefResult>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct MergeSplitCounts {
    vertices: i64,
    edges: i64,
    faces: i64,
    loops: i64,
}

impl BRep {
    /// Merge adjacent faces, then split the retained outer loop.
    ///
    /// This is a transactional topology operator for workflows that first
    /// reclassify two adjacent residual faces into one support region and then
    /// author replacement side-trim topology on that retained region. All split
    /// requests are resolved against the retained face's current outer loop in
    /// sequence. The operator commits only after validating the combined
    /// `dV=0 dE=-1+n dF=-1+n dL=-1+n` delta and structural topology checks.
    pub fn merge_faces_across_edge_then_split(
        &mut self,
        shared_half_edge: HalfEdgeHandle,
        split_requests: &[FaceMergeSplitRequest],
    ) -> Result<FaceMergeSplitResult, CoreError> {
        let before = MergeSplitCounts::capture(self);
        let mut candidate = self.clone();
        let result =
            candidate.apply_merge_faces_across_edge_then_split(shared_half_edge, split_requests)?;
        candidate.verify_merge_split_postconditions(before, &result, split_requests.len())?;
        *self = candidate;
        Ok(result)
    }

    fn apply_merge_faces_across_edge_then_split(
        &mut self,
        shared_half_edge: HalfEdgeHandle,
        split_requests: &[FaceMergeSplitRequest],
    ) -> Result<FaceMergeSplitResult, CoreError> {
        let merge = self.merge_faces_across_edge(shared_half_edge)?;
        let mut splits = Vec::with_capacity(split_requests.len());
        for request in split_requests {
            let retained_loop = self.retained_outer_loop(merge.retained_face)?;
            let start_half_edge =
                self.find_vertex_half_edge_on_loop(retained_loop, request.start)?;
            let end_half_edge = self.find_vertex_half_edge_on_loop(retained_loop, request.end)?;
            splits.push(self.mef(start_half_edge, end_half_edge)?);
        }
        Ok(FaceMergeSplitResult { merge, splits })
    }

    fn verify_merge_split_postconditions(
        &self,
        before: MergeSplitCounts,
        result: &FaceMergeSplitResult,
        split_count: usize,
    ) -> Result<(), CoreError> {
        let actual_delta = MergeSplitCounts::capture(self).delta_from(before);
        let expected_delta = MergeSplitCounts {
            vertices: 0,
            edges: -1 + split_count as i64,
            faces: -1 + split_count as i64,
            loops: -1 + split_count as i64,
        };
        if actual_delta != expected_delta {
            return Err(CoreError::PreconditionFailed(format!(
                "merge split postcondition failed: expected_delta=dV=0 dE={} dF={} dL={} actual_delta=dV={} dE={} dF={} dL={}",
                expected_delta.edges,
                expected_delta.faces,
                expected_delta.loops,
                actual_delta.vertices,
                actual_delta.edges,
                actual_delta.faces,
                actual_delta.loops
            )));
        }
        if self.faces.contains_key(result.merge.removed_face)
            || self.edges.contains_key(result.merge.removed_edge)
        {
            return Err(CoreError::PreconditionFailed(
                "merge split removed merge operands still exist".into(),
            ));
        }
        for split in &result.splits {
            if !self.edges.contains_key(split.edge)
                || !self.faces.contains_key(split.face)
                || !self.loops.contains_key(split.loop_)
            {
                return Err(CoreError::PreconditionFailed(
                    "merge split created split operand is missing".into(),
                ));
            }
        }
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)?;
        Ok(())
    }

    fn retained_outer_loop(&self, face: FaceHandle) -> Result<LoopHandle, CoreError> {
        self.faces
            .get(face)
            .map(|face_data| face_data.outer_loop)
            .ok_or_else(|| CoreError::InvalidHandle("retained face is invalid".into()))
    }

    fn find_vertex_half_edge_on_loop(
        &self,
        loop_handle: LoopHandle,
        vertex: VertexHandle,
    ) -> Result<HalfEdgeHandle, CoreError> {
        let start = self
            .loops
            .get(loop_handle)
            .ok_or_else(|| CoreError::InvalidHandle("retained loop is invalid".into()))?
            .half_edge;
        let mut cursor = start;
        for _ in 0..=self.half_edges.len() {
            let half_edge = self
                .half_edges
                .get(cursor)
                .ok_or_else(|| CoreError::InvalidHandle("retained half-edge is invalid".into()))?;
            if half_edge.loop_ != loop_handle {
                return Err(CoreError::PreconditionFailed(
                    "retained loop contains half-edge from another loop".into(),
                ));
            }
            if half_edge.origin == vertex {
                return Ok(cursor);
            }
            cursor = half_edge.next;
            if cursor == start {
                break;
            }
        }
        Err(CoreError::PreconditionFailed(format!(
            "vertex {vertex:?} is not on retained merge loop {loop_handle:?}"
        )))
    }
}

impl MergeSplitCounts {
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
    use crate::euler::FaceMergeSplitRequest;
    use crate::tolerance::Tolerance;
    use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

    #[test]
    fn merge_faces_across_edge_then_split_commits_transactionally() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        for face in brep.faces.values_mut() {
            face.surface = None;
        }
        let (first, second) = adjacent_face_pair(&brep);
        let shared = shared_half_edge(&brep, first, second);
        let mut probe = brep.clone();
        let merge = probe
            .merge_faces_across_edge(shared)
            .expect("topology-only faces merge");
        let vertices = loop_vertices(&probe, merge.retained_loop);
        let requests = [FaceMergeSplitRequest {
            start: vertices[0],
            end: vertices[2],
        }];
        let before = (
            brep.vertices.len(),
            brep.edges.len(),
            brep.faces.len(),
            brep.loops.len(),
        );

        let result = brep
            .merge_faces_across_edge_then_split(shared, &requests)
            .expect("merge plus retained-loop split commits");

        assert_eq!(result.splits.len(), 1);
        assert_eq!(
            (
                brep.vertices.len() as i64 - before.0 as i64,
                brep.edges.len() as i64 - before.1 as i64,
                brep.faces.len() as i64 - before.2 as i64,
                brep.loops.len() as i64 - before.3 as i64
            ),
            (0, 0, 0, 0)
        );
        brep.validate_half_edge_pairing().unwrap();
        brep.validate_half_edge_cycles().unwrap();
        brep.validate_manifold_edges().unwrap();
    }

    #[test]
    fn merge_faces_across_edge_then_split_rolls_back_for_missing_vertex() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        for face in brep.faces.values_mut() {
            face.surface = None;
        }
        let (first, second) = adjacent_face_pair(&brep);
        let shared = shared_half_edge(&brep, first, second);
        let missing = brep.vertices.insert(crate::topology::Vertex {
            position: crate::geometry::Point3::new(99.0, 99.0, 99.0),
            outgoing: HalfEdgeHandle::default(),
            name: None,
            tolerance: None,
        });
        let before = (brep.edges.len(), brep.faces.len(), brep.loops.len());

        let result = brep.merge_faces_across_edge_then_split(
            shared,
            &[FaceMergeSplitRequest {
                start: missing,
                end: brep.half_edges[shared].origin,
            }],
        );

        assert!(result.is_err());
        assert_eq!(
            before,
            (brep.edges.len(), brep.faces.len(), brep.loops.len())
        );
    }

    fn adjacent_face_pair(brep: &BRep) -> (FaceHandle, FaceHandle) {
        for (half_edge, data) in &brep.half_edges {
            let first = brep.loops[data.loop_].face;
            let second = brep.loops[brep.half_edges[data.twin].loop_].face;
            if first != second {
                return (first, second);
            }
            let _ = half_edge;
        }
        panic!("fixture has adjacent faces")
    }

    fn shared_half_edge(brep: &BRep, first: FaceHandle, second: FaceHandle) -> HalfEdgeHandle {
        for (half_edge, data) in &brep.half_edges {
            let source = brep.loops[data.loop_].face;
            let target = brep.loops[brep.half_edges[data.twin].loop_].face;
            if source == first && target == second {
                return half_edge;
            }
        }
        panic!("fixture has shared half-edge")
    }

    fn loop_vertices(brep: &BRep, loop_handle: LoopHandle) -> Vec<VertexHandle> {
        let start = brep.loops[loop_handle].half_edge;
        let mut cursor = start;
        let mut vertices = Vec::new();
        loop {
            vertices.push(brep.half_edges[cursor].origin);
            cursor = brep.half_edges[cursor].next;
            if cursor == start {
                return vertices;
            }
        }
    }
}
