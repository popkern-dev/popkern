//! Retire whole faces to expose their opposite half-edges as boundary.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle};

/// Result of retiring production faces and exposing their opposite edges.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundaryFaceRetireResult {
    /// Faces removed from the topology.
    pub retired_faces: Vec<FaceHandle>,
    /// Edges whose surviving half-edge became a boundary side.
    pub exposed_edges: Vec<EdgeHandle>,
}

impl BRep {
    /// Retire whole faces and expose their opposite half-edges as boundary edges.
    ///
    /// This is a transaction boundary for boundary-preparation work. The live
    /// BRep is only replaced after the candidate passes structural validation.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if a face or loop handle is stale.
    /// - [`CoreError::PreconditionFailed`] if a face has inner loops, exposes
    ///   no surviving twin side, has inconsistent twin pairing, or would remove
    ///   the last face in its shell.
    pub fn retire_faces_to_boundary(
        &mut self,
        faces: &[FaceHandle],
    ) -> Result<BoundaryFaceRetireResult, CoreError> {
        let mut candidate = self.clone();
        let result = candidate.apply_retire_faces_to_boundary(faces)?;
        candidate.verify_retire_faces_to_boundary_postconditions(&result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_retire_faces_to_boundary(
        &mut self,
        faces: &[FaceHandle],
    ) -> Result<BoundaryFaceRetireResult, CoreError> {
        let mut retired_faces = Vec::new();
        let mut exposed_edges = Vec::new();
        for &face in faces {
            let retiring_loop = self.retirable_outer_loop(face)?;
            let doomed_half_edges = self.loop_half_edges_for_retirement(retiring_loop)?;
            for half_edge in &doomed_half_edges {
                let twin = self.half_edges[*half_edge].twin;
                let edge = self.half_edges[*half_edge].edge;
                let twin_data = self.half_edges.get(twin).ok_or_else(|| {
                    CoreError::PreconditionFailed(
                        "retiring face would expose an edge without a surviving twin".into(),
                    )
                })?;
                if twin_data.twin != *half_edge {
                    return Err(CoreError::PreconditionFailed(
                        "retiring face found inconsistent twin pairing".into(),
                    ));
                }
                self.half_edges[twin].twin = HalfEdgeHandle::default();
                self.edges[edge].half_edge = twin;
                self.boundary_edges.insert(edge);
                exposed_edges.push(edge);
            }
            self.repair_outgoing_handles_before_retirement(&doomed_half_edges);
            for half_edge in doomed_half_edges {
                self.half_edges.remove(half_edge);
            }
            self.repair_shell_face_after_retirement(face)?;
            self.loops.remove(retiring_loop);
            self.faces.remove(face);
            retired_faces.push(face);
        }
        Ok(BoundaryFaceRetireResult {
            retired_faces,
            exposed_edges,
        })
    }

    fn verify_retire_faces_to_boundary_postconditions(
        &self,
        result: &BoundaryFaceRetireResult,
    ) -> Result<(), CoreError> {
        for face in &result.retired_faces {
            if self.faces.contains_key(*face) {
                return Err(CoreError::PreconditionFailed(
                    "retired boundary face still exists after commit".into(),
                ));
            }
        }
        for edge in &result.exposed_edges {
            if !self.boundary_edges.contains(edge) {
                return Err(CoreError::PreconditionFailed(
                    "retired boundary face exposed edge is not boundary".into(),
                ));
            }
            let exposed_half_edge = self
                .edges
                .get(*edge)
                .ok_or_else(|| {
                    CoreError::PreconditionFailed(
                        "retired boundary face exposed edge no longer exists".into(),
                    )
                })?
                .half_edge;
            let exposed = self.half_edges.get(exposed_half_edge).ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "retired boundary face exposed half-edge no longer exists".into(),
                )
            })?;
            if exposed.twin != HalfEdgeHandle::default() {
                return Err(CoreError::PreconditionFailed(
                    "retired boundary face exposed half-edge kept a live twin".into(),
                ));
            }
            self.verify_retired_boundary_surviving_loop(exposed.loop_, exposed_half_edge)?;
        }
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)?;
        Ok(())
    }

    fn verify_retired_boundary_surviving_loop(
        &self,
        loop_: LoopHandle,
        required_half_edge: HalfEdgeHandle,
    ) -> Result<(), CoreError> {
        let loop_data = self.loops.get(loop_).ok_or_else(|| {
            CoreError::PreconditionFailed(
                "retired boundary face exposed half-edge references a retired loop".into(),
            )
        })?;
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut visited = BTreeSet::new();
        for _ in 0..=self.half_edges.len() {
            let half_edge = self.half_edges.get(cursor).ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "retired boundary surviving loop references a missing half-edge".into(),
                )
            })?;
            if half_edge.loop_ != loop_ {
                return Err(CoreError::PreconditionFailed(
                    "retired boundary surviving loop left its loop".into(),
                ));
            }
            if !visited.insert(cursor) {
                return Err(CoreError::PreconditionFailed(
                    "retired boundary surviving loop repeats before closing".into(),
                ));
            }
            cursor = half_edge.next;
            if cursor == start {
                if visited.contains(&required_half_edge) {
                    return Ok(());
                }
                return Err(CoreError::PreconditionFailed(
                    "retired boundary surviving loop does not contain exposed half-edge".into(),
                ));
            }
        }
        Err(CoreError::PreconditionFailed(
            "retired boundary surviving loop did not close".into(),
        ))
    }

    fn retirable_outer_loop(
        &self,
        face: FaceHandle,
    ) -> Result<crate::topology::LoopHandle, CoreError> {
        let face_data = self
            .faces
            .get(face)
            .ok_or_else(|| CoreError::InvalidHandle("retire face handle is invalid".into()))?;
        if !face_data.inner_loops.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "retire_faces_to_boundary does not retire faces with inner loops".into(),
            ));
        }
        Ok(face_data.outer_loop)
    }

    fn loop_half_edges_for_retirement(
        &self,
        loop_handle: crate::topology::LoopHandle,
    ) -> Result<Vec<HalfEdgeHandle>, CoreError> {
        let start = self
            .loops
            .get(loop_handle)
            .ok_or_else(|| CoreError::InvalidHandle("retire loop handle is invalid".into()))?
            .half_edge;
        let mut cursor = start;
        let mut half_edges = Vec::new();
        for _ in 0..=self.half_edges.len() {
            let half_edge = self.half_edges.get(cursor).ok_or_else(|| {
                CoreError::InvalidHandle("retire loop half-edge is invalid".into())
            })?;
            if half_edge.loop_ != loop_handle {
                return Err(CoreError::PreconditionFailed(
                    "retire loop walk left the requested loop".into(),
                ));
            }
            half_edges.push(cursor);
            cursor = half_edge.next;
            if cursor == start {
                return Ok(half_edges);
            }
        }
        Err(CoreError::PreconditionFailed(
            "retire loop walk did not close".into(),
        ))
    }

    fn repair_outgoing_handles_before_retirement(&mut self, doomed_half_edges: &[HalfEdgeHandle]) {
        for &half_edge in doomed_half_edges {
            let origin = self.half_edges[half_edge].origin;
            if self.vertices[origin].outgoing != half_edge {
                continue;
            }
            let replacement = self
                .half_edges
                .iter()
                .find(|(candidate, data)| {
                    data.origin == origin && !doomed_half_edges.contains(candidate)
                })
                .map(|(candidate, _)| candidate)
                .unwrap_or_default();
            self.vertices[origin].outgoing = replacement;
        }
    }

    fn repair_shell_face_after_retirement(&mut self, face: FaceHandle) -> Result<(), CoreError> {
        let shell = self
            .faces
            .get(face)
            .ok_or_else(|| CoreError::InvalidHandle("retire face handle is invalid".into()))?
            .shell;
        if self.shells[shell].face != face {
            return Ok(());
        }
        let replacement = self
            .faces
            .iter()
            .find(|(candidate, data)| *candidate != face && data.shell == shell)
            .map(|(candidate, _)| candidate)
            .ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "retire_faces_to_boundary cannot retire the last face in a shell".into(),
                )
            })?;
        self.shells[shell].face = replacement;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn retire_cube_face_to_boundary_commits_open_shell() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let face = brep.faces.keys().next().expect("cube has a face");
        let result = brep
            .retire_faces_to_boundary(&[face])
            .expect("face retires");

        assert_eq!(result.retired_faces, vec![face]);
        assert_eq!(result.exposed_edges.len(), 4);
        assert_eq!(brep.faces.len(), 5);
        assert_eq!(brep.loops.len(), 5);
        assert_eq!(brep.edges.len(), 12);
        assert_eq!(brep.boundary_edges.len(), 4);
        brep.validate().unwrap();
    }

    #[test]
    fn retire_last_face_rolls_back() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let faces = brep.faces.keys().collect::<Vec<_>>();
        let before = (brep.faces.len(), brep.loops.len(), brep.half_edges.len());
        let err = brep
            .retire_faces_to_boundary(&faces)
            .expect_err("last face in shell cannot retire");

        assert!(
            err.to_string().contains("retir")
                || err.to_string().contains("expose")
                || err.to_string().contains("invalid")
        );
        assert_eq!(
            before,
            (brep.faces.len(), brep.loops.len(), brep.half_edges.len())
        );
    }
}
