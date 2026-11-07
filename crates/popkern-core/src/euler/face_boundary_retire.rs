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
    /// Faces retired in one call may share edges: an edge both of whose
    /// faces go is removed with them, as is a vertex left on no
    /// half-edge, so the faces kept are a sheet with their own boundary.
    /// A solid shell opened this way keeps its region, as a solid being
    /// rebuilt does: a construction that closes it again leaves the
    /// solid whole. A caller whose result is the sheet settles it into
    /// the void ([`BRep::settle_open_shells`]); until then the opening
    /// fails the cell-structure checks of `validate_all`.
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
        // Faces retired together: an edge a face of this call exposed,
        // reached again by another, has lost both sides and goes, as
        // does a vertex left on no half-edge (a box keeping two faces
        // of six; an edge open before the call still declines).
        let mut exposed_here = std::collections::HashSet::new();
        let mut touched_vertices = Vec::new();
        for &face in faces {
            let retiring_loop = self.retirable_outer_loop(face)?;
            let doomed_half_edges = self.loop_half_edges_for_retirement(retiring_loop)?;
            let mut lost = Vec::new();
            for half_edge in &doomed_half_edges {
                let twin = self.half_edges[*half_edge].mate();
                let edge = self.half_edges[*half_edge].edge;
                touched_vertices.push(self.half_edges[*half_edge].origin);
                if twin.is_none() && exposed_here.contains(&edge) {
                    lost.push(edge);
                    continue;
                }
                let Some((twin, twin_data)) =
                    twin.and_then(|t| self.half_edges.get(t).map(|data| (t, data)))
                else {
                    return Err(CoreError::PreconditionFailed(
                        "retiring face would expose an edge without a surviving twin".into(),
                    ));
                };
                if twin_data.mate() != Some(*half_edge) {
                    return Err(CoreError::PreconditionFailed(
                        "retiring face found inconsistent twin pairing".into(),
                    ));
                }
                self.half_edges[twin].clear_mate();
                // The surviving twin becomes canonical; it traverses the
                // curve backward, so the curve turns to run its way (the
                // canonical-forward convention, PK-274, as face-region
                // removal keeps it). A closed edge has no vertex to tell
                // the two directions apart, so without this its pcurve
                // checks against the curve reversed (PK-978: an open
                // tube's top circle read 2R off).
                if self.edges[edge].half_edge == *half_edge
                    && let Some(curve) = self.edges[edge].curve.as_mut()
                {
                    std::mem::swap(&mut curve.t_start, &mut curve.t_end);
                }
                self.edges[edge].half_edge = twin;
                exposed_edges.push(edge);
                exposed_here.insert(edge);
            }
            self.repair_outgoing_handles_before_retirement(&doomed_half_edges);
            for half_edge in doomed_half_edges {
                self.remove_half_edge_reanchored(half_edge);
            }
            for edge in lost {
                self.edges.remove(edge);
                exposed_here.remove(&edge);
                exposed_edges.retain(|&e| e != edge);
            }
            self.repair_shell_face_after_retirement(face)?;
            self.loops.remove(retiring_loop);
            self.faces.remove(face);
            retired_faces.push(face);
        }
        let used: std::collections::HashSet<_> =
            self.half_edges.values().map(|h| h.origin).collect();
        for vertex in touched_vertices {
            if !used.contains(&vertex) {
                self.vertices.remove(vertex);
            }
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
            if !self.is_free_edge(*edge) {
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
            if exposed.mate().is_some() {
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
        self.faces
            .get(face)
            .ok_or_else(|| CoreError::InvalidHandle("retire face handle is invalid".into()))?;
        if !self.move_shell_anchor_off(face) {
            return Err(CoreError::PreconditionFailed(
                "retire_faces_to_boundary cannot retire the last face in a shell".into(),
            ));
        }
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
        assert_eq!(brep.free_edges().len(), 4);
        brep.validate().unwrap();
    }

    /// A retired cap leaves its circle with one half-edge: the
    /// survivor runs the circle backward, so the curve turns to run its
    /// way, or the closed edge's pcurve checks against it reversed
    /// (2R off, PK-978).
    #[test]
    fn retired_caps_leave_an_open_tube_that_validates_once_settled() {
        use crate::geometry::surface::SurfaceGeometry;
        let mut tube = BRep::cylinder(0.3, 2.0, Tolerance::default()).expect("cylinder");
        let caps: Vec<_> = tube
            .faces
            .iter()
            .filter(|(_, f)| {
                matches!(
                    f.surface.as_ref().map(|s| &s.geometry),
                    Some(SurfaceGeometry::Plane { .. })
                )
            })
            .map(|(h, _)| h)
            .collect();
        tube.retire_faces_to_boundary(&caps).expect("caps retire");
        assert_opening_reported_then_settles(&mut tube);
    }

    /// Four faces of a cube retired together: the edges between two
    /// retired faces go with them, as do the two corners on no kept
    /// face, and the two faces kept are a creased sheet of their own
    /// (7 edges, 6 open, 6 vertices).
    #[test]
    fn faces_retired_together_drop_the_edges_between_them() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let first = brep.faces.keys().next().expect("a face");
        let face_of = |b: &BRep, h| b.loops[b.half_edges[h].loop_].face;
        let he = brep.loops[brep.faces[first].outer_loop].half_edge;
        let second = face_of(&brep, brep.half_edges[he].mate().unwrap_or_default());
        let doomed: Vec<_> = brep
            .faces
            .keys()
            .filter(|&f| f != first && f != second)
            .collect();
        brep.retire_faces_to_boundary(&doomed)
            .expect("faces retire together");
        assert_eq!(
            (
                brep.faces.len(),
                brep.edges.len(),
                brep.free_edges().len(),
                brep.vertices.len()
            ),
            (2, 7, 6, 6)
        );
        assert_opening_reported_then_settles(&mut brep);
    }

    /// The opened solid keeps its region, so each free edge's wedge
    /// spans the solid's shell and the outside's (the opening is
    /// reported); settled into the void the sheet validates clean.
    fn assert_opening_reported_then_settles(brep: &mut BRep) {
        let opening = brep.validate_all();
        assert!(
            !opening.is_empty()
                && opening
                    .iter()
                    .all(|e| matches!(e, crate::error::ValidationError::WedgeShellMismatch { .. })),
            "{opening:?}"
        );
        assert_eq!(opening.len(), brep.free_edges().len());
        assert_eq!(brep.settle_open_shells().len(), 1);
        let findings = brep.validate_all();
        assert!(findings.is_empty(), "{findings:?}");
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
