//! Retire construction-residual faces after replacement topology is validated.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, TopologyRole};

/// Result of removing construction-residual topology.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstructionResidualRetireResult {
    /// Residual faces removed from the topology.
    pub retired_faces: Vec<FaceHandle>,
    /// Residual loops removed with those faces.
    pub removed_loops: Vec<LoopHandle>,
    /// Half-edges removed with those faces.
    pub removed_half_edges: Vec<HalfEdgeHandle>,
    /// Edges removed because every incident half-edge was residual.
    pub removed_edges: Vec<EdgeHandle>,
    /// Edges whose surviving side became a production boundary.
    pub exposed_edges: Vec<EdgeHandle>,
}

/// Result of removing unwired construction-residual edges.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnwiredConstructionResidualRetireResult {
    /// Residual edges removed from the topology.
    pub removed_edges: Vec<EdgeHandle>,
    /// Residual half-edges removed with those edges.
    pub removed_half_edges: Vec<HalfEdgeHandle>,
}

impl BRep {
    /// Remove faces already marked [`TopologyRole::ConstructionResidual`].
    ///
    /// The operator is transactional. It deletes residual-only edge records and
    /// exposes surviving production twins as boundary edges, then verifies that
    /// no live topology references the removed handles before commit.
    pub fn retire_construction_residual_faces(
        &mut self,
        faces: &[FaceHandle],
    ) -> Result<ConstructionResidualRetireResult, CoreError> {
        let mut candidate = self.clone();
        let result = candidate.apply_retire_construction_residual_faces(faces)?;
        candidate.verify_construction_residual_retire_postconditions(&result)?;
        *self = candidate;
        Ok(result)
    }

    /// Remove construction-residual edges whose half-edges are not wired into a live loop.
    ///
    /// A live half-edge without an owning loop is invalid for strict Euler
    /// operators. This operator removes only edge records explicitly marked as
    /// construction residue and only when none of their incident half-edges is
    /// owned by a live loop.
    pub fn retire_unwired_construction_residual_edges(
        &mut self,
    ) -> Result<UnwiredConstructionResidualRetireResult, CoreError> {
        let mut candidate = self.clone();
        let result = candidate.apply_retire_unwired_construction_residual_edges()?;
        candidate.verify_unwired_construction_residual_retire_postconditions(&result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_retire_construction_residual_faces(
        &mut self,
        faces: &[FaceHandle],
    ) -> Result<ConstructionResidualRetireResult, CoreError> {
        let mut retired_faces = Vec::new();
        let mut removed_loops = Vec::new();
        let mut removed_half_edges = Vec::new();
        let mut removed_edges = BTreeSet::new();
        let mut exposed_edges = BTreeSet::new();

        let mut doomed_half_edges = BTreeSet::new();
        for &face in faces {
            let loop_handle = self.retirable_residual_outer_loop(face)?;
            for half_edge in self.loop_half_edges_for_residual_retirement(loop_handle)? {
                doomed_half_edges.insert(half_edge);
            }
        }

        for &half_edge in &doomed_half_edges {
            let half_edge_data = self.half_edges.get(half_edge).ok_or_else(|| {
                CoreError::InvalidHandle("residual retire half-edge handle is invalid".into())
            })?;
            let edge = half_edge_data.edge;
            let twin = half_edge_data.twin;
            if doomed_half_edges.contains(&twin) {
                removed_edges.insert(edge);
                continue;
            }
            if let Some(twin_data) = self.half_edges.get(twin) {
                if twin_data.edge != edge || twin_data.twin != half_edge {
                    return Err(CoreError::PreconditionFailed(
                        "residual retire found inconsistent surviving twin pairing".into(),
                    ));
                }
                self.half_edges[twin].twin = HalfEdgeHandle::default();
                self.edges[edge].half_edge = twin;
                self.boundary_edges.insert(edge);
                exposed_edges.insert(edge);
            } else {
                removed_edges.insert(edge);
            }
        }

        self.repair_outgoing_handles_before_residual_retirement(&doomed_half_edges);
        for half_edge in &doomed_half_edges {
            self.remove_half_edge_reanchored(*half_edge);
            removed_half_edges.push(*half_edge);
        }
        for edge in &removed_edges {
            self.boundary_edges.remove(edge);
            self.edges.remove(*edge);
        }
        for &face in faces {
            let loop_handle = self
                .faces
                .get(face)
                .ok_or_else(|| {
                    CoreError::InvalidHandle("residual retire face handle is invalid".into())
                })?
                .outer_loop;
            self.repair_shell_face_after_residual_retirement(face)?;
            self.loops.remove(loop_handle);
            self.faces.remove(face);
            retired_faces.push(face);
            removed_loops.push(loop_handle);
        }

        Ok(ConstructionResidualRetireResult {
            retired_faces,
            removed_loops,
            removed_half_edges,
            removed_edges: removed_edges.into_iter().collect(),
            exposed_edges: exposed_edges.into_iter().collect(),
        })
    }

    fn retirable_residual_outer_loop(&self, face: FaceHandle) -> Result<LoopHandle, CoreError> {
        let face_data = self.faces.get(face).ok_or_else(|| {
            CoreError::InvalidHandle("residual retire face handle is invalid".into())
        })?;
        if face_data.role != TopologyRole::ConstructionResidual {
            return Err(CoreError::PreconditionFailed(
                "residual retire can only remove construction-residual faces".into(),
            ));
        }
        if !face_data.inner_loops.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "residual retire does not remove faces with inner loops".into(),
            ));
        }
        Ok(face_data.outer_loop)
    }

    fn loop_half_edges_for_residual_retirement(
        &self,
        loop_handle: LoopHandle,
    ) -> Result<Vec<HalfEdgeHandle>, CoreError> {
        let start = self
            .loops
            .get(loop_handle)
            .ok_or_else(|| CoreError::InvalidHandle("residual retire loop is invalid".into()))?
            .half_edge;
        let mut cursor = start;
        let mut half_edges = Vec::new();
        for _ in 0..=self.half_edges.len() {
            let half_edge = self.half_edges.get(cursor).ok_or_else(|| {
                CoreError::InvalidHandle("residual retire loop half-edge is invalid".into())
            })?;
            if half_edge.loop_ != loop_handle {
                return Err(CoreError::PreconditionFailed(
                    "residual retire loop walk left the requested loop".into(),
                ));
            }
            half_edges.push(cursor);
            cursor = half_edge.next;
            if cursor == start {
                for (half_edge, data) in &self.half_edges {
                    if data.loop_ == loop_handle && !half_edges.contains(&half_edge) {
                        half_edges.push(half_edge);
                    }
                }
                return Ok(half_edges);
            }
        }
        Err(CoreError::PreconditionFailed(
            "residual retire loop walk did not close".into(),
        ))
    }

    fn repair_outgoing_handles_before_residual_retirement(
        &mut self,
        doomed_half_edges: &BTreeSet<HalfEdgeHandle>,
    ) {
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

    fn repair_shell_face_after_residual_retirement(
        &mut self,
        face: FaceHandle,
    ) -> Result<(), CoreError> {
        let shell = self
            .faces
            .get(face)
            .ok_or_else(|| {
                CoreError::InvalidHandle("residual retire face handle is invalid".into())
            })?
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
                    "residual retire cannot remove the last face in a shell".into(),
                )
            })?;
        self.shells[shell].face = replacement;
        Ok(())
    }

    fn verify_construction_residual_retire_postconditions(
        &self,
        result: &ConstructionResidualRetireResult,
    ) -> Result<(), CoreError> {
        for face in &result.retired_faces {
            if self.faces.contains_key(*face) {
                return Err(CoreError::PreconditionFailed(
                    "residual retire left a retired face live".into(),
                ));
            }
        }
        for loop_handle in &result.removed_loops {
            if self.loops.contains_key(*loop_handle) {
                return Err(CoreError::PreconditionFailed(
                    "residual retire left a retired loop live".into(),
                ));
            }
        }
        for half_edge in &result.removed_half_edges {
            if self.half_edges.contains_key(*half_edge) {
                return Err(CoreError::PreconditionFailed(
                    "residual retire left a retired half-edge live".into(),
                ));
            }
        }
        for edge in &result.removed_edges {
            if self.edges.contains_key(*edge) || self.boundary_edges.contains(edge) {
                return Err(CoreError::PreconditionFailed(
                    "residual retire left a retired edge live".into(),
                ));
            }
        }
        for edge in &result.exposed_edges {
            let half_edge = self
                .edges
                .get(*edge)
                .ok_or_else(|| {
                    CoreError::PreconditionFailed("residual retire exposed stale edge".into())
                })?
                .half_edge;
            let half_edge_data = self.half_edges.get(half_edge).ok_or_else(|| {
                CoreError::PreconditionFailed("residual retire exposed stale half-edge".into())
            })?;
            if half_edge_data.twin != HalfEdgeHandle::default()
                || !self.boundary_edges.contains(edge)
            {
                return Err(CoreError::PreconditionFailed(
                    "residual retire exposed edge is not a boundary side".into(),
                ));
            }
        }
        self.verify_no_live_references_to_removed_residuals(result)?;
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        self.validate_vertex_consistency()
            .map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)?;
        Ok(())
    }

    fn verify_no_live_references_to_removed_residuals(
        &self,
        result: &ConstructionResidualRetireResult,
    ) -> Result<(), CoreError> {
        for face in self.faces.values() {
            if result.removed_loops.contains(&face.outer_loop)
                || face
                    .inner_loops
                    .iter()
                    .any(|loop_handle| result.removed_loops.contains(loop_handle))
            {
                return Err(CoreError::PreconditionFailed(
                    "live face references a retired residual loop".into(),
                ));
            }
        }
        for loop_data in self.loops.values() {
            if result.removed_half_edges.contains(&loop_data.half_edge)
                || result.retired_faces.contains(&loop_data.face)
            {
                return Err(CoreError::PreconditionFailed(
                    "live loop references retired residual topology".into(),
                ));
            }
        }
        for edge in self.edges.values() {
            if result.removed_half_edges.contains(&edge.half_edge) {
                return Err(CoreError::PreconditionFailed(
                    "live edge references a retired residual half-edge".into(),
                ));
            }
        }
        for (live_half_edge, half_edge) in &self.half_edges {
            if result.removed_edges.contains(&half_edge.edge) {
                return Err(CoreError::PreconditionFailed(format!(
                    "live half-edge {live_half_edge:?} references retired residual edge {:?}",
                    half_edge.edge
                )));
            }
            if result.removed_loops.contains(&half_edge.loop_) {
                return Err(CoreError::PreconditionFailed(format!(
                    "live half-edge {live_half_edge:?} references retired residual loop {:?}",
                    half_edge.loop_
                )));
            }
            if result.removed_half_edges.contains(&half_edge.next) {
                return Err(CoreError::PreconditionFailed(format!(
                    "live half-edge {live_half_edge:?} next references retired residual half-edge {:?}",
                    half_edge.next
                )));
            }
            if result.removed_half_edges.contains(&half_edge.prev) {
                return Err(CoreError::PreconditionFailed(format!(
                    "live half-edge {live_half_edge:?} prev references retired residual \
                         half-edge {:?}",
                    half_edge.prev
                )));
            }
        }
        Ok(())
    }

    fn apply_retire_unwired_construction_residual_edges(
        &mut self,
    ) -> Result<UnwiredConstructionResidualRetireResult, CoreError> {
        let candidates = self
            .edges
            .iter()
            .filter_map(|(edge, data)| {
                (data.role == TopologyRole::ConstructionResidual
                    && self.edge_is_unwired_construction_residual(edge))
                .then_some(edge)
            })
            .collect::<Vec<_>>();
        let mut removed_edges = Vec::new();
        let mut removed_half_edges = BTreeSet::new();
        for edge in candidates {
            for half_edge in self
                .half_edges
                .iter()
                .filter_map(|(half_edge, data)| (data.edge == edge).then_some(half_edge))
                .collect::<Vec<_>>()
            {
                removed_half_edges.insert(half_edge);
            }
            self.boundary_edges.remove(&edge);
            self.edges.remove(edge);
            removed_edges.push(edge);
        }
        self.repair_outgoing_handles_before_residual_retirement(&removed_half_edges);
        for half_edge in &removed_half_edges {
            self.remove_half_edge_reanchored(*half_edge);
        }
        Ok(UnwiredConstructionResidualRetireResult {
            removed_edges,
            removed_half_edges: removed_half_edges.into_iter().collect(),
        })
    }

    fn edge_is_unwired_construction_residual(&self, edge: EdgeHandle) -> bool {
        let mut has_half_edge = false;
        for (_, half_edge) in self
            .half_edges
            .iter()
            .filter(|(_, half_edge)| half_edge.edge == edge)
        {
            has_half_edge = true;
            if self.loops.contains_key(half_edge.loop_) {
                return false;
            }
        }
        has_half_edge
    }

    fn verify_unwired_construction_residual_retire_postconditions(
        &self,
        result: &UnwiredConstructionResidualRetireResult,
    ) -> Result<(), CoreError> {
        for edge in &result.removed_edges {
            if self.edges.contains_key(*edge) || self.boundary_edges.contains(edge) {
                return Err(CoreError::PreconditionFailed(
                    "unwired residual retire left a removed edge live".into(),
                ));
            }
        }
        for half_edge in &result.removed_half_edges {
            if self.half_edges.contains_key(*half_edge) {
                return Err(CoreError::PreconditionFailed(
                    "unwired residual retire left a removed half-edge live".into(),
                ));
            }
        }
        for edge in self.edges.values() {
            if result.removed_half_edges.contains(&edge.half_edge) {
                return Err(CoreError::PreconditionFailed(
                    "live edge references an unwired retired residual half-edge".into(),
                ));
            }
        }
        for vertex in self.vertices.values() {
            if result.removed_half_edges.contains(&vertex.outgoing) {
                return Err(CoreError::PreconditionFailed(
                    "live vertex references an unwired retired residual half-edge".into(),
                ));
            }
        }
        for loop_data in self.loops.values() {
            if result.removed_half_edges.contains(&loop_data.half_edge) {
                return Err(CoreError::PreconditionFailed(
                    "live loop references an unwired retired residual half-edge".into(),
                ));
            }
        }
        for (live_half_edge, half_edge) in &self.half_edges {
            if result.removed_edges.contains(&half_edge.edge) {
                return Err(CoreError::PreconditionFailed(format!(
                    "live half-edge {live_half_edge:?} references an unwired retired residual edge"
                )));
            }
            if result.removed_half_edges.contains(&half_edge.next)
                || result.removed_half_edges.contains(&half_edge.prev)
                || result.removed_half_edges.contains(&half_edge.twin)
            {
                return Err(CoreError::PreconditionFailed(format!(
                    "live half-edge {live_half_edge:?} references an unwired retired residual half-edge"
                )));
            }
        }
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        self.validate_vertex_consistency()
            .map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)?;
        Ok(())
    }
}

#[cfg(test)]
#[path = "construction_residual_retire_tests.rs"]
mod tests;
