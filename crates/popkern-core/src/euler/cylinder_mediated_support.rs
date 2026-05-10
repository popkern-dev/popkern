//! Preflight contract for cylinder-mediated support continuation.

use crate::error::CoreError;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

/// Request for a cylinder-mediated outside/outside support-continuation plan.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CylinderMediatedSupportContinuationRequest {
    /// Face that should own the support side after the continuation.
    pub support_face: FaceHandle,
    /// Existing outside/cap face containing the generator endpoints.
    pub outside_face: FaceHandle,
    /// Cylinder face mediating the outside-to-support continuation.
    pub cylinder_face: FaceHandle,
    /// Generator start vertex on the outside loop.
    pub start: VertexHandle,
    /// Generator end vertex on the outside loop.
    pub end: VertexHandle,
    /// Outside loop containing both generator endpoints.
    pub outside_loop: LoopHandle,
    /// Half-edge on `outside_loop` whose origin is `start`.
    pub outside_start_half_edge: HalfEdgeHandle,
    /// Half-edge on `outside_loop` whose origin is `end`.
    pub outside_end_half_edge: HalfEdgeHandle,
    /// Support loop adjacent to the cylinder.
    pub support_loop: LoopHandle,
    /// Support-side half-edge adjacent to the cylinder.
    pub support_cylinder_half_edge: HalfEdgeHandle,
    /// Outside-side half-edge adjacent to the cylinder.
    pub outside_cylinder_half_edge: HalfEdgeHandle,
    /// Directed cylinder interval from outside adjacency to support adjacency.
    pub cylinder_interval: Vec<HalfEdgeHandle>,
}

#[cfg(test)]
#[path = "cylinder_mediated_support_tests.rs"]
mod tests;

/// Preflighted plan for a cylinder-mediated support continuation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CylinderMediatedSupportContinuationPlan {
    /// Validated request operands.
    pub request: CylinderMediatedSupportContinuationRequest,
    /// Topology delta is intentionally unknown until the mutating operator is
    /// implemented and proven.
    pub expected_delta: CylinderMediatedSupportContinuationDelta,
}

/// Expected topology delta status for the continuation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CylinderMediatedSupportContinuationDelta {
    /// The mutating Euler sequence has not been accepted yet.
    Unknown {
        /// Stable diagnostic reason.
        reason: &'static str,
    },
}

impl BRep {
    /// Validate operands for a future cylinder-mediated support continuation.
    ///
    /// This preflight is read-only. It deliberately does not claim execution
    /// readiness because the Euler sequence and topology delta are not yet
    /// accepted.
    pub fn preflight_cylinder_mediated_support_continuation(
        &self,
        request: CylinderMediatedSupportContinuationRequest,
    ) -> Result<CylinderMediatedSupportContinuationPlan, CoreError> {
        self.require_face(request.support_face, "support face")?;
        self.require_face(request.outside_face, "outside face")?;
        self.require_face(request.cylinder_face, "cylinder face")?;
        self.require_loop_on_face(request.outside_loop, request.outside_face, "outside loop")?;
        self.require_loop_on_face(request.support_loop, request.support_face, "support loop")?;
        self.require_vertex(request.start, "start vertex")?;
        self.require_vertex(request.end, "end vertex")?;
        self.require_half_edge_origin_on_loop(
            request.outside_start_half_edge,
            request.outside_loop,
            request.start,
            "outside start half-edge",
        )?;
        self.require_half_edge_origin_on_loop(
            request.outside_end_half_edge,
            request.outside_loop,
            request.end,
            "outside end half-edge",
        )?;
        self.require_adjacency(
            request.support_cylinder_half_edge,
            request.support_loop,
            request.cylinder_face,
            "support/cylinder adjacency",
        )?;
        self.require_adjacency(
            request.outside_cylinder_half_edge,
            request.outside_loop,
            request.cylinder_face,
            "outside/cylinder adjacency",
        )?;
        self.require_directed_cylinder_interval(
            request.cylinder_face,
            request.outside_cylinder_half_edge,
            request.support_cylinder_half_edge,
            &request.cylinder_interval,
        )?;
        self.require_same_shell(
            request.support_face,
            request.outside_face,
            request.cylinder_face,
        )?;
        Ok(CylinderMediatedSupportContinuationPlan {
            request,
            expected_delta: CylinderMediatedSupportContinuationDelta::Unknown {
                reason: "cylinder-mediated-support-continuation-operator-required",
            },
        })
    }

    fn require_face(&self, face: FaceHandle, label: &str) -> Result<(), CoreError> {
        self.faces
            .contains_key(face)
            .then_some(())
            .ok_or_else(|| CoreError::InvalidHandle(format!("{label} {face:?} is not live")))
    }

    fn require_loop_on_face(
        &self,
        loop_handle: LoopHandle,
        face: FaceHandle,
        label: &str,
    ) -> Result<(), CoreError> {
        let data = self.loops.get(loop_handle).ok_or_else(|| {
            CoreError::InvalidHandle(format!("{label} {loop_handle:?} is not live"))
        })?;
        if data.face != face {
            return Err(CoreError::PreconditionFailed(format!(
                "{label} {loop_handle:?} belongs to {:?}, expected {face:?}",
                data.face
            )));
        }
        Ok(())
    }

    fn require_vertex(&self, vertex: VertexHandle, label: &str) -> Result<(), CoreError> {
        self.vertices
            .contains_key(vertex)
            .then_some(())
            .ok_or_else(|| CoreError::InvalidHandle(format!("{label} {vertex:?} is not live")))
    }

    fn require_half_edge_origin_on_loop(
        &self,
        half_edge: HalfEdgeHandle,
        loop_handle: LoopHandle,
        origin: VertexHandle,
        label: &str,
    ) -> Result<(), CoreError> {
        let data = self.half_edges.get(half_edge).ok_or_else(|| {
            CoreError::InvalidHandle(format!("{label} {half_edge:?} is not live"))
        })?;
        if data.loop_ != loop_handle {
            return Err(CoreError::PreconditionFailed(format!(
                "{label} {half_edge:?} belongs to {:?}, expected {loop_handle:?}",
                data.loop_
            )));
        }
        if data.origin != origin {
            return Err(CoreError::PreconditionFailed(format!(
                "{label} {half_edge:?} starts at {:?}, expected {origin:?}",
                data.origin
            )));
        }
        Ok(())
    }

    fn require_adjacency(
        &self,
        half_edge: HalfEdgeHandle,
        source_loop: LoopHandle,
        adjacent_face: FaceHandle,
        label: &str,
    ) -> Result<(), CoreError> {
        let data = self.half_edges.get(half_edge).ok_or_else(|| {
            CoreError::InvalidHandle(format!("{label} {half_edge:?} is not live"))
        })?;
        if data.loop_ != source_loop {
            return Err(CoreError::PreconditionFailed(format!(
                "{label} {half_edge:?} belongs to {:?}, expected {source_loop:?}",
                data.loop_
            )));
        }
        let twin = self.half_edges.get(data.twin).ok_or_else(|| {
            CoreError::PreconditionFailed(format!("{label} {half_edge:?} has no live twin"))
        })?;
        let twin_face = self.loops[twin.loop_].face;
        if twin_face != adjacent_face {
            return Err(CoreError::PreconditionFailed(format!(
                "{label} {half_edge:?} is adjacent to {twin_face:?}, expected {adjacent_face:?}"
            )));
        }
        Ok(())
    }

    fn require_directed_cylinder_interval(
        &self,
        cylinder_face: FaceHandle,
        outside_cylinder_half_edge: HalfEdgeHandle,
        support_cylinder_half_edge: HalfEdgeHandle,
        interval: &[HalfEdgeHandle],
    ) -> Result<(), CoreError> {
        if interval.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "cylinder-mediated support continuation requires a non-empty cylinder interval"
                    .into(),
            ));
        }
        let outside_cylinder_twin = self.half_edges[outside_cylinder_half_edge].twin;
        let support_cylinder_twin = self.half_edges[support_cylinder_half_edge].twin;
        if interval.first().copied() != Some(outside_cylinder_twin)
            || interval.last().copied() != Some(support_cylinder_twin)
        {
            return Err(CoreError::PreconditionFailed(format!(
                "cylinder interval endpoints mismatch: expected {outside_cylinder_twin:?}->{support_cylinder_twin:?}, got {:?}->{:?}",
                interval.first(),
                interval.last()
            )));
        }
        let cylinder_loop = self.half_edges[outside_cylinder_twin].loop_;
        if self.loops[cylinder_loop].face != cylinder_face {
            return Err(CoreError::PreconditionFailed(format!(
                "outside cylinder half-edge {outside_cylinder_twin:?} belongs to loop {cylinder_loop:?} on {:?}, expected {cylinder_face:?}",
                self.loops[cylinder_loop].face
            )));
        }
        for &half_edge in interval {
            let data = self.half_edges.get(half_edge).ok_or_else(|| {
                CoreError::InvalidHandle(format!(
                    "cylinder interval half-edge {half_edge:?} is not live"
                ))
            })?;
            if data.loop_ != cylinder_loop {
                return Err(CoreError::PreconditionFailed(format!(
                    "cylinder interval half-edge {half_edge:?} leaves loop {cylinder_loop:?}"
                )));
            }
        }
        for pair in interval.windows(2) {
            if self.half_edges[pair[0]].next != pair[1] {
                return Err(CoreError::PreconditionFailed(format!(
                    "cylinder interval is not directed at {:?}->{:?}",
                    pair[0], pair[1]
                )));
            }
        }
        Ok(())
    }

    fn require_same_shell(
        &self,
        support_face: FaceHandle,
        outside_face: FaceHandle,
        cylinder_face: FaceHandle,
    ) -> Result<(), CoreError> {
        let shell = self.faces[support_face].shell;
        if self.faces[outside_face].shell != shell || self.faces[cylinder_face].shell != shell {
            return Err(CoreError::PreconditionFailed(format!(
                "cylinder-mediated support continuation requires one shell: support={:?} outside={:?} cylinder={:?}",
                self.faces[support_face].shell,
                self.faces[outside_face].shell,
                self.faces[cylinder_face].shell
            )));
        }
        Ok(())
    }
}
