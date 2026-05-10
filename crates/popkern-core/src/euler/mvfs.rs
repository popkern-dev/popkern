//! mvfs — Make Vertex, Face, Shell
//! kvfs — Kill Vertex, Face, Shell

use crate::error::CoreError;
use crate::euler::contract::{EulerDelta, EulerTopologyCounts};
use crate::geometry::Point3;
use crate::topology::{
    BRep, Face, FaceHandle, Loop, LoopHandle, Shell, ShellHandle, Solid, SolidHandle, Vertex,
    VertexHandle,
};

/// Result of the mvfs (Make Vertex, Face, Shell) operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MvfsResult {
    /// The newly created vertex.
    pub vertex: VertexHandle,
    /// The newly created face.
    pub face: FaceHandle,
    /// The newly created loop (the face's outer boundary).
    pub loop_: LoopHandle,
    /// The newly created shell.
    pub shell: ShellHandle,
    /// The newly created solid.
    pub solid: SolidHandle,
}

/// Explicit inverse request for a seed created by [`BRep::mvfs`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KvfsRequest {
    /// Isolated vertex to remove.
    pub vertex: VertexHandle,
    /// Isolated face to remove.
    pub face: FaceHandle,
    /// Empty outer loop to remove.
    pub loop_: LoopHandle,
    /// Single-face shell to remove.
    pub shell: ShellHandle,
    /// Solid owned by the shell.
    pub solid: SolidHandle,
}

/// Result of the kvfs (Kill Vertex, Face, Shell) operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KvfsResult {
    /// Removed vertex.
    pub vertex: VertexHandle,
    /// Removed face.
    pub face: FaceHandle,
    /// Removed loop.
    pub loop_: LoopHandle,
    /// Removed shell.
    pub shell: ShellHandle,
    /// Removed solid.
    pub solid: SolidHandle,
}

impl From<MvfsResult> for KvfsRequest {
    fn from(value: MvfsResult) -> Self {
        Self {
            vertex: value.vertex,
            face: value.face,
            loop_: value.loop_,
            shell: value.shell,
            solid: value.solid,
        }
    }
}

impl BRep {
    /// Make Vertex, Face, Shell — the seed operation for a new solid.
    ///
    /// Creates one vertex, face, loop, shell, and solid, all wired together.
    /// The loop and vertex have no half-edges yet.
    pub fn mvfs(&mut self, position: Point3) -> Result<MvfsResult, CoreError> {
        let before = EulerTopologyCounts::from_brep(self);
        let mut candidate = self.clone();
        let result = candidate.apply_mvfs(position);
        candidate.verify_mvfs_postconditions(before, result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_mvfs(&mut self, position: Point3) -> MvfsResult {
        let vertex = self.vertices.insert(Vertex {
            position,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });

        // We need face, loop, shell, solid — they are mutually referential.
        // Insert with defaults and patch.
        let face = self.faces.insert(Face {
            outer_loop: Default::default(),
            inner_loops: vec![],
            shell: Default::default(),
            surface: None,
            name: None,
            tolerance: None,
            role: Default::default(),
        });

        let loop_ = self.loops.insert(Loop {
            half_edge: Default::default(),
            face,
        });

        let shell = self.shells.insert(Shell {
            face,
            solid: Default::default(),
        });

        let solid = self.solids.insert(Solid { shell });

        // Patch the face with the real loop and shell handles.
        let face_data = &mut self.faces[face];
        face_data.outer_loop = loop_;
        face_data.shell = shell;

        // Patch the shell with the real solid handle.
        self.shells[shell].solid = solid;

        MvfsResult {
            vertex,
            face,
            loop_,
            shell,
            solid,
        }
    }

    /// Kill Vertex, Face, Shell — inverse of mvfs.
    ///
    /// Removes the solid that contains the given vertex, along with its shell,
    /// face, loop, and the vertex itself.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if the vertex handle is invalid.
    /// - [`CoreError::PreconditionFailed`] if the vertex has outgoing edges,
    ///   the shell has more than one face, the face has inner loops, or edges
    ///   still exist in the model for this solid.
    pub fn kvfs(&mut self, vertex: VertexHandle) -> Result<KvfsResult, CoreError> {
        if !self.vertices.contains_key(vertex) {
            return Err(CoreError::InvalidHandle("vertex handle is invalid".into()));
        }
        let mut candidates = self.mvfs_seed_candidates();
        if candidates.len() != 1 {
            return Err(CoreError::PreconditionFailed(format!(
                "kvfs(vertex) requires exactly one mvfs seed candidate, found {}; use kvfs_seed \
                 with explicit face/loop/shell/solid handles",
                candidates.len()
            )));
        }
        let mut request = candidates.remove(0);
        request.vertex = vertex;
        self.kvfs_seed(request)
    }

    /// Strict transactional inverse of [`BRep::mvfs`].
    pub fn kvfs_seed(&mut self, request: KvfsRequest) -> Result<KvfsResult, CoreError> {
        self.verify_kvfs_preconditions(request)?;
        let before = EulerTopologyCounts::from_brep(self);
        let mut candidate = self.clone();
        let result = candidate.apply_kvfs_seed(request);
        candidate.verify_kvfs_postconditions(before, result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_kvfs_seed(&mut self, request: KvfsRequest) -> KvfsResult {
        self.solids.remove(request.solid);
        self.shells.remove(request.shell);
        self.faces.remove(request.face);
        self.loops.remove(request.loop_);
        self.vertices.remove(request.vertex);

        KvfsResult {
            vertex: request.vertex,
            face: request.face,
            loop_: request.loop_,
            shell: request.shell,
            solid: request.solid,
        }
    }

    fn verify_mvfs_postconditions(
        &self,
        before: EulerTopologyCounts,
        result: MvfsResult,
    ) -> Result<(), CoreError> {
        EulerDelta::verify(
            "mvfs",
            before,
            EulerTopologyCounts::from_brep(self),
            EulerDelta::mvfs(),
        )?;
        let vertex = self.vertices.get(result.vertex).ok_or_else(|| {
            CoreError::PreconditionFailed("mvfs created vertex is not live".into())
        })?;
        if self.half_edges.contains_key(vertex.outgoing) {
            return Err(CoreError::PreconditionFailed(
                "mvfs created vertex has a live outgoing half-edge".into(),
            ));
        }
        let face = self
            .faces
            .get(result.face)
            .ok_or_else(|| CoreError::PreconditionFailed("mvfs created face is not live".into()))?;
        if face.outer_loop != result.loop_ || face.shell != result.shell {
            return Err(CoreError::PreconditionFailed(
                "mvfs face does not own the created loop and shell".into(),
            ));
        }
        if !face.inner_loops.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "mvfs created face unexpectedly has inner loops".into(),
            ));
        }
        let loop_data = self
            .loops
            .get(result.loop_)
            .ok_or_else(|| CoreError::PreconditionFailed("mvfs created loop is not live".into()))?;
        if loop_data.face != result.face || self.half_edges.contains_key(loop_data.half_edge) {
            return Err(CoreError::PreconditionFailed(
                "mvfs loop does not point to the created empty face".into(),
            ));
        }
        if self
            .shells
            .get(result.shell)
            .is_none_or(|shell| shell.face != result.face || shell.solid != result.solid)
        {
            return Err(CoreError::PreconditionFailed(
                "mvfs shell does not point to the created face and solid".into(),
            ));
        }
        if self
            .solids
            .get(result.solid)
            .is_none_or(|solid| solid.shell != result.shell)
        {
            return Err(CoreError::PreconditionFailed(
                "mvfs solid does not point to the created shell".into(),
            ));
        }
        self.validate_euler_poincare().map_err(CoreError::from)?;
        self.validate_vertex_consistency().map_err(CoreError::from)
    }

    fn verify_kvfs_preconditions(&self, request: KvfsRequest) -> Result<(), CoreError> {
        let vertex = self
            .vertices
            .get(request.vertex)
            .ok_or_else(|| CoreError::InvalidHandle("kvfs vertex handle is invalid".into()))?;
        if self.half_edges.contains_key(vertex.outgoing) {
            return Err(CoreError::PreconditionFailed(
                "kvfs vertex still has outgoing half-edges".into(),
            ));
        }
        let face = self
            .faces
            .get(request.face)
            .ok_or_else(|| CoreError::InvalidHandle("kvfs face handle is invalid".into()))?;
        if face.outer_loop != request.loop_
            || face.shell != request.shell
            || !face.inner_loops.is_empty()
        {
            return Err(CoreError::PreconditionFailed(
                "kvfs face is not an isolated mvfs face".into(),
            ));
        }
        let loop_data = self
            .loops
            .get(request.loop_)
            .ok_or_else(|| CoreError::InvalidHandle("kvfs loop handle is invalid".into()))?;
        if loop_data.face != request.face || self.half_edges.contains_key(loop_data.half_edge) {
            return Err(CoreError::PreconditionFailed(
                "kvfs loop is not the empty outer loop of the requested face".into(),
            ));
        }
        let shell = self
            .shells
            .get(request.shell)
            .ok_or_else(|| CoreError::InvalidHandle("kvfs shell handle is invalid".into()))?;
        if shell.face != request.face || shell.solid != request.solid {
            return Err(CoreError::PreconditionFailed(
                "kvfs shell does not own the requested face and solid".into(),
            ));
        }
        if self
            .solids
            .get(request.solid)
            .is_none_or(|solid| solid.shell != request.shell)
        {
            return Err(CoreError::InvalidHandle(
                "kvfs solid handle is invalid for the requested shell".into(),
            ));
        }
        let shell_face_count = self
            .faces
            .values()
            .filter(|candidate| candidate.shell == request.shell)
            .count();
        if shell_face_count != 1 {
            return Err(CoreError::PreconditionFailed(
                "kvfs shell must contain exactly one face".into(),
            ));
        }
        Ok(())
    }

    fn verify_kvfs_postconditions(
        &self,
        before: EulerTopologyCounts,
        result: KvfsResult,
    ) -> Result<(), CoreError> {
        EulerDelta::verify(
            "kvfs",
            before,
            EulerTopologyCounts::from_brep(self),
            EulerDelta::kvfs(),
        )?;
        if self.vertices.contains_key(result.vertex)
            || self.faces.contains_key(result.face)
            || self.loops.contains_key(result.loop_)
            || self.shells.contains_key(result.shell)
            || self.solids.contains_key(result.solid)
        {
            return Err(CoreError::PreconditionFailed(
                "kvfs left removed seed topology live".into(),
            ));
        }
        for face in self.faces.values() {
            if face.outer_loop == result.loop_
                || face.inner_loops.contains(&result.loop_)
                || face.shell == result.shell
            {
                return Err(CoreError::PreconditionFailed(
                    "kvfs left a live face referencing removed topology".into(),
                ));
            }
        }
        for loop_data in self.loops.values() {
            if loop_data.face == result.face {
                return Err(CoreError::PreconditionFailed(
                    "kvfs left a live loop referencing the removed face".into(),
                ));
            }
        }
        for shell in self.shells.values() {
            if shell.face == result.face || shell.solid == result.solid {
                return Err(CoreError::PreconditionFailed(
                    "kvfs left a live shell referencing removed topology".into(),
                ));
            }
        }
        for solid in self.solids.values() {
            if solid.shell == result.shell {
                return Err(CoreError::PreconditionFailed(
                    "kvfs left a live solid referencing the removed shell".into(),
                ));
            }
        }
        self.validate_euler_poincare().map_err(CoreError::from)?;
        self.validate_vertex_consistency().map_err(CoreError::from)
    }

    fn mvfs_seed_candidates(&self) -> Vec<KvfsRequest> {
        self.faces
            .iter()
            .filter_map(|(face_handle, face)| {
                if !face.inner_loops.is_empty() {
                    return None;
                }
                let loop_data = self.loops.get(face.outer_loop)?;
                if loop_data.face != face_handle
                    || self.half_edges.contains_key(loop_data.half_edge)
                {
                    return None;
                }
                let shell = self.shells.get(face.shell)?;
                if shell.face != face_handle {
                    return None;
                }
                let solid = self.solids.get(shell.solid)?;
                if solid.shell != face.shell {
                    return None;
                }
                let shell_face_count = self
                    .faces
                    .values()
                    .filter(|candidate| candidate.shell == face.shell)
                    .count();
                (shell_face_count == 1).then_some(KvfsRequest {
                    vertex: Default::default(),
                    face: face_handle,
                    loop_: face.outer_loop,
                    shell: face.shell,
                    solid: shell.solid,
                })
            })
            .collect()
    }
}

#[cfg(test)]
#[path = "mvfs_tests.rs"]
mod tests;
