//! make_acorn / kill_acorn: an isolated vertex in a region.

use super::EulerError;
use crate::geometry::Point3;
use crate::topology::{BRep, RegionHandle, ShellHandle, Vertex, VertexHandle};

/// What [`BRep::make_acorn`] made.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AcornResult {
    /// The vertex.
    pub vertex: VertexHandle,
    /// The shell of the region holding it.
    pub shell: ShellHandle,
}

impl BRep {
    /// Make an isolated vertex at `position` in `region`, held by a new
    /// shell of the region. Counts: vertices +1, shells +1, acorns +1.
    ///
    /// # Errors
    ///
    /// [`EulerError::UnknownRegion`].
    pub fn make_acorn(
        &mut self,
        region: RegionHandle,
        position: Point3,
    ) -> Result<AcornResult, EulerError> {
        if !self.regions.contains_key(region) {
            return Err(EulerError::UnknownRegion(region));
        }
        let vertex = self.vertices.insert(Vertex {
            position,
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let shell = self
            .add_acorn_shell(region, vertex)
            .ok_or(EulerError::UnknownRegion(region))?;
        Ok(AcornResult { vertex, shell })
    }

    /// Remove acorn `vertex` and its shell, the inverse of
    /// [`BRep::make_acorn`]. Counts: vertices −1, shells −1, acorns −1.
    ///
    /// # Errors
    ///
    /// [`EulerError::UnknownVertex`], or [`EulerError::NotAnAcorn`] when
    /// no shell holds the vertex as an acorn.
    pub fn kill_acorn(&mut self, vertex: VertexHandle) -> Result<(), EulerError> {
        if !self.vertices.contains_key(vertex) {
            return Err(EulerError::UnknownVertex(vertex));
        }
        let shell = self
            .acorn_shell(vertex)
            .ok_or(EulerError::NotAnAcorn { vertex })?;
        self.shells.remove(shell);
        self.vertices.remove(vertex);
        Ok(())
    }

    /// The shell holding `vertex` as an acorn.
    pub fn acorn_shell(&self, vertex: VertexHandle) -> Option<ShellHandle> {
        self.shells
            .iter()
            .find(|(_, s)| s.acorn == Some(vertex))
            .map(|(h, _)| h)
    }
}
