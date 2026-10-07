//! Shared contracts for transactional Euler operators.

use crate::error::CoreError;
use crate::topology::{BRep, HalfEdgeHandle, LoopHandle};

/// Snapshot of topology arena sizes used for Euler delta checks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct EulerTopologyCounts {
    pub(crate) vertices: usize,
    pub(crate) edges: usize,
    pub(crate) half_edges: usize,
    pub(crate) loops: usize,
    pub(crate) faces: usize,
    pub(crate) shells: usize,
    pub(crate) solids: usize,
    pub(crate) boundary_edges: usize,
}

/// Signed topology count delta expected from an Euler operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct EulerDelta {
    pub(crate) vertices: isize,
    pub(crate) edges: isize,
    pub(crate) half_edges: isize,
    pub(crate) loops: isize,
    pub(crate) faces: isize,
    pub(crate) shells: isize,
    pub(crate) solids: isize,
    pub(crate) boundary_edges: isize,
}

impl EulerTopologyCounts {
    pub(crate) fn from_brep(brep: &BRep) -> Self {
        Self {
            vertices: brep.vertices.len(),
            edges: brep.edges.len(),
            half_edges: brep.half_edges.len(),
            loops: brep.loops.len(),
            faces: brep.faces.len(),
            shells: brep.shells.len(),
            solids: brep.solids.len(),
            boundary_edges: brep.boundary_edges.len(),
        }
    }

    pub(crate) fn delta_to(self, after: Self) -> EulerDelta {
        EulerDelta {
            vertices: count_delta(self.vertices, after.vertices),
            edges: count_delta(self.edges, after.edges),
            half_edges: count_delta(self.half_edges, after.half_edges),
            loops: count_delta(self.loops, after.loops),
            faces: count_delta(self.faces, after.faces),
            shells: count_delta(self.shells, after.shells),
            solids: count_delta(self.solids, after.solids),
            boundary_edges: count_delta(self.boundary_edges, after.boundary_edges),
        }
    }
}

impl EulerDelta {
    pub(crate) const fn mvfs() -> Self {
        Self {
            vertices: 1,
            edges: 0,
            half_edges: 0,
            loops: 1,
            faces: 1,
            shells: 1,
            solids: 1,
            boundary_edges: 0,
        }
    }

    pub(crate) const fn kvfs() -> Self {
        Self {
            vertices: -1,
            edges: 0,
            half_edges: 0,
            loops: -1,
            faces: -1,
            shells: -1,
            solids: -1,
            boundary_edges: 0,
        }
    }

    pub(crate) const fn mev() -> Self {
        Self {
            vertices: 1,
            edges: 1,
            half_edges: 2,
            loops: 0,
            faces: 0,
            shells: 0,
            solids: 0,
            boundary_edges: 0,
        }
    }

    pub(crate) const fn kev() -> Self {
        Self {
            vertices: -1,
            edges: -1,
            half_edges: -2,
            loops: 0,
            faces: 0,
            shells: 0,
            solids: 0,
            boundary_edges: 0,
        }
    }

    pub(crate) const fn mef() -> Self {
        Self {
            vertices: 0,
            edges: 1,
            half_edges: 2,
            loops: 1,
            faces: 1,
            shells: 0,
            solids: 0,
            boundary_edges: 0,
        }
    }

    pub(crate) const fn kef() -> Self {
        Self {
            vertices: 0,
            edges: -1,
            half_edges: -2,
            loops: -1,
            faces: -1,
            shells: 0,
            solids: 0,
            boundary_edges: 0,
        }
    }

    pub(crate) const fn kemr() -> Self {
        Self {
            vertices: 0,
            edges: -1,
            half_edges: -2,
            loops: 1,
            faces: 0,
            shells: 0,
            solids: 0,
            boundary_edges: 0,
        }
    }

    pub(crate) const fn mekr() -> Self {
        Self {
            vertices: 0,
            edges: 1,
            half_edges: 2,
            loops: -1,
            faces: 0,
            shells: 0,
            solids: 0,
            boundary_edges: 0,
        }
    }

    pub(crate) const fn kfmrh() -> Self {
        Self {
            vertices: 0,
            edges: 0,
            half_edges: 0,
            loops: 0,
            faces: -1,
            shells: -1,
            solids: -1,
            boundary_edges: 0,
        }
    }

    pub(crate) const fn mfkrh() -> Self {
        Self {
            vertices: 0,
            edges: 0,
            half_edges: 0,
            loops: 0,
            faces: 1,
            shells: 1,
            solids: 1,
            boundary_edges: 0,
        }
    }

    /// This delta applied `times` times.
    pub(crate) const fn times(self, times: isize) -> Self {
        Self {
            vertices: self.vertices * times,
            edges: self.edges * times,
            half_edges: self.half_edges * times,
            loops: self.loops * times,
            faces: self.faces * times,
            shells: self.shells * times,
            solids: self.solids * times,
            boundary_edges: self.boundary_edges * times,
        }
    }

    /// This delta followed by `other`.
    pub(crate) const fn then(self, other: Self) -> Self {
        Self {
            vertices: self.vertices + other.vertices,
            edges: self.edges + other.edges,
            half_edges: self.half_edges + other.half_edges,
            loops: self.loops + other.loops,
            faces: self.faces + other.faces,
            shells: self.shells + other.shells,
            solids: self.solids + other.solids,
            boundary_edges: self.boundary_edges + other.boundary_edges,
        }
    }

    pub(crate) fn verify(
        op: &'static str,
        before: EulerTopologyCounts,
        after: EulerTopologyCounts,
        expected: Self,
    ) -> Result<(), CoreError> {
        let actual = before.delta_to(after);
        if actual != expected {
            return Err(CoreError::PreconditionFailed(format!(
                "{op} postcondition failed: expected_delta={expected:?} actual_delta={actual:?}"
            )));
        }
        Ok(())
    }
}

fn count_delta(before: usize, after: usize) -> isize {
    after as isize - before as isize
}

pub(crate) fn verify_loop_contains_half_edges(
    brep: &BRep,
    loop_: LoopHandle,
    required: &[HalfEdgeHandle],
) -> Result<(), CoreError> {
    let loop_data = brep
        .loops
        .get(loop_)
        .ok_or_else(|| CoreError::PreconditionFailed("loop is not live".into()))?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    let mut visited = std::collections::BTreeSet::new();
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor).ok_or_else(|| {
            CoreError::PreconditionFailed("loop references a missing half-edge".into())
        })?;
        if half_edge.loop_ != loop_ {
            return Err(CoreError::PreconditionFailed(
                "loop walk left the expected loop".into(),
            ));
        }
        if !visited.insert(cursor) {
            return Err(CoreError::PreconditionFailed(
                "loop repeats before closing".into(),
            ));
        }
        cursor = half_edge.next;
        if cursor == start {
            if required.iter().all(|half_edge| visited.contains(half_edge)) {
                return Ok(());
            }
            return Err(CoreError::PreconditionFailed(
                "loop is missing required half-edge".into(),
            ));
        }
    }
    Err(CoreError::PreconditionFailed(
        "loop did not close during verification".into(),
    ))
}
