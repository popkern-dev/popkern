//! Expose a same-loop edge pair as two boundary sides.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::topology::{BRep, Edge, EdgeHandle, HalfEdgeHandle};

/// Result of converting a same-loop edge into two boundary edge records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SameLoopBoundaryExposureResult {
    /// Original edge record kept by the requested half-edge.
    pub first_boundary_edge: EdgeHandle,
    /// New edge record assigned to the former twin half-edge.
    pub second_boundary_edge: EdgeHandle,
    /// Requested half-edge, retained in its original loop.
    pub first_half_edge: HalfEdgeHandle,
    /// Former twin half-edge, retained in its original loop.
    pub second_half_edge: HalfEdgeHandle,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SameLoopBoundaryCounts {
    vertices: i64,
    edges: i64,
    faces: i64,
    loops: i64,
}

impl BRep {
    /// Convert a same-loop twin pair into two boundary edge records.
    ///
    /// This is for bridge exposure cases where an existing cap/support loop
    /// contains both sides of one edge and both half-edges are protected bridge
    /// operands. Unlike `kemr` or `dissolve_same_loop_bridge`, this operator
    /// keeps both half-edges in the operand loop. It only splits the shared edge
    /// identity so each half-edge can act as an unowned boundary side for later
    /// bridge construction.
    pub fn expose_same_loop_edge_as_boundary(
        &mut self,
        half_edge: HalfEdgeHandle,
    ) -> Result<SameLoopBoundaryExposureResult, CoreError> {
        let before = SameLoopBoundaryCounts::capture(self);
        let mut candidate = self.clone();
        let result = candidate.apply_expose_same_loop_edge_as_boundary(half_edge)?;
        candidate.verify_same_loop_boundary_postconditions(before, result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_expose_same_loop_edge_as_boundary(
        &mut self,
        half_edge: HalfEdgeHandle,
    ) -> Result<SameLoopBoundaryExposureResult, CoreError> {
        let data = self
            .half_edges
            .get(half_edge)
            .ok_or_else(|| CoreError::InvalidHandle("same-loop half-edge is invalid".into()))?;
        let twin = data.twin;
        let edge = data.edge;
        if self.boundary_edges.contains(&edge) {
            return Err(CoreError::PreconditionFailed(
                "same-loop edge is already boundary".into(),
            ));
        }
        let twin_data = self
            .half_edges
            .get(twin)
            .ok_or_else(|| CoreError::InvalidHandle("same-loop twin is invalid".into()))?;
        if twin_data.twin != half_edge || twin_data.edge != edge {
            return Err(CoreError::PreconditionFailed(
                "same-loop edge has inconsistent twin pairing".into(),
            ));
        }
        if data.loop_ != twin_data.loop_ {
            return Err(CoreError::PreconditionFailed(
                "same-loop boundary exposure requires both half-edges on one loop".into(),
            ));
        }
        let loop_ = data.loop_;
        if !self.loops.contains_key(loop_) {
            return Err(CoreError::PreconditionFailed(
                "same-loop boundary exposure requires a live operand loop".into(),
            ));
        }

        let original_edge = self
            .edges
            .get(edge)
            .ok_or_else(|| CoreError::InvalidHandle("same-loop edge is invalid".into()))?
            .clone();
        let second_edge = self.edges.insert(Edge {
            half_edge: twin,
            curve: original_edge.curve,
            name: original_edge.name,
            tolerance: original_edge.tolerance,
            canonical: original_edge.canonical,
        });
        self.edges[edge].half_edge = half_edge;
        self.half_edges[half_edge].twin = HalfEdgeHandle::default();
        self.half_edges[twin].twin = HalfEdgeHandle::default();
        self.half_edges[twin].edge = second_edge;
        self.boundary_edges.insert(edge);
        self.boundary_edges.insert(second_edge);
        self.loops[loop_].half_edge = half_edge;
        Ok(SameLoopBoundaryExposureResult {
            first_boundary_edge: edge,
            second_boundary_edge: second_edge,
            first_half_edge: half_edge,
            second_half_edge: twin,
        })
    }

    fn verify_same_loop_boundary_postconditions(
        &self,
        before: SameLoopBoundaryCounts,
        result: SameLoopBoundaryExposureResult,
    ) -> Result<(), CoreError> {
        let actual_delta = SameLoopBoundaryCounts::capture(self).delta_from(before);
        let expected_delta = SameLoopBoundaryCounts {
            vertices: 0,
            edges: 1,
            faces: 0,
            loops: 0,
        };
        if actual_delta != expected_delta {
            return Err(CoreError::PreconditionFailed(format!(
                "same-loop boundary exposure postcondition failed: expected_delta=dV=0 dE=1 dF=0 dL=0 actual_delta=dV={} dE={} dF={} dL={}",
                actual_delta.vertices, actual_delta.edges, actual_delta.faces, actual_delta.loops
            )));
        }
        for edge in [result.first_boundary_edge, result.second_boundary_edge] {
            if !self.boundary_edges.contains(&edge) {
                return Err(CoreError::PreconditionFailed(
                    "same-loop exposed edge is not boundary".into(),
                ));
            }
        }
        let first = self.half_edges.get(result.first_half_edge).ok_or_else(|| {
            CoreError::PreconditionFailed(
                "same-loop exposed requested half-edge no longer exists".into(),
            )
        })?;
        let second = self
            .half_edges
            .get(result.second_half_edge)
            .ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "same-loop exposed twin half-edge no longer exists".into(),
                )
            })?;
        if first.twin != HalfEdgeHandle::default() || second.twin != HalfEdgeHandle::default() {
            return Err(CoreError::PreconditionFailed(
                "same-loop exposed half-edge kept a live twin".into(),
            ));
        }
        if first.edge != result.first_boundary_edge || second.edge != result.second_boundary_edge {
            return Err(CoreError::PreconditionFailed(
                "same-loop exposed half-edge edge ownership is wrong".into(),
            ));
        }
        self.verify_same_loop_boundary_loop_cycle(result)?;
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)?;
        Ok(())
    }

    fn verify_same_loop_boundary_loop_cycle(
        &self,
        result: SameLoopBoundaryExposureResult,
    ) -> Result<(), CoreError> {
        let first = self.half_edges.get(result.first_half_edge).ok_or_else(|| {
            CoreError::PreconditionFailed(
                "same-loop exposed requested half-edge is missing from touched loop".into(),
            )
        })?;
        let second = self
            .half_edges
            .get(result.second_half_edge)
            .ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "same-loop exposed twin half-edge is missing from touched loop".into(),
                )
            })?;
        if first.loop_ != second.loop_ {
            return Err(CoreError::PreconditionFailed(
                "same-loop exposed half-edges no longer share the touched loop".into(),
            ));
        }
        let loop_ = first.loop_;
        let loop_data = self.loops.get(loop_).ok_or_else(|| {
            CoreError::PreconditionFailed(
                "same-loop exposed half-edges reference a retired loop".into(),
            )
        })?;
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut visited = BTreeSet::new();
        for _ in 0..=self.half_edges.len() {
            let half_edge = self.half_edges.get(cursor).ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "same-loop touched loop references a missing half-edge".into(),
                )
            })?;
            if half_edge.loop_ != loop_ {
                return Err(CoreError::PreconditionFailed(
                    "same-loop touched cycle left the operand loop".into(),
                ));
            }
            if !visited.insert(cursor) {
                return Err(CoreError::PreconditionFailed(
                    "same-loop touched cycle repeats before closing".into(),
                ));
            }
            cursor = half_edge.next;
            if cursor == start {
                if visited.contains(&result.first_half_edge)
                    && visited.contains(&result.second_half_edge)
                {
                    return Ok(());
                }
                return Err(CoreError::PreconditionFailed(
                    "same-loop touched cycle does not contain both exposed half-edges".into(),
                ));
            }
        }
        Err(CoreError::PreconditionFailed(
            "same-loop touched cycle did not close".into(),
        ))
    }
}

impl SameLoopBoundaryCounts {
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
    use crate::geometry::Point3;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn expose_same_loop_edge_as_boundary_preserves_operand_half_edges() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let vertex = brep.vertices.keys().next().expect("cube has vertices");
        let outgoing = brep.vertices[vertex].outgoing;
        let spur = brep
            .mev(outgoing, Point3::new(0.5, 0.5, 0.5))
            .expect("same-loop spur is created");
        let loop_ = brep.half_edges[spur.half_edge].loop_;
        assert_eq!(loop_, brep.half_edges[spur.half_edge_twin].loop_);
        let before = (
            brep.vertices.len(),
            brep.edges.len(),
            brep.faces.len(),
            brep.loops.len(),
            brep.half_edges.len(),
        );

        let result = brep
            .expose_same_loop_edge_as_boundary(spur.half_edge)
            .expect("same-loop edge exposes as two boundary sides");

        assert_eq!(
            (
                brep.vertices.len() as i64 - before.0 as i64,
                brep.edges.len() as i64 - before.1 as i64,
                brep.faces.len() as i64 - before.2 as i64,
                brep.loops.len() as i64 - before.3 as i64,
                brep.half_edges.len() as i64 - before.4 as i64,
            ),
            (0, 1, 0, 0, 0)
        );
        assert_eq!(result.first_half_edge, spur.half_edge);
        assert_eq!(result.second_half_edge, spur.half_edge_twin);
        assert_eq!(brep.half_edges[spur.half_edge].loop_, loop_);
        assert_eq!(brep.half_edges[spur.half_edge_twin].loop_, loop_);
        assert!(brep.boundary_edges.contains(&result.first_boundary_edge));
        assert!(brep.boundary_edges.contains(&result.second_boundary_edge));
        assert_eq!(brep.half_edges[spur.half_edge].twin, Default::default());
        assert_eq!(
            brep.half_edges[spur.half_edge_twin].twin,
            Default::default()
        );
        brep.validate_half_edge_pairing().unwrap();
        brep.validate_half_edge_cycles().unwrap();
        brep.validate_manifold_edges().unwrap();
    }

    #[test]
    fn expose_same_loop_edge_as_boundary_rolls_back_for_cross_face_edge() {
        let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let edge = brep.edges.keys().next().expect("cube has edges");
        let half_edge = brep.edges[edge].half_edge;
        let before = (brep.edges.len(), brep.boundary_edges.len());

        let err = brep
            .expose_same_loop_edge_as_boundary(half_edge)
            .expect_err("cross-face edge is not a same-loop operand");

        assert!(err.to_string().contains("one loop"));
        assert_eq!(before, (brep.edges.len(), brep.boundary_edges.len()));
    }
}
