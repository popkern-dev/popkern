//! Aggregate bridge-face construction for disjoint boundary cycles.
use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::topology::{
    BRep, Edge, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle,
};
/// One bridge face cycle in an aggregate face-family bridge.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaceFamilyBridgeCycleRequest {
    /// First source loop contributing a boundary interval.
    pub first_loop: LoopHandle,
    /// Directed interval on `first_loop`.
    pub first_interval: Vec<HalfEdgeHandle>,
    /// Second source loop contributing the opposite boundary interval.
    pub second_loop: LoopHandle,
    /// Directed interval on `second_loop`, oriented opposite to the intended
    /// new bridge-face traversal.
    pub second_interval: Vec<HalfEdgeHandle>,
}

/// Validated aggregate bridge plan.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaceFamilyBridgePlan {
    /// Bridge cycles to author in one transaction.
    pub cycles: Vec<FaceFamilyBridgeCycleRequest>,
    /// Expected topology delta for the aggregate mutation.
    pub expected_delta: FaceFamilyBridgeTopologyDelta,
}

/// Expected topology count delta for an aggregate bridge mutation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FaceFamilyBridgeTopologyDelta {
    /// Vertex count delta.
    pub vertices: i64,
    /// Edge count delta.
    pub edges: i64,
    /// Face count delta.
    pub faces: i64,
    /// Loop count delta.
    pub loops: i64,
}

/// Result of one committed aggregate bridge cycle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FaceFamilyBridgeCycleResult {
    /// First connector edge.
    pub first_connector_edge: EdgeHandle,
    /// Second connector edge.
    pub second_connector_edge: EdgeHandle,
    /// First connector half-edge in the new bridge face.
    pub first_connector_half_edge: HalfEdgeHandle,
    /// Second connector half-edge in the new bridge face.
    pub second_connector_half_edge: HalfEdgeHandle,
    /// New bridge face.
    pub bridge_face: FaceHandle,
    /// New bridge loop.
    pub bridge_loop: LoopHandle,
}

/// Result of a committed aggregate bridge mutation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaceFamilyBridgeResult {
    /// Per-cycle committed bridge faces.
    pub cycles: Vec<FaceFamilyBridgeCycleResult>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TopologyCounts {
    vertices: i64,
    edges: i64,
    faces: i64,
    loops: i64,
}

impl BRep {
    /// Validate disjoint boundary cycles for an aggregate bridge-face mutation.
    ///
    /// The aggregate operator consumes each source boundary edge at most once.
    /// Shared-loop pairwise plans that repeat a long outside interval are
    /// rejected here because one manifold half-edge side cannot bound multiple
    /// new bridge faces.
    pub fn preflight_face_family_bridge(
        &self,
        cycles: Vec<FaceFamilyBridgeCycleRequest>,
    ) -> Result<FaceFamilyBridgePlan, CoreError> {
        if cycles.is_empty() {
            return Err(CoreError::PreconditionFailed(
                "face-family bridge requires at least one cycle".into(),
            ));
        }
        let mut consumed_edges = BTreeSet::new();
        for (index, cycle) in cycles.iter().enumerate() {
            self.require_bridge_cycle_interval(index, cycle.first_loop, &cycle.first_interval)?;
            self.require_bridge_cycle_interval(index, cycle.second_loop, &cycle.second_interval)?;
            if cycle.first_loop == cycle.second_loop {
                return Err(CoreError::PreconditionFailed(format!(
                    "face-family bridge cycle {index} requires distinct source loops"
                )));
            }
            let first_face = self.loops[cycle.first_loop].face;
            let second_face = self.loops[cycle.second_loop].face;
            if self.faces[first_face].shell != self.faces[second_face].shell {
                return Err(CoreError::PreconditionFailed(format!(
                    "face-family bridge cycle {index} source faces are in different shells"
                )));
            }
            for &half_edge in half_open_path(&cycle.first_interval)?
                .iter()
                .chain(half_open_path(&cycle.second_interval)?.iter())
            {
                let edge = self.half_edges[half_edge].edge;
                if !consumed_edges.insert(edge) {
                    return Err(CoreError::PreconditionFailed(format!(
                        "face-family bridge cycles overlap on boundary edge {edge:?}"
                    )));
                }
            }
        }
        let cycle_count = cycles.len() as i64;
        Ok(FaceFamilyBridgePlan {
            cycles,
            expected_delta: FaceFamilyBridgeTopologyDelta {
                vertices: 0,
                edges: cycle_count * 2,
                faces: cycle_count,
                loops: cycle_count,
            },
        })
    }

    /// Transactionally author aggregate bridge faces from a validated plan.
    pub fn face_family_bridge(
        &mut self,
        plan: &FaceFamilyBridgePlan,
    ) -> Result<FaceFamilyBridgeResult, CoreError> {
        let before = TopologyCounts::capture(self);
        let mut candidate = self.clone();
        let result = candidate.apply_face_family_bridge_candidate(plan)?;
        candidate.verify_face_family_bridge_postconditions(plan, before, &result)?;
        *self = candidate;
        Ok(result)
    }

    fn require_bridge_cycle_interval(
        &self,
        index: usize,
        loop_handle: LoopHandle,
        interval: &[HalfEdgeHandle],
    ) -> Result<(), CoreError> {
        let path = half_open_path(interval)?;
        for &half_edge in interval {
            let data = self.half_edges.get(half_edge).ok_or_else(|| {
                CoreError::InvalidHandle(format!(
                    "face-family bridge cycle {index} has stale interval half-edge"
                ))
            })?;
            if data.loop_ != loop_handle {
                return Err(CoreError::PreconditionFailed(format!(
                    "face-family bridge cycle {index} interval leaves requested loop"
                )));
            }
        }
        for pair in interval.windows(2) {
            if self.half_edges[pair[0]].next != pair[1] {
                return Err(CoreError::PreconditionFailed(format!(
                    "face-family bridge cycle {index} interval is not directed through the requested loop"
                )));
            }
        }
        for &half_edge in path {
            let data = &self.half_edges[half_edge];
            if !self.boundary_edges.contains(&data.edge) || self.half_edges.contains_key(data.twin)
            {
                return Err(CoreError::PreconditionFailed(format!(
                    "face-family bridge cycle {index} path half-edge {half_edge:?} is not an unowned boundary side"
                )));
            }
        }
        Ok(())
    }

    fn apply_face_family_bridge_candidate(
        &mut self,
        plan: &FaceFamilyBridgePlan,
    ) -> Result<FaceFamilyBridgeResult, CoreError> {
        let mut results = Vec::with_capacity(plan.cycles.len());
        for cycle in &plan.cycles {
            let first_path = self.create_bridge_boundary_twins(&cycle.first_interval)?;
            let second_path = self.create_bridge_boundary_twins(&cycle.second_interval)?;
            let bridge_loop = self.insert_bridge_cycle_face(cycle.first_loop)?;
            let result =
                self.splice_aggregate_bridge_cycle(bridge_loop, &first_path, &second_path)?;
            results.push(result);
        }
        Ok(FaceFamilyBridgeResult { cycles: results })
    }

    fn create_bridge_boundary_twins(
        &mut self,
        interval: &[HalfEdgeHandle],
    ) -> Result<Vec<HalfEdgeHandle>, CoreError> {
        let mut twins = Vec::new();
        for &half_edge in half_open_path(interval)? {
            let data = &self.half_edges[half_edge];
            let edge = data.edge;
            let origin = self.half_edges[data.next].origin;
            if !self.boundary_edges.contains(&edge) || self.half_edges.contains_key(data.twin) {
                return Err(CoreError::PreconditionFailed(format!(
                    "face-family bridge path half-edge {half_edge:?} is not an unowned boundary side"
                )));
            }
            let twin = self.half_edges.insert(HalfEdge {
                origin,
                twin: half_edge,
                next: Default::default(),
                prev: Default::default(),
                edge,
                loop_: Default::default(),
                pcurve: None,
            });
            self.half_edges[half_edge].twin = twin;
            self.boundary_edges.remove(&edge);
            twins.push(twin);
        }
        twins.reverse();
        Ok(twins)
    }

    fn insert_bridge_cycle_face(
        &mut self,
        source_loop: LoopHandle,
    ) -> Result<LoopHandle, CoreError> {
        let source_face = self.loops[source_loop].face;
        let source = self.faces[source_face].clone();
        let bridge_face = self.faces.insert(Face {
            outer_loop: Default::default(),
            inner_loops: vec![],
            shell: source.shell,
            surface: None,
            name: None,
            tolerance: source.tolerance,
            role: Default::default(),
        });
        let bridge_loop = self.loops.insert(Loop {
            half_edge: Default::default(),
            face: bridge_face,
        });
        self.faces[bridge_face].outer_loop = bridge_loop;
        Ok(bridge_loop)
    }

    fn splice_aggregate_bridge_cycle(
        &mut self,
        bridge_loop: LoopHandle,
        first_path: &[HalfEdgeHandle],
        second_path: &[HalfEdgeHandle],
    ) -> Result<FaceFamilyBridgeCycleResult, CoreError> {
        let first_start = required_path_handle(first_path.first(), "first", "start")?;
        let first_end = required_path_handle(first_path.last(), "first", "end")?;
        let second_start = required_path_handle(second_path.first(), "second", "start")?;
        let second_end = required_path_handle(second_path.last(), "second", "end")?;

        let first_connector = self
            .insert_boundary_connector(self.bridge_path_terminal_vertex(first_end)?, bridge_loop);
        let second_connector = self
            .insert_boundary_connector(self.bridge_path_terminal_vertex(second_end)?, bridge_loop);
        self.wire_aggregate_bridge_path(first_path);
        self.wire_aggregate_bridge_path(second_path);
        self.half_edges[first_end].next = first_connector.1;
        self.half_edges[first_connector.1].prev = first_end;
        self.half_edges[first_connector.1].next = second_start;
        self.half_edges[second_start].prev = first_connector.1;
        self.half_edges[second_end].next = second_connector.1;
        self.half_edges[second_connector.1].prev = second_end;
        self.half_edges[second_connector.1].next = first_start;
        self.half_edges[first_start].prev = second_connector.1;

        for half_edge in first_path
            .iter()
            .chain(second_path.iter())
            .chain([first_connector.1, second_connector.1].iter())
        {
            self.half_edges[*half_edge].loop_ = bridge_loop;
        }
        self.loops[bridge_loop].half_edge = first_start;
        let bridge_face = self.loops[bridge_loop].face;
        Ok(FaceFamilyBridgeCycleResult {
            first_connector_edge: first_connector.0,
            second_connector_edge: second_connector.0,
            first_connector_half_edge: first_connector.1,
            second_connector_half_edge: second_connector.1,
            bridge_face,
            bridge_loop,
        })
    }

    fn insert_boundary_connector(
        &mut self,
        origin: crate::topology::VertexHandle,
        loop_: LoopHandle,
    ) -> (EdgeHandle, HalfEdgeHandle) {
        let edge = self.insert_edge_tracked(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });
        let half_edge = self.half_edges.insert(HalfEdge {
            origin,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge,
            loop_,
            pcurve: None,
        });
        self.edges[edge].half_edge = half_edge;
        self.boundary_edges.insert(edge);
        (edge, half_edge)
    }

    fn wire_aggregate_bridge_path(&mut self, path: &[HalfEdgeHandle]) {
        for pair in path.windows(2) {
            self.half_edges[pair[0]].next = pair[1];
            self.half_edges[pair[1]].prev = pair[0];
        }
    }

    fn bridge_path_terminal_vertex(
        &self,
        path_half_edge: HalfEdgeHandle,
    ) -> Result<crate::topology::VertexHandle, CoreError> {
        let twin = self
            .half_edges
            .get(path_half_edge)
            .ok_or_else(|| CoreError::InvalidHandle("bridge path half-edge is invalid".into()))?
            .twin;
        self.half_edges
            .get(twin)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "face-family bridge path half-edge has no source-side twin".into(),
                )
            })
    }

    fn verify_face_family_bridge_postconditions(
        &self,
        plan: &FaceFamilyBridgePlan,
        before: TopologyCounts,
        result: &FaceFamilyBridgeResult,
    ) -> Result<(), CoreError> {
        let expected = TopologyCounts {
            vertices: plan.expected_delta.vertices,
            edges: plan.expected_delta.edges,
            faces: plan.expected_delta.faces,
            loops: plan.expected_delta.loops,
        };
        let actual = TopologyCounts::capture(self).delta_from(before);
        if actual != expected {
            return Err(CoreError::PreconditionFailed(format!(
                "face-family bridge postcondition failed: expected_delta=dV={} dE={} dF={} dL={} actual_delta=dV={} dE={} dF={} dL={}",
                expected.vertices,
                expected.edges,
                expected.faces,
                expected.loops,
                actual.vertices,
                actual.edges,
                actual.faces,
                actual.loops
            )));
        }
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)?;
        for cycle in &result.cycles {
            self.verify_bridge_loop_cycle(cycle.bridge_loop)?;
            self.verify_bridge_loop_vertex_continuity(cycle.bridge_loop)?;
            if !self.boundary_edges.contains(&cycle.first_connector_edge)
                || !self.boundary_edges.contains(&cycle.second_connector_edge)
            {
                return Err(CoreError::PreconditionFailed(
                    "face-family bridge connector edge lost boundary ownership".into(),
                ));
            }
        }
        Ok(())
    }

    fn verify_bridge_loop_cycle(&self, loop_: LoopHandle) -> Result<(), CoreError> {
        let start = self.loops[loop_].half_edge;
        let mut cursor = start;
        let mut visited = BTreeSet::new();
        for _ in 0..=self.half_edges.len() {
            let half_edge = &self.half_edges[cursor];
            if half_edge.loop_ != loop_ {
                return Err(CoreError::PreconditionFailed(
                    "face-family bridge cycle left its loop".into(),
                ));
            }
            if !visited.insert(cursor) {
                return Err(CoreError::PreconditionFailed(
                    "face-family bridge cycle repeats before closing".into(),
                ));
            }
            cursor = half_edge.next;
            if cursor == start {
                return Ok(());
            }
        }
        Err(CoreError::PreconditionFailed(
            "face-family bridge cycle did not close".into(),
        ))
    }

    fn verify_bridge_loop_vertex_continuity(&self, loop_: LoopHandle) -> Result<(), CoreError> {
        let start = self.loops[loop_].half_edge;
        let mut cursor = start;
        for _ in 0..=self.half_edges.len() {
            let half_edge = &self.half_edges[cursor];
            let next = self.half_edges.get(half_edge.next).ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "face-family bridge cycle has stale next half-edge".into(),
                )
            })?;
            let expected_next_origin = if self.boundary_edges.contains(&half_edge.edge) {
                next.origin
            } else {
                self.half_edges
                    .get(half_edge.twin)
                    .map(|twin| twin.origin)
                    .ok_or_else(|| {
                        CoreError::PreconditionFailed(
                            "face-family bridge non-boundary half-edge has stale twin".into(),
                        )
                    })?
            };
            if expected_next_origin != next.origin {
                return Err(CoreError::PreconditionFailed(format!(
                    "face-family bridge cycle has discontinuous vertex adjacency at {cursor:?}"
                )));
            }
            cursor = half_edge.next;
            if cursor == start {
                return Ok(());
            }
        }
        Err(CoreError::PreconditionFailed(
            "face-family bridge cycle continuity check did not close".into(),
        ))
    }
}

impl TopologyCounts {
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

fn half_open_path(interval: &[HalfEdgeHandle]) -> Result<&[HalfEdgeHandle], CoreError> {
    if interval.len() < 2 {
        return Err(CoreError::PreconditionFailed(
            "face-family bridge path needs distinct endpoint origins".into(),
        ));
    }
    Ok(&interval[..interval.len() - 1])
}

fn required_path_handle(
    handle: Option<&HalfEdgeHandle>,
    path: &'static str,
    end: &'static str,
) -> Result<HalfEdgeHandle, CoreError> {
    handle.copied().ok_or_else(|| {
        CoreError::PreconditionFailed(format!("{path} aggregate bridge path {end} is empty"))
    })
}

#[cfg(test)]
#[path = "face_family_bridge_tests.rs"]
mod tests;
