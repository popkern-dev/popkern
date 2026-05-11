//! Same-loop bridge-face construction for disjoint boundary intervals.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::topology::{
    BRep, Edge, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle,
};

/// One bridge face between two boundary intervals on the same source loop.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SameLoopFaceBridgeRequest {
    /// Source loop containing both boundary intervals.
    pub loop_: LoopHandle,
    /// First directed interval, including the after-last half-edge.
    pub first_interval: Vec<HalfEdgeHandle>,
    /// Second directed interval, including the after-last half-edge.
    pub second_interval: Vec<HalfEdgeHandle>,
}

/// Validated same-loop bridge-face plan.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SameLoopFaceBridgePlan {
    /// Source loop containing both boundary intervals.
    pub loop_: LoopHandle,
    /// Source face owning the source loop.
    pub source_face: FaceHandle,
    /// First directed interval, including the after-last half-edge.
    pub first_interval: Vec<HalfEdgeHandle>,
    /// Second directed interval, including the after-last half-edge.
    pub second_interval: Vec<HalfEdgeHandle>,
    /// Expected topology delta.
    pub expected_delta: SameLoopFaceBridgeTopologyDelta,
}

/// Expected topology count delta for a same-loop bridge-face mutation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SameLoopFaceBridgeTopologyDelta {
    /// Vertex count delta.
    pub vertices: i64,
    /// Edge count delta.
    pub edges: i64,
    /// Face count delta.
    pub faces: i64,
    /// Loop count delta.
    pub loops: i64,
}

/// Result of a committed same-loop bridge-face mutation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SameLoopFaceBridgeResult {
    /// First connector edge.
    pub first_connector_edge: EdgeHandle,
    /// Second connector edge.
    pub second_connector_edge: EdgeHandle,
    /// First connector half-edge.
    pub first_connector_half_edge: HalfEdgeHandle,
    /// Second connector half-edge.
    pub second_connector_half_edge: HalfEdgeHandle,
    /// New bridge face.
    pub bridge_face: FaceHandle,
    /// New bridge loop.
    pub bridge_loop: LoopHandle,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SameLoopBridgeCounts {
    vertices: i64,
    edges: i64,
    faces: i64,
    loops: i64,
}

impl BRep {
    /// Validate two disjoint boundary intervals on one loop as a bridge face.
    pub fn preflight_same_loop_face_bridge(
        &self,
        request: SameLoopFaceBridgeRequest,
    ) -> Result<SameLoopFaceBridgePlan, CoreError> {
        let source_face = self
            .loops
            .get(request.loop_)
            .map(|loop_data| loop_data.face)
            .ok_or_else(|| CoreError::InvalidHandle("same-loop bridge loop is invalid".into()))?;
        self.require_same_loop_bridge_interval(request.loop_, &request.first_interval, "first")?;
        self.require_same_loop_bridge_interval(request.loop_, &request.second_interval, "second")?;
        require_disjoint_same_loop_bridge_paths(
            self,
            &request.first_interval,
            &request.second_interval,
        )?;
        require_non_degenerate_same_loop_connectors(
            self,
            &request.first_interval,
            &request.second_interval,
        )?;
        Ok(SameLoopFaceBridgePlan {
            loop_: request.loop_,
            source_face,
            first_interval: request.first_interval,
            second_interval: request.second_interval,
            expected_delta: SameLoopFaceBridgeTopologyDelta {
                vertices: 0,
                edges: 2,
                faces: 1,
                loops: 1,
            },
        })
    }

    /// Transactionally author a bridge face between two intervals on one loop.
    pub fn same_loop_face_bridge(
        &mut self,
        plan: &SameLoopFaceBridgePlan,
    ) -> Result<SameLoopFaceBridgeResult, CoreError> {
        let before = SameLoopBridgeCounts::capture(self);
        let mut candidate = self.clone();
        let result = candidate.apply_same_loop_face_bridge_candidate(plan)?;
        candidate.verify_same_loop_face_bridge_postconditions(plan, before, result)?;
        *self = candidate;
        Ok(result)
    }

    fn require_same_loop_bridge_interval(
        &self,
        loop_: LoopHandle,
        interval: &[HalfEdgeHandle],
        label: &'static str,
    ) -> Result<(), CoreError> {
        let path = same_loop_bridge_path(interval)?;
        for &half_edge in interval {
            let data = self.half_edges.get(half_edge).ok_or_else(|| {
                CoreError::InvalidHandle(format!(
                    "same-loop bridge {label} interval has stale half-edge"
                ))
            })?;
            if data.loop_ != loop_ {
                return Err(CoreError::PreconditionFailed(format!(
                    "same-loop bridge {label} interval leaves requested loop"
                )));
            }
        }
        for pair in interval.windows(2) {
            if self.half_edges[pair[0]].next != pair[1] {
                return Err(CoreError::PreconditionFailed(format!(
                    "same-loop bridge {label} interval is not directed"
                )));
            }
        }
        for &half_edge in path {
            let data = &self.half_edges[half_edge];
            if !self.boundary_edges.contains(&data.edge) || self.half_edges.contains_key(data.twin)
            {
                return Err(CoreError::PreconditionFailed(format!(
                    "same-loop bridge {label} path half-edge {half_edge:?} is not an unowned boundary side"
                )));
            }
        }
        Ok(())
    }

    fn apply_same_loop_face_bridge_candidate(
        &mut self,
        plan: &SameLoopFaceBridgePlan,
    ) -> Result<SameLoopFaceBridgeResult, CoreError> {
        if self.loops[plan.loop_].face != plan.source_face {
            return Err(CoreError::PreconditionFailed(
                "same-loop bridge plan is stale".into(),
            ));
        }
        let first_path = self.create_same_loop_bridge_boundary_twins(&plan.first_interval)?;
        let second_path = self.create_same_loop_bridge_boundary_twins(&plan.second_interval)?;
        let bridge_loop = self.insert_same_loop_bridge_face_loop(plan.source_face)?;
        self.splice_same_loop_bridge_cycle(bridge_loop, &first_path, &second_path)
    }

    fn create_same_loop_bridge_boundary_twins(
        &mut self,
        interval: &[HalfEdgeHandle],
    ) -> Result<Vec<HalfEdgeHandle>, CoreError> {
        let mut twins = Vec::new();
        for &half_edge in same_loop_bridge_path(interval)? {
            let data = &self.half_edges[half_edge];
            let edge = data.edge;
            let origin = self.half_edges[data.next].origin;
            if !self.boundary_edges.contains(&edge) || self.half_edges.contains_key(data.twin) {
                return Err(CoreError::PreconditionFailed(format!(
                    "same-loop bridge path half-edge {half_edge:?} is not an unowned boundary side"
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

    fn insert_same_loop_bridge_face_loop(
        &mut self,
        source_face: FaceHandle,
    ) -> Result<LoopHandle, CoreError> {
        let source = self
            .faces
            .get(source_face)
            .ok_or_else(|| CoreError::InvalidHandle("same-loop bridge face is invalid".into()))?
            .clone();
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

    fn splice_same_loop_bridge_cycle(
        &mut self,
        bridge_loop: LoopHandle,
        first_path: &[HalfEdgeHandle],
        second_path: &[HalfEdgeHandle],
    ) -> Result<SameLoopFaceBridgeResult, CoreError> {
        let first_start = required_same_loop_path_handle(first_path.first(), "first", "start")?;
        let first_end = required_same_loop_path_handle(first_path.last(), "first", "end")?;
        let second_start = required_same_loop_path_handle(second_path.first(), "second", "start")?;
        let second_end = required_same_loop_path_handle(second_path.last(), "second", "end")?;

        let first_connector = self.insert_same_loop_bridge_connector(
            self.same_loop_bridge_path_terminal_vertex(first_end)?,
            bridge_loop,
        );
        let second_connector = self.insert_same_loop_bridge_connector(
            self.same_loop_bridge_path_terminal_vertex(second_end)?,
            bridge_loop,
        );
        wire_same_loop_bridge_path(self, first_path);
        wire_same_loop_bridge_path(self, second_path);
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
        Ok(SameLoopFaceBridgeResult {
            first_connector_edge: first_connector.0,
            second_connector_edge: second_connector.0,
            first_connector_half_edge: first_connector.1,
            second_connector_half_edge: second_connector.1,
            bridge_face,
            bridge_loop,
        })
    }

    fn insert_same_loop_bridge_connector(
        &mut self,
        origin: crate::topology::VertexHandle,
        loop_: LoopHandle,
    ) -> (EdgeHandle, HalfEdgeHandle) {
        let edge = self.edges.insert(Edge {
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

    fn same_loop_bridge_path_terminal_vertex(
        &self,
        path_half_edge: HalfEdgeHandle,
    ) -> Result<crate::topology::VertexHandle, CoreError> {
        let twin = self
            .half_edges
            .get(path_half_edge)
            .ok_or_else(|| {
                CoreError::InvalidHandle("same-loop bridge path half-edge is invalid".into())
            })?
            .twin;
        self.half_edges
            .get(twin)
            .map(|half_edge| half_edge.origin)
            .ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "same-loop bridge path half-edge has no source-side twin".into(),
                )
            })
    }

    fn verify_same_loop_face_bridge_postconditions(
        &self,
        plan: &SameLoopFaceBridgePlan,
        before: SameLoopBridgeCounts,
        result: SameLoopFaceBridgeResult,
    ) -> Result<(), CoreError> {
        let actual = SameLoopBridgeCounts::capture(self).delta_from(before);
        let expected = SameLoopBridgeCounts {
            vertices: plan.expected_delta.vertices,
            edges: plan.expected_delta.edges,
            faces: plan.expected_delta.faces,
            loops: plan.expected_delta.loops,
        };
        if actual != expected {
            return Err(CoreError::PreconditionFailed(format!(
                "same-loop bridge postcondition failed: expected_delta=dV={} dE={} dF={} dL={} actual_delta=dV={} dE={} dF={} dL={}",
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
        verify_same_loop_bridge_loop_cycle(self, result.bridge_loop)?;
        if !self.boundary_edges.contains(&result.first_connector_edge)
            || !self.boundary_edges.contains(&result.second_connector_edge)
        {
            return Err(CoreError::PreconditionFailed(
                "same-loop bridge connector edge lost boundary ownership".into(),
            ));
        }
        Ok(())
    }
}

fn same_loop_bridge_path(interval: &[HalfEdgeHandle]) -> Result<&[HalfEdgeHandle], CoreError> {
    if interval.len() < 2 {
        return Err(CoreError::PreconditionFailed(
            "same-loop bridge path needs distinct endpoint origins".into(),
        ));
    }
    Ok(&interval[..interval.len() - 1])
}

fn require_disjoint_same_loop_bridge_paths(
    brep: &BRep,
    first_interval: &[HalfEdgeHandle],
    second_interval: &[HalfEdgeHandle],
) -> Result<(), CoreError> {
    let mut consumed_edges = BTreeSet::new();
    for &half_edge in same_loop_bridge_path(first_interval)?
        .iter()
        .chain(same_loop_bridge_path(second_interval)?.iter())
    {
        let edge = brep.half_edges[half_edge].edge;
        if !consumed_edges.insert(edge) {
            return Err(CoreError::PreconditionFailed(format!(
                "same-loop bridge intervals overlap on boundary edge {edge:?}"
            )));
        }
    }
    Ok(())
}

fn require_non_degenerate_same_loop_connectors(
    brep: &BRep,
    first_interval: &[HalfEdgeHandle],
    second_interval: &[HalfEdgeHandle],
) -> Result<(), CoreError> {
    let first_path = same_loop_bridge_path(first_interval)?;
    let second_path = same_loop_bridge_path(second_interval)?;
    let first_start = required_same_loop_path_handle(first_path.first(), "first", "start")?;
    let first_end = required_same_loop_path_handle(first_path.last(), "first", "end")?;
    let second_start = required_same_loop_path_handle(second_path.first(), "second", "start")?;
    let second_end = required_same_loop_path_handle(second_path.last(), "second", "end")?;
    let first_end_vertex = brep.half_edges[brep.half_edges[first_end].next].origin;
    let second_end_vertex = brep.half_edges[brep.half_edges[second_end].next].origin;
    if first_end_vertex == brep.half_edges[second_start].origin
        || second_end_vertex == brep.half_edges[first_start].origin
    {
        return Err(CoreError::PreconditionFailed(
            "same-loop bridge connector would be degenerate".into(),
        ));
    }
    Ok(())
}

fn required_same_loop_path_handle(
    handle: Option<&HalfEdgeHandle>,
    path: &'static str,
    end: &'static str,
) -> Result<HalfEdgeHandle, CoreError> {
    handle.copied().ok_or_else(|| {
        CoreError::PreconditionFailed(format!("{path} same-loop bridge path {end} is empty"))
    })
}

fn wire_same_loop_bridge_path(brep: &mut BRep, path: &[HalfEdgeHandle]) {
    for pair in path.windows(2) {
        brep.half_edges[pair[0]].next = pair[1];
        brep.half_edges[pair[1]].prev = pair[0];
    }
}

fn verify_same_loop_bridge_loop_cycle(brep: &BRep, loop_: LoopHandle) -> Result<(), CoreError> {
    let start = brep.loops[loop_].half_edge;
    let mut cursor = start;
    let mut visited = BTreeSet::new();
    for _ in 0..=brep.half_edges.len() {
        let half_edge = &brep.half_edges[cursor];
        if half_edge.loop_ != loop_ {
            return Err(CoreError::PreconditionFailed(
                "same-loop bridge cycle left its loop".into(),
            ));
        }
        if !visited.insert(cursor) {
            return Err(CoreError::PreconditionFailed(
                "same-loop bridge cycle repeats before closing".into(),
            ));
        }
        cursor = half_edge.next;
        if cursor == start {
            return Ok(());
        }
    }
    Err(CoreError::PreconditionFailed(
        "same-loop bridge cycle did not close".into(),
    ))
}

impl SameLoopBridgeCounts {
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
#[path = "same_loop_face_bridge_tests.rs"]
mod tests;
