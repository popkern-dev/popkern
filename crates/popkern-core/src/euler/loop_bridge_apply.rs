//! Transactional two-loop bridge mutation.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::euler::TwoLoopBridgePlan;
use crate::topology::{
    BRep, Edge, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle,
};

/// Result of a committed two-loop bridge mutation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TwoLoopBridgeResult {
    /// New edge connecting the first-loop start endpoint to the second-loop start endpoint.
    pub start_bridge_edge: EdgeHandle,
    /// New edge connecting the first-loop end endpoint to the second-loop end endpoint.
    pub end_bridge_edge: EdgeHandle,
    /// Half-edge from the first-loop start endpoint to the second-loop start endpoint.
    pub start_bridge_half_edge: HalfEdgeHandle,
    /// Twin half-edge from the second-loop start endpoint to the first-loop start endpoint.
    pub start_bridge_half_edge_twin: HalfEdgeHandle,
    /// Half-edge from the first-loop end endpoint to the second-loop end endpoint.
    pub end_bridge_half_edge: HalfEdgeHandle,
    /// Twin half-edge from the second-loop end endpoint to the first-loop end endpoint.
    pub end_bridge_half_edge_twin: HalfEdgeHandle,
    /// New face owning the first planned bridge cycle.
    pub first_bridge_face: FaceHandle,
    /// New loop owning the first planned bridge cycle.
    pub first_bridge_loop: LoopHandle,
    /// New face owning the second planned bridge cycle.
    pub second_bridge_face: FaceHandle,
    /// New loop owning the second planned bridge cycle.
    pub second_bridge_loop: LoopHandle,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TopologyCounts {
    vertices: i64,
    edges: i64,
    faces: i64,
    loops: i64,
}

impl BRep {
    /// Transactionally bridge two loops using a validated preflight plan.
    ///
    /// The operator owns the live mutation boundary: it builds a candidate on a
    /// cloned BRep, verifies the requested Euler delta and touched half-edge
    /// ownership, then swaps the candidate into place only after those
    /// postconditions pass.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if any handle captured by the plan is stale.
    /// - [`CoreError::PreconditionFailed`] if the plan is internally inconsistent
    ///   with the live model, or the candidate fails postcondition checks.
    pub fn two_loop_bridge(
        &mut self,
        plan: &TwoLoopBridgePlan,
    ) -> Result<TwoLoopBridgeResult, CoreError> {
        let before = TopologyCounts::capture(self);
        let mut candidate = self.clone();
        let result = candidate.apply_two_loop_bridge_candidate(plan)?;
        candidate.verify_two_loop_bridge_postconditions(plan, before, result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_two_loop_bridge_candidate(
        &mut self,
        plan: &TwoLoopBridgePlan,
    ) -> Result<TwoLoopBridgeResult, CoreError> {
        let first_start = *plan.first_forward_interval.first().ok_or_else(|| {
            CoreError::PreconditionFailed("first bridge interval is empty".into())
        })?;
        let first_end = *plan.first_forward_interval.last().ok_or_else(|| {
            CoreError::PreconditionFailed("first bridge interval is empty".into())
        })?;
        let second_start = *plan.second_forward_interval.first().ok_or_else(|| {
            CoreError::PreconditionFailed("second bridge interval is empty".into())
        })?;
        let second_end = *plan.second_forward_interval.last().ok_or_else(|| {
            CoreError::PreconditionFailed("second bridge interval is empty".into())
        })?;

        let first_loop = self.half_edge_loop(first_start, "first_start")?;
        let second_loop = self.half_edge_loop(second_start, "second_start")?;
        if self.half_edge_loop(first_end, "first_end")? != first_loop
            || self.half_edge_loop(second_end, "second_end")? != second_loop
        {
            return Err(CoreError::PreconditionFailed(
                "two-loop bridge endpoints are not on their planned loops".into(),
            ));
        }
        if self.faces[plan.first_face].shell != self.faces[plan.second_face].shell {
            return Err(CoreError::PreconditionFailed(
                "two-loop bridge requires operand faces in the same shell".into(),
            ));
        }

        let start_bridge_edge = self.insert_unshaped_bridge_edge();
        let end_bridge_edge = self.insert_unshaped_bridge_edge();

        let first_start_vertex = self.half_edges[first_start].origin;
        let first_end_vertex = self.half_edges[first_end].origin;
        let second_start_vertex = self.half_edges[second_start].origin;
        let second_end_vertex = self.half_edges[second_end].origin;

        let start_bridge_half_edge = self.half_edges.insert(HalfEdge {
            origin: first_start_vertex,
            twin: Default::default(),
            next: second_start,
            prev: first_end,
            edge: start_bridge_edge,
            loop_: first_loop,
            pcurve: None,
        });
        let start_bridge_half_edge_twin = self.half_edges.insert(HalfEdge {
            origin: second_start_vertex,
            twin: start_bridge_half_edge,
            next: first_start,
            prev: second_end,
            edge: start_bridge_edge,
            loop_: second_loop,
            pcurve: None,
        });
        self.half_edges[start_bridge_half_edge].twin = start_bridge_half_edge_twin;

        let end_bridge_half_edge = self.half_edges.insert(HalfEdge {
            origin: first_end_vertex,
            twin: Default::default(),
            next: second_end,
            prev: first_start,
            edge: end_bridge_edge,
            loop_: second_loop,
            pcurve: None,
        });
        let end_bridge_half_edge_twin = self.half_edges.insert(HalfEdge {
            origin: second_end_vertex,
            twin: end_bridge_half_edge,
            next: first_end,
            prev: second_start,
            edge: end_bridge_edge,
            loop_: first_loop,
            pcurve: None,
        });
        self.half_edges[end_bridge_half_edge].twin = end_bridge_half_edge_twin;

        self.edges[start_bridge_edge].half_edge = start_bridge_half_edge;
        self.edges[end_bridge_edge].half_edge = end_bridge_half_edge;

        let first_bridge_loop = self.insert_bridge_face_loop(plan.first_face)?;
        let second_bridge_loop = self.insert_bridge_face_loop(plan.second_face)?;

        let first_forward_twins =
            self.create_boundary_twins_for_bridge_path(&plan.first_forward_interval)?;
        let first_reverse_twins =
            self.create_boundary_twins_for_bridge_path(&plan.first_reverse_interval)?;
        let second_forward_twins =
            self.create_boundary_twins_for_bridge_path(&plan.second_forward_interval)?;
        let second_reverse_twins =
            self.create_boundary_twins_for_bridge_path(&plan.second_reverse_interval)?;

        self.splice_bridge_face_cycle(
            first_bridge_loop,
            &first_forward_twins,
            start_bridge_half_edge,
            &second_reverse_twins,
            end_bridge_half_edge_twin,
        )?;
        self.splice_bridge_face_cycle(
            second_bridge_loop,
            &second_forward_twins,
            start_bridge_half_edge_twin,
            &first_reverse_twins,
            end_bridge_half_edge,
        )?;

        self.loops[first_loop].half_edge = plan.first_forward_interval[0];
        self.loops[second_loop].half_edge = plan.second_forward_interval[0];

        Ok(TwoLoopBridgeResult {
            start_bridge_edge,
            end_bridge_edge,
            start_bridge_half_edge,
            start_bridge_half_edge_twin,
            end_bridge_half_edge,
            end_bridge_half_edge_twin,
            first_bridge_face: self.loops[first_bridge_loop].face,
            first_bridge_loop,
            second_bridge_face: self.loops[second_bridge_loop].face,
            second_bridge_loop,
        })
    }

    fn insert_unshaped_bridge_edge(&mut self) -> EdgeHandle {
        self.edges.insert(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
        })
    }

    fn insert_bridge_face_loop(
        &mut self,
        source_face: FaceHandle,
    ) -> Result<LoopHandle, CoreError> {
        let source = self
            .faces
            .get(source_face)
            .ok_or_else(|| CoreError::InvalidHandle("source face is invalid".into()))?;
        let shell = source.shell;
        let tolerance = source.tolerance;
        let bridge_face = self.faces.insert(Face {
            outer_loop: Default::default(),
            inner_loops: vec![],
            shell,
            surface: None,
            name: None,
            tolerance,
        });
        let bridge_loop = self.loops.insert(Loop {
            half_edge: Default::default(),
            face: bridge_face,
        });
        self.faces[bridge_face].outer_loop = bridge_loop;
        Ok(bridge_loop)
    }

    fn create_boundary_twins_for_bridge_path(
        &mut self,
        interval: &[HalfEdgeHandle],
    ) -> Result<Vec<HalfEdgeHandle>, CoreError> {
        let path = half_open_path(interval)?;
        let mut twins = Vec::with_capacity(path.len());
        for &half_edge in path {
            twins.push(self.create_boundary_twin_for_bridge(half_edge)?);
        }
        twins.reverse();
        Ok(twins)
    }

    fn create_boundary_twin_for_bridge(
        &mut self,
        half_edge: HalfEdgeHandle,
    ) -> Result<HalfEdgeHandle, CoreError> {
        let data = self
            .half_edges
            .get(half_edge)
            .ok_or_else(|| CoreError::InvalidHandle("bridge path half-edge is invalid".into()))?;
        let edge = data.edge;
        let next = data.next;
        let existing_twin = data.twin;
        if !self.boundary_edges.contains(&edge) || self.half_edges.contains_key(existing_twin) {
            return Err(CoreError::PreconditionFailed(format!(
                "two-loop bridge path half-edge {half_edge:?} is not an unowned boundary side"
            )));
        }
        let origin = self.half_edges[next].origin;
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
        Ok(twin)
    }

    fn splice_bridge_face_cycle(
        &mut self,
        loop_handle: LoopHandle,
        first_path: &[HalfEdgeHandle],
        first_bridge: HalfEdgeHandle,
        second_path: &[HalfEdgeHandle],
        second_bridge: HalfEdgeHandle,
    ) -> Result<(), CoreError> {
        let first_path_start = required_path_handle(first_path.first(), "first", "start")?;
        let second_path_start = required_path_handle(second_path.first(), "second", "start")?;
        let first_path_end = required_path_handle(first_path.last(), "first", "end")?;
        let second_path_end = required_path_handle(second_path.last(), "second", "end")?;

        self.wire_bridge_path(first_path);
        self.wire_bridge_path(second_path);
        self.half_edges[first_path_end].next = first_bridge;
        self.half_edges[first_bridge].prev = first_path_end;
        self.half_edges[first_bridge].next = second_path_start;
        self.half_edges[second_path_start].prev = first_bridge;
        self.half_edges[second_path_end].next = second_bridge;
        self.half_edges[second_bridge].prev = second_path_end;
        self.half_edges[second_bridge].next = first_path_start;
        self.half_edges[first_path_start].prev = second_bridge;

        for half_edge in first_path
            .iter()
            .chain(second_path.iter())
            .chain([first_bridge, second_bridge].iter())
        {
            self.half_edges[*half_edge].loop_ = loop_handle;
        }
        self.loops[loop_handle].half_edge = first_path_start;
        Ok(())
    }

    fn wire_bridge_path(&mut self, path: &[HalfEdgeHandle]) {
        for pair in path.windows(2) {
            self.half_edges[pair[0]].next = pair[1];
            self.half_edges[pair[1]].prev = pair[0];
        }
    }

    fn verify_two_loop_bridge_postconditions(
        &self,
        plan: &TwoLoopBridgePlan,
        before: TopologyCounts,
        result: TwoLoopBridgeResult,
    ) -> Result<(), CoreError> {
        let actual_delta = TopologyCounts::capture(self).delta_from(before);
        let expected_delta = TopologyCounts {
            vertices: plan.expected_delta.vertices,
            edges: plan.expected_delta.edges,
            faces: plan.expected_delta.faces,
            loops: plan.expected_delta.loops,
        };
        if actual_delta != expected_delta {
            return Err(CoreError::PreconditionFailed(format!(
                "two-loop bridge postcondition failed: expected_delta=dV={} dE={} dF={} dL={} actual_delta=dV={} dE={} dF={} dL={}",
                expected_delta.vertices,
                expected_delta.edges,
                expected_delta.faces,
                expected_delta.loops,
                actual_delta.vertices,
                actual_delta.edges,
                actual_delta.faces,
                actual_delta.loops
            )));
        }
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)?;
        self.verify_two_loop_bridge_loop_cycle(
            self.half_edges[plan.first_forward_interval[0]].loop_,
            &[
                plan.first_forward_interval[0],
                *plan.first_forward_interval.last().ok_or_else(|| {
                    CoreError::PreconditionFailed("first bridge interval is empty".into())
                })?,
            ],
        )?;
        self.verify_two_loop_bridge_loop_cycle(
            self.half_edges[plan.second_forward_interval[0]].loop_,
            &[
                plan.second_forward_interval[0],
                *plan.second_forward_interval.last().ok_or_else(|| {
                    CoreError::PreconditionFailed("second bridge interval is empty".into())
                })?,
            ],
        )?;
        self.verify_two_loop_bridge_loop_cycle(
            result.first_bridge_loop,
            &[
                result.start_bridge_half_edge,
                result.end_bridge_half_edge_twin,
            ],
        )?;
        self.verify_two_loop_bridge_loop_cycle(
            result.second_bridge_loop,
            &[
                result.start_bridge_half_edge_twin,
                result.end_bridge_half_edge,
            ],
        )?;
        self.verify_bridge_cycle_ownership(&plan.boundary_cycles[0].first_interval)?;
        self.verify_bridge_cycle_ownership(&plan.boundary_cycles[0].second_interval)?;
        self.verify_bridge_cycle_ownership(&plan.boundary_cycles[1].first_interval)?;
        self.verify_bridge_cycle_ownership(&plan.boundary_cycles[1].second_interval)?;
        Ok(())
    }

    fn verify_two_loop_bridge_loop_cycle(
        &self,
        loop_: LoopHandle,
        required_half_edges: &[HalfEdgeHandle],
    ) -> Result<(), CoreError> {
        let loop_data = self.loops.get(loop_).ok_or_else(|| {
            CoreError::PreconditionFailed(
                "two-loop bridge touched half-edge references a missing loop".into(),
            )
        })?;
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut visited = BTreeSet::new();
        for _ in 0..=self.half_edges.len() {
            let half_edge = self.half_edges.get(cursor).ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "two-loop bridge touched loop references a missing half-edge".into(),
                )
            })?;
            if half_edge.loop_ != loop_ {
                return Err(CoreError::PreconditionFailed(
                    "two-loop bridge touched cycle left its loop".into(),
                ));
            }
            if !visited.insert(cursor) {
                return Err(CoreError::PreconditionFailed(
                    "two-loop bridge touched cycle repeats before closing".into(),
                ));
            }
            cursor = half_edge.next;
            if cursor == start {
                if required_half_edges
                    .iter()
                    .all(|half_edge| visited.contains(half_edge))
                {
                    return Ok(());
                }
                return Err(CoreError::PreconditionFailed(
                    "two-loop bridge touched cycle is missing a required half-edge".into(),
                ));
            }
        }
        Err(CoreError::PreconditionFailed(
            "two-loop bridge touched cycle did not close".into(),
        ))
    }

    fn verify_bridge_cycle_ownership(&self, interval: &[HalfEdgeHandle]) -> Result<(), CoreError> {
        for &half_edge in half_open_path(interval)? {
            let twin = self.half_edges[half_edge].twin;
            let loop_handle = self.half_edges[twin].loop_;
            let face = self.loops[loop_handle].face;
            if face == self.loops[self.half_edges[half_edge].loop_].face {
                return Err(CoreError::PreconditionFailed(
                    "two-loop bridge twin remained on the source face".into(),
                ));
            }
        }
        Ok(())
    }

    fn half_edge_loop(
        &self,
        half_edge: HalfEdgeHandle,
        label: &'static str,
    ) -> Result<LoopHandle, CoreError> {
        self.half_edges
            .get(half_edge)
            .map(|data| data.loop_)
            .ok_or_else(|| CoreError::InvalidHandle(format!("{label} half-edge is invalid")))
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
            "two-loop bridge path needs distinct endpoint origins".into(),
        ));
    }
    Ok(&interval[..interval.len() - 1])
}

fn required_path_handle(
    handle: Option<&HalfEdgeHandle>,
    path: &'static str,
    end: &'static str,
) -> Result<HalfEdgeHandle, CoreError> {
    handle
        .copied()
        .ok_or_else(|| CoreError::PreconditionFailed(format!("{path} bridge path {end} is empty")))
}

#[cfg(test)]
#[path = "loop_bridge_apply_tests.rs"]
mod tests;
