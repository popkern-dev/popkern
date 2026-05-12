//! Boundary-side separator bridge mutation.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::topology::{
    BRep, Edge, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle,
    VertexHandle,
};

use super::boundary_side_separator_bridge_apply_types::{
    CycleSide, SeparatorCounts, SeparatorHalfEdges, required_separator_path_handle,
    separator_half_open_path, wire_separator_path,
};
use super::{BoundarySideSeparatorBridgePlan, FaceFamilyBridgeCycleRequest};

/// Result of a committed boundary-side separator bridge mutation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundarySideSeparatorBridgeResult {
    /// Separator edge allocated between the duplicate bridge consumers.
    pub separator_edge: EdgeHandle,
    /// New bridge faces, one per planned cycle.
    pub bridge_faces: Vec<FaceHandle>,
    /// New bridge loops, one per planned cycle.
    pub bridge_loops: Vec<LoopHandle>,
}

impl BRep {
    /// Transactionally author a bridge family with one duplicate boundary side.
    pub fn boundary_side_separator_bridge(
        &mut self,
        plan: &BoundarySideSeparatorBridgePlan,
    ) -> Result<BoundarySideSeparatorBridgeResult, CoreError> {
        let before = SeparatorCounts::capture(self);
        let mut candidate = self.clone();
        let result = candidate.apply_boundary_side_separator_bridge_candidate(plan)?;
        candidate.verify_boundary_side_separator_postconditions(plan, before, &result)?;
        *self = candidate;
        Ok(result)
    }

    fn apply_boundary_side_separator_bridge_candidate(
        &mut self,
        plan: &BoundarySideSeparatorBridgePlan,
    ) -> Result<BoundarySideSeparatorBridgeResult, CoreError> {
        let duplicate = duplicate_cycle_index(plan)?;
        let separator = self.insert_separator_edge(plan.separator_source_half_edge)?;
        let mut bridge_faces = Vec::with_capacity(plan.cycles.len());
        let mut bridge_loops = Vec::with_capacity(plan.cycles.len());

        for (index, cycle) in plan.cycles.iter().enumerate() {
            let loop_ = self.insert_separator_bridge_face_loop(cycle.first_loop)?;
            let first_path = self.separator_bridge_path(
                plan,
                separator,
                cycle,
                CycleSide::First,
                index == duplicate,
            )?;
            let second_path = self.separator_bridge_path(
                plan,
                separator,
                cycle,
                CycleSide::Second,
                index == duplicate,
            )?;
            self.splice_separator_bridge_cycle(loop_, &first_path, &second_path)?;
            bridge_faces.push(self.loops[loop_].face);
            bridge_loops.push(loop_);
        }
        self.finalize_separator_edge(separator)?;

        Ok(BoundarySideSeparatorBridgeResult {
            separator_edge: separator.edge,
            bridge_faces,
            bridge_loops,
        })
    }

    fn insert_separator_edge(
        &mut self,
        source: HalfEdgeHandle,
    ) -> Result<SeparatorHalfEdges, CoreError> {
        let source_data = self.half_edges.get(source).ok_or_else(|| {
            CoreError::InvalidHandle("boundary-side separator source half-edge is stale".into())
        })?;
        let start = source_data.origin;
        let end = self
            .half_edges
            .get(source_data.next)
            .map(|next| next.origin)
            .ok_or_else(|| {
                CoreError::InvalidHandle(
                    "boundary-side separator source after-last half-edge is stale".into(),
                )
            })?;
        let edge = self.edges.insert(Edge {
            half_edge: Default::default(),
            curve: self.edges[source_data.edge].curve.clone(),
            name: self.edges[source_data.edge].name,
            tolerance: self.edges[source_data.edge].tolerance,
            canonical: self.edges[source_data.edge].canonical,
            role: self.edges[source_data.edge].role,
        });
        let forward = self.insert_separator_half_edge(edge, start);
        let reverse = self.insert_separator_half_edge(edge, end);
        self.half_edges[forward].twin = reverse;
        self.half_edges[reverse].twin = forward;
        self.edges[edge].half_edge = forward;
        Ok(SeparatorHalfEdges {
            edge,
            forward,
            reverse,
        })
    }

    fn insert_separator_half_edge(
        &mut self,
        edge: EdgeHandle,
        origin: VertexHandle,
    ) -> HalfEdgeHandle {
        self.half_edges.insert(HalfEdge {
            origin,
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge,
            loop_: Default::default(),
            pcurve: None,
        })
    }

    fn separator_bridge_path(
        &mut self,
        plan: &BoundarySideSeparatorBridgePlan,
        separator: SeparatorHalfEdges,
        cycle: &FaceFamilyBridgeCycleRequest,
        side: CycleSide,
        duplicate_cycle: bool,
    ) -> Result<Vec<HalfEdgeHandle>, CoreError> {
        let interval = match side {
            CycleSide::First => &cycle.first_interval,
            CycleSide::Second => &cycle.second_interval,
        };
        if duplicate_cycle && interval_uses_separator_source(plan, interval)? {
            return match side {
                CycleSide::First => self.create_separator_bridge_boundary_twins(interval),
                CycleSide::Second => self.create_separator_bridge_path_with_substitution(
                    plan,
                    interval,
                    separator.forward,
                ),
            };
        }
        if !duplicate_cycle && interval_uses_separator_source(plan, interval)? {
            return self.create_separator_bridge_path_with_substitution(
                plan,
                interval,
                separator.reverse,
            );
        }
        self.create_separator_bridge_boundary_twins(interval)
    }

    fn create_separator_bridge_path_with_substitution(
        &mut self,
        plan: &BoundarySideSeparatorBridgePlan,
        interval: &[HalfEdgeHandle],
        substitute: HalfEdgeHandle,
    ) -> Result<Vec<HalfEdgeHandle>, CoreError> {
        let mut path = Vec::new();
        let mut substituted = false;
        for &half_edge in separator_half_open_path(interval)? {
            if half_edge == plan.separator_source_half_edge {
                path.push(substitute);
                substituted = true;
            } else {
                let mut twin_path = self.create_separator_bridge_boundary_twins(&[
                    half_edge,
                    self.half_edges[half_edge].next,
                ])?;
                path.append(&mut twin_path);
            }
        }
        if !substituted {
            return Err(CoreError::PreconditionFailed(
                "boundary-side separator substitution did not consume the source side".into(),
            ));
        }
        path.reverse();
        Ok(path)
    }

    fn create_separator_bridge_boundary_twins(
        &mut self,
        interval: &[HalfEdgeHandle],
    ) -> Result<Vec<HalfEdgeHandle>, CoreError> {
        let mut twins = Vec::new();
        for &half_edge in separator_half_open_path(interval)? {
            let data = &self.half_edges[half_edge];
            let edge = data.edge;
            let origin = self.half_edges[data.next].origin;
            if !self.boundary_edges.contains(&edge) || self.half_edges.contains_key(data.twin) {
                return Err(CoreError::PreconditionFailed(format!(
                    "boundary-side separator bridge path half-edge {half_edge:?} is not an unowned boundary side"
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

    fn insert_separator_bridge_face_loop(
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

    fn splice_separator_bridge_cycle(
        &mut self,
        bridge_loop: LoopHandle,
        first_path: &[HalfEdgeHandle],
        second_path: &[HalfEdgeHandle],
    ) -> Result<(), CoreError> {
        let first_start = required_separator_path_handle(first_path.first(), "first", "start")?;
        let first_end = required_separator_path_handle(first_path.last(), "first", "end")?;
        let second_start = required_separator_path_handle(second_path.first(), "second", "start")?;
        let second_end = required_separator_path_handle(second_path.last(), "second", "end")?;
        let first_connector = self.insert_separator_connector(
            self.separator_path_terminal_vertex(first_end)?,
            bridge_loop,
        );
        let second_connector = self.insert_separator_connector(
            self.separator_path_terminal_vertex(second_end)?,
            bridge_loop,
        );
        wire_separator_path(self, first_path);
        wire_separator_path(self, second_path);
        self.half_edges[first_end].next = first_connector;
        self.half_edges[first_connector].prev = first_end;
        self.half_edges[first_connector].next = second_start;
        self.half_edges[second_start].prev = first_connector;
        self.half_edges[second_end].next = second_connector;
        self.half_edges[second_connector].prev = second_end;
        self.half_edges[second_connector].next = first_start;
        self.half_edges[first_start].prev = second_connector;
        for half_edge in first_path
            .iter()
            .chain(second_path.iter())
            .chain([first_connector, second_connector].iter())
        {
            self.half_edges[*half_edge].loop_ = bridge_loop;
        }
        self.loops[bridge_loop].half_edge = first_start;
        Ok(())
    }

    fn insert_separator_connector(
        &mut self,
        origin: VertexHandle,
        loop_: LoopHandle,
    ) -> HalfEdgeHandle {
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
        half_edge
    }

    fn finalize_separator_edge(&mut self, separator: SeparatorHalfEdges) -> Result<(), CoreError> {
        let forward_used = self
            .half_edges
            .get(separator.forward)
            .is_some_and(|half_edge| half_edge.loop_ != LoopHandle::default());
        let reverse_used = self
            .half_edges
            .get(separator.reverse)
            .is_some_and(|half_edge| half_edge.loop_ != LoopHandle::default());
        match (forward_used, reverse_used) {
            (true, true) => Ok(()),
            (true, false) => {
                self.half_edges[separator.forward].twin = HalfEdgeHandle::default();
                self.edges[separator.edge].half_edge = separator.forward;
                self.half_edges.remove(separator.reverse);
                self.boundary_edges.insert(separator.edge);
                Ok(())
            }
            (false, true) => {
                self.half_edges[separator.reverse].twin = HalfEdgeHandle::default();
                self.edges[separator.edge].half_edge = separator.reverse;
                self.half_edges.remove(separator.forward);
                self.boundary_edges.insert(separator.edge);
                Ok(())
            }
            (false, false) => Err(CoreError::PreconditionFailed(
                "boundary-side separator edge was not consumed by any bridge cycle".into(),
            )),
        }
    }

    fn separator_path_terminal_vertex(
        &self,
        half_edge: HalfEdgeHandle,
    ) -> Result<VertexHandle, CoreError> {
        let data = self.half_edges.get(half_edge).ok_or_else(|| {
            CoreError::InvalidHandle("separator bridge path half-edge is stale".into())
        })?;
        if self.boundary_edges.contains(&data.edge) {
            return self
                .half_edges
                .get(data.next)
                .map(|next| next.origin)
                .ok_or_else(|| {
                    CoreError::InvalidHandle(
                        "separator bridge boundary path has stale next half-edge".into(),
                    )
                });
        }
        self.half_edges
            .get(data.twin)
            .map(|twin| twin.origin)
            .ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "separator bridge path half-edge has no terminal twin".into(),
                )
            })
    }

    fn verify_boundary_side_separator_postconditions(
        &self,
        plan: &BoundarySideSeparatorBridgePlan,
        before: SeparatorCounts,
        result: &BoundarySideSeparatorBridgeResult,
    ) -> Result<(), CoreError> {
        let expected = SeparatorCounts {
            vertices: plan.expected_delta.vertices,
            edges: plan.expected_delta.edges,
            faces: plan.expected_delta.faces,
            loops: plan.expected_delta.loops,
        };
        let actual = SeparatorCounts::capture(self).delta_from(before);
        if actual != expected {
            return Err(CoreError::PreconditionFailed(format!(
                "boundary-side separator bridge postcondition failed: expected_delta=dV={} dE={} dF={} dL={} actual_delta=dV={} dE={} dF={} dL={}",
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
        self.verify_separator_edge_finalized(result.separator_edge)?;
        self.validate_half_edge_pairing().map_err(CoreError::from)?;
        self.validate_half_edge_cycles().map_err(CoreError::from)?;
        self.validate_manifold_edges().map_err(CoreError::from)?;
        for &loop_ in &result.bridge_loops {
            self.verify_separator_bridge_loop(loop_)?;
        }
        Ok(())
    }

    fn verify_separator_edge_finalized(&self, edge: EdgeHandle) -> Result<(), CoreError> {
        let half_edge = self
            .edges
            .get(edge)
            .map(|edge_data| edge_data.half_edge)
            .ok_or_else(|| {
                CoreError::PreconditionFailed("boundary-side separator edge is missing".into())
            })?;
        let data = self.half_edges.get(half_edge).ok_or_else(|| {
            CoreError::PreconditionFailed(
                "boundary-side separator edge references a missing half-edge".into(),
            )
        })?;
        if self.boundary_edges.contains(&edge) {
            if data.twin == HalfEdgeHandle::default() {
                return Ok(());
            }
            return Err(CoreError::PreconditionFailed(
                "boundary-side separator boundary edge kept a twin".into(),
            ));
        }
        let twin = self.half_edges.get(data.twin).ok_or_else(|| {
            CoreError::PreconditionFailed(
                "boundary-side separator paired edge has missing twin".into(),
            )
        })?;
        if twin.twin == half_edge && twin.edge == edge {
            return Ok(());
        }
        Err(CoreError::PreconditionFailed(
            "boundary-side separator paired edge is inconsistent".into(),
        ))
    }

    fn verify_separator_bridge_loop(&self, loop_: LoopHandle) -> Result<(), CoreError> {
        let start = self.loops[loop_].half_edge;
        let mut cursor = start;
        let mut visited = BTreeSet::new();
        for _ in 0..=self.half_edges.len() {
            let half_edge = self.half_edges.get(cursor).ok_or_else(|| {
                CoreError::PreconditionFailed(
                    "boundary-side separator bridge loop has stale half-edge".into(),
                )
            })?;
            if half_edge.loop_ != loop_ {
                return Err(CoreError::PreconditionFailed(
                    "boundary-side separator bridge cycle left its loop".into(),
                ));
            }
            if !visited.insert(cursor) {
                return Err(CoreError::PreconditionFailed(
                    "boundary-side separator bridge cycle repeats before closing".into(),
                ));
            }
            cursor = half_edge.next;
            if cursor == start {
                return Ok(());
            }
        }
        Err(CoreError::PreconditionFailed(
            "boundary-side separator bridge cycle did not close".into(),
        ))
    }
}

fn duplicate_cycle_index(plan: &BoundarySideSeparatorBridgePlan) -> Result<usize, CoreError> {
    let matches = plan
        .cycles
        .iter()
        .enumerate()
        .filter(|(_, cycle)| {
            cycle.first_loop == cycle.second_loop
                && interval_uses_separator_source(plan, &cycle.first_interval).unwrap_or(false)
                && interval_uses_separator_source(plan, &cycle.second_interval).unwrap_or(false)
        })
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    if matches.len() == 1 {
        return Ok(matches[0]);
    }
    Err(CoreError::PreconditionFailed(format!(
        "boundary-side separator bridge expected one duplicate same-loop cycle, got {}",
        matches.len()
    )))
}

fn interval_uses_separator_source(
    plan: &BoundarySideSeparatorBridgePlan,
    interval: &[HalfEdgeHandle],
) -> Result<bool, CoreError> {
    Ok(separator_half_open_path(interval)?.contains(&plan.separator_source_half_edge))
}
