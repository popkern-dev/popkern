//! Two-loop bridge preflight for cross-face topology operations.

use std::collections::BTreeSet;

use crate::error::CoreError;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

/// Operand set for a two-loop bridge between existing loop intervals.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TwoLoopBridgeRequest {
    /// First loop to bridge.
    pub first_loop: LoopHandle,
    /// Half-edge whose origin is the first endpoint on `first_loop`.
    pub first_start: HalfEdgeHandle,
    /// Half-edge whose origin is the second endpoint on `first_loop`.
    pub first_end: HalfEdgeHandle,
    /// Second loop to bridge.
    pub second_loop: LoopHandle,
    /// Half-edge whose origin is the first endpoint on `second_loop`.
    pub second_start: HalfEdgeHandle,
    /// Half-edge whose origin is the second endpoint on `second_loop`.
    pub second_end: HalfEdgeHandle,
}

/// Validated read-only plan for a two-loop bridge.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TwoLoopBridgePlan {
    /// First face owning `first_loop`.
    pub first_face: FaceHandle,
    /// Second face owning `second_loop`.
    pub second_face: FaceHandle,
    /// Length of the first loop.
    pub first_loop_len: usize,
    /// Length of the second loop.
    pub second_loop_len: usize,
    /// First-loop interval from `first_start` to `first_end`, inclusive.
    pub first_forward_interval: Vec<HalfEdgeHandle>,
    /// First-loop interval from `first_end` to `first_start`, inclusive.
    pub first_reverse_interval: Vec<HalfEdgeHandle>,
    /// Second-loop interval from `second_start` to `second_end`, inclusive.
    pub second_forward_interval: Vec<HalfEdgeHandle>,
    /// Second-loop interval from `second_end` to `second_start`, inclusive.
    pub second_reverse_interval: Vec<HalfEdgeHandle>,
    /// The two complementary boundary cycles implied by the bridge.
    pub boundary_cycles: [TwoLoopBridgeBoundaryCycle; 2],
    /// Expected topology delta for the mutating bridge operation.
    pub expected_delta: TwoLoopBridgeTopologyDelta,
}

/// One closed boundary cycle implied by pairing opposite intervals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TwoLoopBridgeBoundaryCycle {
    /// Interval from the first loop.
    pub first_interval: Vec<HalfEdgeHandle>,
    /// Opposite interval from the second loop.
    pub second_interval: Vec<HalfEdgeHandle>,
}

/// Expected topology count delta for the future two-loop bridge mutation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TwoLoopBridgeTopologyDelta {
    /// Vertex count delta.
    pub vertices: i64,
    /// Edge count delta.
    pub edges: i64,
    /// Face count delta.
    pub faces: i64,
    /// Loop count delta.
    pub loops: i64,
}

impl BRep {
    /// Validate two existing loops as operands for a future cross-face bridge.
    ///
    /// This is intentionally read-only. The mutating bridge operator must consume
    /// this complete plan and preserve its interval coverage postconditions.
    ///
    /// # Errors
    ///
    /// - [`CoreError::InvalidHandle`] if any loop or half-edge is stale.
    /// - [`CoreError::PreconditionFailed`] if endpoint half-edges do not belong
    ///   to the stated loops, the loops are not distinct, or either directed
    ///   interval pair is not complementary coverage of its loop.
    pub fn preflight_two_loop_bridge(
        &self,
        request: TwoLoopBridgeRequest,
    ) -> Result<TwoLoopBridgePlan, CoreError> {
        if request.first_loop == request.second_loop {
            return Err(CoreError::PreconditionFailed(
                "two-loop bridge requires distinct loops".into(),
            ));
        }
        let first_face = self.loop_face(request.first_loop)?;
        let second_face = self.loop_face(request.second_loop)?;
        if first_face == second_face {
            return Err(CoreError::PreconditionFailed(
                "two-loop bridge is for loops on distinct faces; use mekr for loops on one face"
                    .into(),
            ));
        }

        let first_loop_len = self.loop_half_edge_count(request.first_loop, request.first_start)?;
        let second_loop_len =
            self.loop_half_edge_count(request.second_loop, request.second_start)?;
        self.require_disjoint_bridge_loop_edges(request.first_loop, request.second_loop)?;
        let first_forward_interval = self.loop_interval_between_origins(
            request.first_loop,
            request.first_start,
            request.first_end,
        )?;
        let first_reverse_interval = self.loop_interval_between_origins(
            request.first_loop,
            request.first_end,
            request.first_start,
        )?;
        let second_forward_interval = self.loop_interval_between_origins(
            request.second_loop,
            request.second_start,
            request.second_end,
        )?;
        let second_reverse_interval = self.loop_interval_between_origins(
            request.second_loop,
            request.second_end,
            request.second_start,
        )?;
        require_complementary_intervals(
            first_loop_len,
            first_forward_interval.len(),
            first_reverse_interval.len(),
            "first",
        )?;
        require_complementary_intervals(
            second_loop_len,
            second_forward_interval.len(),
            second_reverse_interval.len(),
            "second",
        )?;

        let boundary_cycles = [
            TwoLoopBridgeBoundaryCycle {
                first_interval: first_forward_interval.clone(),
                second_interval: second_reverse_interval.clone(),
            },
            TwoLoopBridgeBoundaryCycle {
                first_interval: first_reverse_interval.clone(),
                second_interval: second_forward_interval.clone(),
            },
        ];
        Ok(TwoLoopBridgePlan {
            first_face,
            second_face,
            first_loop_len,
            second_loop_len,
            first_forward_interval,
            first_reverse_interval,
            second_forward_interval,
            second_reverse_interval,
            boundary_cycles,
            expected_delta: TwoLoopBridgeTopologyDelta {
                vertices: 0,
                edges: 2,
                faces: 2,
                loops: 2,
            },
        })
    }

    fn loop_face(&self, loop_handle: LoopHandle) -> Result<FaceHandle, CoreError> {
        self.loops
            .get(loop_handle)
            .map(|loop_data| loop_data.face)
            .ok_or_else(|| CoreError::InvalidHandle("loop handle is invalid".into()))
    }

    fn loop_half_edge_count(
        &self,
        loop_handle: LoopHandle,
        start: HalfEdgeHandle,
    ) -> Result<usize, CoreError> {
        let mut cursor = start;
        let mut count = 0;
        for _ in 0..=self.half_edges.len() {
            let half_edge = self
                .half_edges
                .get(cursor)
                .ok_or_else(|| CoreError::InvalidHandle("half-edge handle is invalid".into()))?;
            if half_edge.loop_ != loop_handle {
                return Err(CoreError::PreconditionFailed(
                    "half-edge does not belong to requested loop".into(),
                ));
            }
            count += 1;
            cursor = half_edge.next;
            if cursor == start {
                return Ok(count);
            }
        }
        Err(CoreError::PreconditionFailed(
            "loop walk did not close while counting half-edges".into(),
        ))
    }

    fn loop_interval_between_origins(
        &self,
        loop_handle: LoopHandle,
        start: HalfEdgeHandle,
        end: HalfEdgeHandle,
    ) -> Result<Vec<HalfEdgeHandle>, CoreError> {
        let mut cursor = start;
        let mut interval = Vec::new();
        for _ in 0..=self.half_edges.len() {
            let half_edge = self
                .half_edges
                .get(cursor)
                .ok_or_else(|| CoreError::InvalidHandle("half-edge handle is invalid".into()))?;
            if half_edge.loop_ != loop_handle {
                return Err(CoreError::PreconditionFailed(
                    "half-edge interval left requested loop".into(),
                ));
            }
            interval.push(cursor);
            if cursor == end {
                return Ok(interval);
            }
            cursor = half_edge.next;
            if cursor == start {
                break;
            }
        }
        Err(CoreError::PreconditionFailed(
            "loop interval endpoints are not reachable in requested loop".into(),
        ))
    }

    fn require_disjoint_bridge_loop_edges(
        &self,
        first_loop: LoopHandle,
        second_loop: LoopHandle,
    ) -> Result<(), CoreError> {
        let first_edges = self.loop_edge_set(first_loop)?;
        let second_edges = self.loop_edge_set(second_loop)?;
        let shared = first_edges
            .intersection(&second_edges)
            .copied()
            .collect::<Vec<_>>();
        if shared.is_empty() {
            return Ok(());
        }
        Err(CoreError::PreconditionFailed(format!(
            "two-loop bridge operands must not share existing edges: shared_edges={shared:?}"
        )))
    }

    fn loop_edge_set(
        &self,
        loop_handle: LoopHandle,
    ) -> Result<BTreeSet<crate::topology::EdgeHandle>, CoreError> {
        let start = self
            .loops
            .get(loop_handle)
            .ok_or_else(|| CoreError::InvalidHandle("loop handle is invalid".into()))?
            .half_edge;
        let mut cursor = start;
        let mut edges = BTreeSet::new();
        for _ in 0..=self.half_edges.len() {
            let half_edge = self
                .half_edges
                .get(cursor)
                .ok_or_else(|| CoreError::InvalidHandle("half-edge handle is invalid".into()))?;
            if half_edge.loop_ != loop_handle {
                return Err(CoreError::PreconditionFailed(
                    "loop edge-set walk left requested loop".into(),
                ));
            }
            edges.insert(half_edge.edge);
            cursor = half_edge.next;
            if cursor == start {
                return Ok(edges);
            }
        }
        Err(CoreError::PreconditionFailed(
            "loop edge-set walk did not close".into(),
        ))
    }
}

fn require_complementary_intervals(
    loop_len: usize,
    forward_len: usize,
    reverse_len: usize,
    label: &'static str,
) -> Result<(), CoreError> {
    if forward_len + reverse_len == loop_len + 2 {
        return Ok(());
    }
    Err(CoreError::PreconditionFailed(format!(
        "{label} loop bridge intervals are not complementary: loop_len={loop_len} \
         forward_len={forward_len} reverse_len={reverse_len}"
    )))
}

#[cfg(test)]
mod tests {
    use crate::euler::TwoLoopBridgeRequest;
    use crate::tolerance::Tolerance;
    use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

    #[test]
    fn preflight_two_loop_bridge_accepts_distinct_cube_face_loops() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let (first_face, second_face) = disjoint_cube_face_pair(&brep);
        let first_loop = brep.faces[first_face].outer_loop;
        let second_loop = brep.faces[second_face].outer_loop;
        let first_start = brep.loops[first_loop].half_edge;
        let first_end = brep.half_edges[first_start].next;
        let second_start = brep.loops[second_loop].half_edge;
        let second_end = brep.half_edges[second_start].next;

        let plan = brep
            .preflight_two_loop_bridge(TwoLoopBridgeRequest {
                first_loop,
                first_start,
                first_end,
                second_loop,
                second_start,
                second_end,
            })
            .expect("two cube loops can be preflighted as bridge operands");

        assert_eq!(plan.first_face, first_face);
        assert_eq!(plan.second_face, second_face);
        assert_eq!(
            plan.first_forward_interval.len() + plan.first_reverse_interval.len(),
            plan.first_loop_len + 2
        );
        assert_eq!(
            plan.second_forward_interval.len() + plan.second_reverse_interval.len(),
            plan.second_loop_len + 2
        );
        assert_eq!(plan.boundary_cycles.len(), 2);
        assert_eq!(plan.expected_delta.vertices, 0);
        assert_eq!(plan.expected_delta.edges, 2);
        assert_eq!(plan.expected_delta.faces, 2);
        assert_eq!(plan.expected_delta.loops, 2);
    }

    #[test]
    fn preflight_two_loop_bridge_refuses_adjacent_faces_with_shared_edge() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let (first_face, second_face) = adjacent_cube_face_pair(&brep);
        let first_loop = brep.faces[first_face].outer_loop;
        let second_loop = brep.faces[second_face].outer_loop;
        let first_start = brep.loops[first_loop].half_edge;
        let first_end = brep.half_edges[first_start].next;
        let second_start = brep.loops[second_loop].half_edge;
        let second_end = brep.half_edges[second_start].next;

        let err = brep
            .preflight_two_loop_bridge(TwoLoopBridgeRequest {
                first_loop,
                first_start,
                first_end,
                second_loop,
                second_start,
                second_end,
            })
            .expect_err("adjacent cube faces share an existing edge");

        assert!(err.to_string().contains("must not share existing edges"));
    }

    #[test]
    fn preflight_two_loop_bridge_refuses_same_loop() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let face = brep.faces.keys().next().expect("cube has a face");
        let loop_handle = brep.faces[face].outer_loop;
        let start = brep.loops[loop_handle].half_edge;
        let end = brep.half_edges[start].next;

        let result = brep.preflight_two_loop_bridge(TwoLoopBridgeRequest {
            first_loop: loop_handle,
            first_start: start,
            first_end: end,
            second_loop: loop_handle,
            second_start: start,
            second_end: end,
        });

        assert!(result.is_err());
    }

    #[test]
    fn preflight_two_loop_bridge_refuses_half_edge_outside_loop() {
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
        let faces = brep.faces.keys().collect::<Vec<_>>();
        let first_loop = brep.faces[faces[0]].outer_loop;
        let second_loop = brep.faces[faces[1]].outer_loop;
        let first_start = brep.loops[first_loop].half_edge;
        let second_start = brep.loops[second_loop].half_edge;
        let second_end = brep.half_edges[second_start].next;

        let result = brep.preflight_two_loop_bridge(TwoLoopBridgeRequest {
            first_loop,
            first_start,
            first_end: second_start,
            second_loop,
            second_start,
            second_end,
        });

        assert!(result.is_err());
    }

    #[allow(dead_code)]
    fn _assert_handles_are_public(_: FaceHandle, _: LoopHandle, _: HalfEdgeHandle) {}

    fn disjoint_cube_face_pair(brep: &BRep) -> (FaceHandle, FaceHandle) {
        cube_face_pair_by_edge_overlap(brep, false)
    }

    fn adjacent_cube_face_pair(brep: &BRep) -> (FaceHandle, FaceHandle) {
        cube_face_pair_by_edge_overlap(brep, true)
    }

    fn cube_face_pair_by_edge_overlap(
        brep: &BRep,
        want_shared_edge: bool,
    ) -> (FaceHandle, FaceHandle) {
        let faces = brep.faces.keys().collect::<Vec<_>>();
        for (i, first) in faces.iter().enumerate() {
            for second in faces.iter().skip(i + 1) {
                if loop_edges_share_any(brep, *first, *second) == want_shared_edge {
                    return (*first, *second);
                }
            }
        }
        panic!("cube fixture does not contain requested face pair")
    }

    fn loop_edges_share_any(brep: &BRep, first: FaceHandle, second: FaceHandle) -> bool {
        let first_loop = brep.faces[first].outer_loop;
        let second_loop = brep.faces[second].outer_loop;
        let mut first_edges = std::collections::BTreeSet::new();
        let mut cursor = brep.loops[first_loop].half_edge;
        loop {
            first_edges.insert(brep.half_edges[cursor].edge);
            cursor = brep.half_edges[cursor].next;
            if cursor == brep.loops[first_loop].half_edge {
                break;
            }
        }
        let mut cursor = brep.loops[second_loop].half_edge;
        loop {
            if first_edges.contains(&brep.half_edges[cursor].edge) {
                return true;
            }
            cursor = brep.half_edges[cursor].next;
            if cursor == brep.loops[second_loop].half_edge {
                break;
            }
        }
        false
    }
}
