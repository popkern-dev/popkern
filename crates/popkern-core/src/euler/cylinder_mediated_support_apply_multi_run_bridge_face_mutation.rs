//! Mutation for multi-run bridge-face plans.

use std::collections::{BTreeSet, HashSet};

use crate::error::CoreError;
use crate::topology::{
    BRep, Edge, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle,
    VertexHandle,
};

use super::super::SplitResourceMergeOperand;
use super::{MultiRunBridgeFacePlan, multi_run_bridge_face_plan};

#[derive(Debug, Clone, PartialEq, Eq)]
pub(in crate::euler::cylinder_mediated_support_apply) struct MultiRunBridgeFaceMutationResult {
    pub(in crate::euler::cylinder_mediated_support_apply) bridge_faces: Vec<FaceHandle>,
    pub(in crate::euler::cylinder_mediated_support_apply) bridge_loops: Vec<LoopHandle>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct MultiRunBridgeFaceCounts {
    edges: i64,
    faces: i64,
    loops: i64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct MultiRunBridgePath {
    half_edges: Vec<HalfEdgeHandle>,
    terminal: VertexHandle,
    separator_edges: usize,
}

pub(super) fn apply_multi_run_bridge_face_plans(
    brep: &mut BRep,
    operands: &[SplitResourceMergeOperand],
) -> Result<Option<MultiRunBridgeFaceMutationResult>, CoreError> {
    let plans = operands
        .iter()
        .filter_map(multi_run_bridge_face_plan)
        .collect::<Vec<_>>();
    if plans.is_empty() {
        return Ok(None);
    }
    let before = MultiRunBridgeFaceCounts::capture(brep);
    let mut expected_edges = 0i64;
    let mut bridge_faces = Vec::with_capacity(plans.len());
    let mut bridge_loops = Vec::with_capacity(plans.len());
    for (operand, plan) in operands
        .iter()
        .filter_map(|operand| multi_run_bridge_face_plan(operand).map(|plan| (operand, plan)))
    {
        let loop_ = insert_multi_run_bridge_face_loop(brep, operand.outside_loop)?;
        let paths = multi_run_bridge_paths(brep, operand, &plan)?;
        expected_edges += paths.len() as i64
            + paths
                .iter()
                .map(|path| path.separator_edges as i64)
                .sum::<i64>();
        splice_multi_run_bridge_cycle(brep, loop_, &paths)?;
        bridge_faces.push(brep.loops[loop_].face);
        bridge_loops.push(loop_);
    }
    verify_multi_run_bridge_face_postconditions(
        brep,
        before,
        expected_edges,
        &plans,
        &bridge_loops,
    )?;
    Ok(Some(MultiRunBridgeFaceMutationResult {
        bridge_faces,
        bridge_loops,
    }))
}

fn multi_run_bridge_paths(
    brep: &mut BRep,
    operand: &SplitResourceMergeOperand,
    plan: &MultiRunBridgeFacePlan,
) -> Result<Vec<MultiRunBridgePath>, CoreError> {
    let mut paths = Vec::with_capacity(plan.cylinder_runs + 1);
    let mut consumed_edges = HashSet::new();
    paths.push(create_multi_run_bridge_boundary_twins(
        brep,
        half_open_path(&operand.outside_interval)?,
        &mut consumed_edges,
    )?);
    for &run_index in &plan.ordered_run_indices {
        paths.push(create_multi_run_bridge_boundary_twins(
            brep,
            &operand.cylinder_runs[run_index].half_edges,
            &mut consumed_edges,
        )?);
    }
    Ok(paths)
}

fn create_multi_run_bridge_boundary_twins(
    brep: &mut BRep,
    half_edges: &[HalfEdgeHandle],
    consumed_edges: &mut HashSet<EdgeHandle>,
) -> Result<MultiRunBridgePath, CoreError> {
    // Owner exposure may have CONSUMED some path half-edges outright
    // (retiring a slit fragment removes both of its half-edges and
    // re-stitches the loop — PK-092 paired cross-owner batch). The
    // path is built over the surviving chain.
    let half_edges: Vec<HalfEdgeHandle> = half_edges
        .iter()
        .copied()
        .filter(|half_edge| brep.half_edges.contains_key(*half_edge))
        .collect();
    if half_edges.is_empty() {
        return Err(CoreError::PreconditionFailed(
            "multi-run bridge face path is empty".into(),
        ));
    }
    let terminal = brep
        .half_edges
        .get(half_edges[0])
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| {
            CoreError::InvalidHandle("multi-run bridge path first half-edge is stale".into())
        })?;
    let mut twins = Vec::with_capacity(half_edges.len());
    let mut separator_edges = 0usize;
    for &half_edge in &half_edges {
        let data = brep
            .half_edges
            .get(half_edge)
            .ok_or_else(|| {
                CoreError::InvalidHandle("multi-run bridge path half-edge is stale".into())
            })?
            .clone();
        let edge = data.edge;
        let origin = brep
            .half_edges
            .get(data.next)
            .map(|next| next.origin)
            .ok_or_else(|| {
                CoreError::InvalidHandle(
                    "multi-run bridge path after-last half-edge is stale".into(),
                )
            })?;
        if !consumed_edges.insert(edge) || brep.half_edges.contains_key(data.twin) {
            twins.push(insert_multi_run_bridge_separator(brep, half_edge, origin)?);
            separator_edges += 1;
            continue;
        }
        if !brep.boundary_edges.contains(&edge) {
            return Err(CoreError::PreconditionFailed(format!(
                "multi-run bridge path half-edge {half_edge:?} is not an unowned boundary side"
            )));
        }
        let twin = brep.half_edges.insert(HalfEdge {
            origin,
            twin: half_edge,
            next: Default::default(),
            prev: Default::default(),
            edge,
            loop_: Default::default(),
            pcurve: None,
        });
        brep.half_edges[half_edge].twin = twin;
        brep.boundary_edges.remove(&edge);
        twins.push(twin);
    }
    twins.reverse();
    Ok(MultiRunBridgePath {
        half_edges: twins,
        terminal,
        separator_edges,
    })
}

fn insert_multi_run_bridge_separator(
    brep: &mut BRep,
    source: HalfEdgeHandle,
    origin: VertexHandle,
) -> Result<HalfEdgeHandle, CoreError> {
    let (source_edge_handle, source_pcurve) = {
        let source_data = brep.half_edges.get(source).ok_or_else(|| {
            CoreError::InvalidHandle("multi-run bridge separator source is stale".into())
        })?;
        (source_data.edge, source_data.pcurve.clone())
    };
    let (curve, name, tolerance, canonical, role) = {
        let source_edge = brep.edges.get(source_edge_handle).ok_or_else(|| {
            CoreError::InvalidHandle("multi-run bridge separator source edge is stale".into())
        })?;
        (
            source_edge.curve.clone(),
            source_edge.name,
            source_edge.tolerance,
            source_edge.canonical,
            source_edge.role,
        )
    };
    let edge = brep.insert_edge_tracked(Edge {
        half_edge: Default::default(),
        curve,
        name,
        tolerance,
        canonical,
        role,
    });
    let half_edge = brep.half_edges.insert(HalfEdge {
        origin,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge,
        loop_: Default::default(),
        pcurve: source_pcurve,
    });
    brep.edges[edge].half_edge = half_edge;
    brep.boundary_edges.insert(edge);
    Ok(half_edge)
}

fn insert_multi_run_bridge_face_loop(
    brep: &mut BRep,
    source_loop: LoopHandle,
) -> Result<LoopHandle, CoreError> {
    let source_face = brep
        .loops
        .get(source_loop)
        .map(|loop_| loop_.face)
        .ok_or_else(|| CoreError::InvalidHandle("multi-run bridge source loop is stale".into()))?;
    let source = brep
        .faces
        .get(source_face)
        .ok_or_else(|| CoreError::InvalidHandle("multi-run bridge source face is stale".into()))?
        .clone();
    let bridge_face = brep.faces.insert(Face {
        outer_loop: Default::default(),
        inner_loops: vec![],
        shell: source.shell,
        surface: None,
        name: None,
        tolerance: source.tolerance,
        role: Default::default(),
    });
    let bridge_loop = brep.loops.insert(Loop {
        half_edge: Default::default(),
        face: bridge_face,
    });
    brep.faces[bridge_face].outer_loop = bridge_loop;
    Ok(bridge_loop)
}

fn splice_multi_run_bridge_cycle(
    brep: &mut BRep,
    bridge_loop: LoopHandle,
    paths: &[MultiRunBridgePath],
) -> Result<(), CoreError> {
    if paths.len() < 3 {
        return Err(CoreError::PreconditionFailed(
            "multi-run bridge face requires an outside path and at least two cylinder runs".into(),
        ));
    }
    for path in paths {
        wire_multi_run_bridge_path(brep, &path.half_edges)?;
    }
    let connectors = paths
        .iter()
        .map(|path| {
            required_path_handle(path.half_edges.last())?;
            Ok(insert_multi_run_bridge_connector(
                brep,
                path.terminal,
                bridge_loop,
            ))
        })
        .collect::<Result<Vec<_>, CoreError>>()?;
    for index in 0..paths.len() {
        let current_end = required_path_handle(paths[index].half_edges.last())?;
        let next_start = required_path_handle(paths[(index + 1) % paths.len()].half_edges.first())?;
        let connector = connectors[index];
        brep.half_edges[current_end].next = connector;
        brep.half_edges[connector].prev = current_end;
        brep.half_edges[connector].next = next_start;
        brep.half_edges[next_start].prev = connector;
    }
    let first_start = required_path_handle(paths[0].half_edges.first())?;
    for half_edge in paths
        .iter()
        .flat_map(|path| path.half_edges.iter())
        .copied()
        .chain(connectors.iter().copied())
    {
        brep.half_edges[half_edge].loop_ = bridge_loop;
    }
    brep.loops[bridge_loop].half_edge = first_start;
    Ok(())
}

fn insert_multi_run_bridge_connector(
    brep: &mut BRep,
    origin: VertexHandle,
    loop_: LoopHandle,
) -> HalfEdgeHandle {
    let edge = brep.insert_edge_tracked(Edge {
        half_edge: Default::default(),
        curve: None,
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: Default::default(),
    });
    let half_edge = brep.half_edges.insert(HalfEdge {
        origin,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge,
        loop_,
        pcurve: None,
    });
    brep.edges[edge].half_edge = half_edge;
    brep.boundary_edges.insert(edge);
    half_edge
}

fn wire_multi_run_bridge_path(brep: &mut BRep, path: &[HalfEdgeHandle]) -> Result<(), CoreError> {
    if path.is_empty() {
        return Err(CoreError::PreconditionFailed(
            "multi-run bridge path is empty".into(),
        ));
    }
    for pair in path.windows(2) {
        brep.half_edges[pair[0]].next = pair[1];
        brep.half_edges[pair[1]].prev = pair[0];
    }
    Ok(())
}

fn verify_multi_run_bridge_face_postconditions(
    brep: &BRep,
    before: MultiRunBridgeFaceCounts,
    expected_edges: i64,
    plans: &[MultiRunBridgeFacePlan],
    bridge_loops: &[LoopHandle],
) -> Result<(), CoreError> {
    let expected = MultiRunBridgeFaceCounts {
        edges: expected_edges,
        faces: plans.len() as i64,
        loops: plans.len() as i64,
    };
    let actual = MultiRunBridgeFaceCounts::capture(brep).delta_from(before);
    if actual != expected {
        return Err(CoreError::PreconditionFailed(format!(
            "multi-run bridge face postcondition failed: expected_delta=dE={} dF={} dL={} actual_delta=dE={} dF={} dL={}",
            expected.edges,
            expected.faces,
            expected.loops,
            actual.edges,
            actual.faces,
            actual.loops
        )));
    }
    brep.validate_half_edge_pairing().map_err(CoreError::from)?;
    brep.validate_half_edge_cycles().map_err(CoreError::from)?;
    brep.validate_manifold_edges().map_err(CoreError::from)?;
    for &loop_ in bridge_loops {
        verify_multi_run_bridge_loop(brep, loop_)?;
    }
    Ok(())
}

fn verify_multi_run_bridge_loop(brep: &BRep, loop_: LoopHandle) -> Result<(), CoreError> {
    let start = brep
        .loops
        .get(loop_)
        .map(|loop_| loop_.half_edge)
        .ok_or_else(|| CoreError::InvalidHandle("multi-run bridge loop is stale".into()))?;
    let mut cursor = start;
    let mut visited = BTreeSet::new();
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor).ok_or_else(|| {
            CoreError::PreconditionFailed("multi-run bridge loop references stale half-edge".into())
        })?;
        if half_edge.loop_ != loop_ {
            return Err(CoreError::PreconditionFailed(
                "multi-run bridge cycle left its loop".into(),
            ));
        }
        if !visited.insert(cursor) {
            return Err(CoreError::PreconditionFailed(
                "multi-run bridge cycle repeats before closing".into(),
            ));
        }
        cursor = half_edge.next;
        if cursor == start {
            return Ok(());
        }
    }
    Err(CoreError::PreconditionFailed(
        "multi-run bridge cycle did not close".into(),
    ))
}

impl MultiRunBridgeFaceCounts {
    fn capture(brep: &BRep) -> Self {
        Self {
            edges: brep.edges.len() as i64,
            faces: brep.faces.len() as i64,
            loops: brep.loops.len() as i64,
        }
    }

    fn delta_from(self, before: Self) -> Self {
        Self {
            edges: self.edges - before.edges,
            faces: self.faces - before.faces,
            loops: self.loops - before.loops,
        }
    }
}

fn half_open_path(interval: &[HalfEdgeHandle]) -> Result<&[HalfEdgeHandle], CoreError> {
    if interval.len() < 2 {
        return Err(CoreError::PreconditionFailed(
            "multi-run bridge interval needs distinct endpoint origins".into(),
        ));
    }
    Ok(&interval[..interval.len() - 1])
}

fn required_path_handle(handle: Option<&HalfEdgeHandle>) -> Result<HalfEdgeHandle, CoreError> {
    handle.copied().ok_or_else(|| {
        CoreError::PreconditionFailed("multi-run bridge path handle is missing".into())
    })
}
