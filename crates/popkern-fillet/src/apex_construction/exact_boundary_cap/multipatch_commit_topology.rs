//! Live topology materialization for exact-boundary multi-patch commits.
#![allow(clippy::result_large_err)]

use popkern_core::geometry::{CurveEval, CurveOnEdge, NurbsCurve, Point3};
use popkern_core::topology::{
    BRep, Edge, EdgeHandle, Face, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, TopologyRole, Vertex,
    VertexHandle,
};

use super::super::ExactBoundaryMultiPatchCandidate;
use super::super::inventory::ExactBoundaryOuterStripPatchInventory;
use super::super::report::ExactBoundaryMultiPatchDiagnostic;
use crate::apex_construction::exact_boundary_cap::sector::ExactBoundaryCapSectorInput;

#[derive(Debug, Clone, Copy)]
pub(super) struct MaterializedPatchFaces {
    pub(super) strip_face: FaceHandle,
    pub(super) inner_face: FaceHandle,
    pub(super) interface_edge: EdgeHandle,
    pub(super) interface_strip_half_edge: HalfEdgeHandle,
    pub(super) interface_inner_half_edge: HalfEdgeHandle,
    pub(super) right_connector_edge: EdgeHandle,
    pub(super) left_connector_edge: EdgeHandle,
}

pub(super) fn materialize_patch_faces(
    brep: &mut BRep,
    sector: &ExactBoundaryCapSectorInput,
    candidate: &ExactBoundaryMultiPatchCandidate,
) -> Result<MaterializedPatchFaces, ExactBoundaryMultiPatchDiagnostic> {
    let strip = &candidate.strip_inventory;
    let source_loop = brep
        .faces
        .get(sector.cap_face)
        .map(|face| face.outer_loop)
        .ok_or_else(|| materialization_failed(sector, strip, "cap_face_missing", None))?;
    let loop_half_edges = collect_loop_half_edges(brep, source_loop)
        .ok_or_else(|| materialization_failed(sector, strip, "source_loop_invalid", None))?;
    if !loop_half_edges.contains(&strip.cap_half_edge) {
        return Err(materialization_failed(
            sector,
            strip,
            "outer_half_edge_not_on_sector_loop",
            None,
        ));
    }

    let rest_path = rest_path_after(loop_half_edges, strip.cap_half_edge);
    if rest_path.is_empty() {
        return Err(materialization_failed(
            sector,
            strip,
            "source_loop_has_no_inner_region",
            None,
        ));
    }

    let outer_start = brep.half_edges[strip.cap_half_edge].origin;
    let outer_end = half_edge_destination_vertex(brep, strip.cap_half_edge).ok_or_else(|| {
        materialization_failed(sector, strip, "outer_half_edge_destination_missing", None)
    })?;
    let (interface_start, interface_end) = interface_endpoints(&candidate.interface.curve);
    let interface_start_vertex = insert_generated_vertex(brep, interface_start);
    let interface_end_vertex = insert_generated_vertex(brep, interface_end);

    let (right_connector_edge, right_connector, right_connector_twin) =
        insert_generated_edge_pair(brep, outer_end, interface_end_vertex);
    let (interface_edge, interface_inner, interface_strip) =
        insert_generated_edge_pair(brep, interface_start_vertex, interface_end_vertex);
    let (left_connector_edge, left_connector, left_connector_twin) =
        insert_generated_edge_pair(brep, interface_start_vertex, outer_start);

    brep.edges[right_connector_edge].curve = CurveOnEdge::line_segment(
        brep.vertices[outer_end].position,
        interface_end,
        brep.tolerance.linear(),
    );
    brep.edges[interface_edge].curve = Some(CurveOnEdge {
        geometry: popkern_core::geometry::CurveGeometry::Nurbs(candidate.interface.curve.clone()),
        t_start: candidate.interface.curve.parameter_range().0,
        t_end: candidate.interface.curve.parameter_range().1,
    });
    brep.edges[left_connector_edge].curve = CurveOnEdge::line_segment(
        interface_start,
        brep.vertices[outer_start].position,
        brep.tolerance.linear(),
    );

    let inner_face = insert_inner_face_like(brep, sector.cap_face)
        .ok_or_else(|| materialization_failed(sector, strip, "inner_face_insert_failed", None))?;
    let inner_loop = brep.faces[inner_face].outer_loop;
    wire_loop(
        brep,
        source_loop,
        sector.cap_face,
        &[
            strip.cap_half_edge,
            right_connector,
            interface_strip,
            left_connector,
        ],
        sector,
        strip,
    )?;
    let mut inner_cycle = vec![interface_inner, right_connector_twin];
    inner_cycle.extend(rest_path);
    inner_cycle.push(left_connector_twin);
    wire_loop(brep, inner_loop, inner_face, &inner_cycle, sector, strip)?;

    Ok(MaterializedPatchFaces {
        strip_face: sector.cap_face,
        inner_face,
        interface_edge,
        interface_strip_half_edge: interface_strip,
        interface_inner_half_edge: interface_inner,
        right_connector_edge,
        left_connector_edge,
    })
}

fn materialization_failed(
    sector: &ExactBoundaryCapSectorInput,
    strip: &ExactBoundaryOuterStripPatchInventory,
    detail_code: &'static str,
    residual: Option<f64>,
) -> ExactBoundaryMultiPatchDiagnostic {
    ExactBoundaryMultiPatchDiagnostic::commit_materialization_failed(
        sector,
        strip,
        detail_code,
        residual,
    )
}

fn collect_loop_half_edges(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let start = brep.loops.get(loop_handle)?.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            return Some(out);
        }
    }
    None
}

fn rest_path_after(
    loop_half_edges: Vec<HalfEdgeHandle>,
    selected: HalfEdgeHandle,
) -> Vec<HalfEdgeHandle> {
    let Some(index) = loop_half_edges
        .iter()
        .position(|half_edge| *half_edge == selected)
    else {
        return Vec::new();
    };
    let mut rest = Vec::with_capacity(loop_half_edges.len().saturating_sub(1));
    for offset in 1..loop_half_edges.len() {
        rest.push(loop_half_edges[(index + offset) % loop_half_edges.len()]);
    }
    rest
}

fn half_edge_destination_vertex(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<VertexHandle> {
    let data = brep.half_edges.get(half_edge)?;
    brep.half_edges
        .get(data.twin)
        .map(|twin| twin.origin)
        .or_else(|| brep.half_edges.get(data.next).map(|next| next.origin))
}

fn interface_endpoints(curve: &NurbsCurve) -> (Point3, Point3) {
    let (t_start, t_end) = curve.parameter_range();
    (curve.point_at(t_start), curve.point_at(t_end))
}

fn insert_generated_vertex(brep: &mut BRep, position: Point3) -> VertexHandle {
    brep.vertices.insert(Vertex {
        position,
        outgoing: HalfEdgeHandle::default(),
        name: None,
        tolerance: None,
    })
}

fn insert_generated_edge_pair(
    brep: &mut BRep,
    start: VertexHandle,
    end: VertexHandle,
) -> (EdgeHandle, HalfEdgeHandle, HalfEdgeHandle) {
    let edge = brep.edges.insert(Edge {
        half_edge: HalfEdgeHandle::default(),
        curve: None,
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: TopologyRole::Production,
    });
    let half_edge = brep.half_edges.insert(HalfEdge {
        origin: start,
        twin: HalfEdgeHandle::default(),
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge,
        loop_: Default::default(),
        pcurve: None,
    });
    let twin = brep.half_edges.insert(HalfEdge {
        origin: end,
        twin: half_edge,
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge,
        loop_: Default::default(),
        pcurve: None,
    });
    brep.half_edges[half_edge].twin = twin;
    brep.edges[edge].half_edge = half_edge;
    if brep.vertices[start].outgoing == HalfEdgeHandle::default() {
        brep.vertices[start].outgoing = half_edge;
    }
    if brep.vertices[end].outgoing == HalfEdgeHandle::default() {
        brep.vertices[end].outgoing = twin;
    }
    (edge, half_edge, twin)
}

fn insert_inner_face_like(brep: &mut BRep, template_face: FaceHandle) -> Option<FaceHandle> {
    let template = brep.faces.get(template_face)?.clone();
    let face = brep.faces.insert(Face {
        outer_loop: Default::default(),
        inner_loops: Vec::new(),
        shell: template.shell,
        surface: None,
        name: None,
        tolerance: template.tolerance,
        role: TopologyRole::Production,
    });
    let loop_handle = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face,
    });
    brep.faces[face].outer_loop = loop_handle;
    Some(face)
}

fn wire_loop(
    brep: &mut BRep,
    loop_handle: popkern_core::topology::LoopHandle,
    face: FaceHandle,
    half_edges: &[HalfEdgeHandle],
    sector: &ExactBoundaryCapSectorInput,
    strip: &ExactBoundaryOuterStripPatchInventory,
) -> Result<(), ExactBoundaryMultiPatchDiagnostic> {
    let n = half_edges.len();
    if n == 0 {
        return Err(materialization_failed(sector, strip, "empty_loop", None));
    }
    for i in 0..n {
        let current = half_edges[i];
        let next = half_edges[(i + 1) % n];
        let prev = half_edges[(i + n - 1) % n];
        let half_edge = brep
            .half_edges
            .get_mut(current)
            .ok_or_else(|| materialization_failed(sector, strip, "stale_half_edge", None))?;
        half_edge.next = next;
        half_edge.prev = prev;
        half_edge.loop_ = loop_handle;
    }
    brep.loops[loop_handle].half_edge = half_edges[0];
    brep.loops[loop_handle].face = face;
    Ok(())
}
