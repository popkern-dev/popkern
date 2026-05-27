//! Support-face seam materialization for exact-boundary commit.

use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{
    BRep, FaceHandle, HalfEdge, HalfEdgeHandle, LoopHandle, TopologyRole,
};

use crate::apex_construction::repeated_apex_partition_topology::{
    insert_boundary_half_edge, insert_face_with_outer_loop,
};

/// Support-owned seam authored for an exact-boundary cap/support edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct ExactBoundarySupportSubdivision {
    /// New support-side half-edge twinned to the cap half-edge.
    pub(super) support_half_edge: HalfEdgeHandle,
    /// New support patch face carrying the prepared support surface.
    pub(super) support_face: FaceHandle,
}

/// Materialize a support-owned seam when the current cap twin belongs to a
/// bridge/residual loop whose whole loop cannot carry the prepared support
/// surface.
pub(super) fn materialize_support_subdivision(
    brep: &mut BRep,
    cap_half_edge: HalfEdgeHandle,
    support_candidate_half_edge: HalfEdgeHandle,
    template_support_face: FaceHandle,
    support_surface: SurfaceGeometry,
) -> Option<ExactBoundarySupportSubdivision> {
    let cap = brep.half_edges.get(cap_half_edge)?;
    let cap_edge = cap.edge;
    let cap_start = cap.origin;
    let previous_cap_twin = cap.twin;
    let cap_end = brep.half_edges.get(previous_cap_twin)?.origin;
    let candidate_face = brep
        .half_edges
        .get(support_candidate_half_edge)
        .and_then(|half_edge| brep.loops.get(half_edge.loop_))
        .map(|loop_| loop_.face)?;
    if candidate_face != template_support_face {
        return None;
    }
    let start_point = brep.vertices.get(cap_start)?.position;
    let end_point = brep.vertices.get(cap_end)?.position;
    if !endpoints_lie_on_support(
        &support_surface,
        start_point,
        end_point,
        brep.tolerance.linear(),
    ) {
        return None;
    }

    let rollback = brep.clone();
    let result = (|| {
        preserve_previous_cap_neighbor_side(brep, previous_cap_twin, cap_start, cap_end)?;
        let support_half_edge = brep.half_edges.insert(HalfEdge {
            origin: cap_end,
            twin: cap_half_edge,
            next: HalfEdgeHandle::default(),
            prev: HalfEdgeHandle::default(),
            edge: cap_edge,
            loop_: LoopHandle::default(),
            pcurve: None,
        });
        brep.half_edges[cap_half_edge].twin = support_half_edge;
        brep.edges[cap_edge].half_edge = cap_half_edge;

        let closing = author_support_patch_closing_side(brep, cap_start, cap_end)?;
        let support_face =
            insert_face_with_outer_loop(brep, template_support_face, &[support_half_edge, closing])
                .ok()?;
        brep.faces[support_face].surface = Some(SurfaceOnFace {
            geometry: support_surface,
            orientation: true,
        });
        validate_support_subdivision(brep, cap_half_edge, support_half_edge, support_face)?;
        Some(ExactBoundarySupportSubdivision {
            support_half_edge,
            support_face,
        })
    })();
    match result {
        Some(subdivision) => Some(subdivision),
        None => {
            *brep = rollback;
            None
        }
    }
}

fn preserve_previous_cap_neighbor_side(
    brep: &mut BRep,
    previous_cap_twin: HalfEdgeHandle,
    cap_start: popkern_core::topology::VertexHandle,
    cap_end: popkern_core::topology::VertexHandle,
) -> Option<()> {
    let old_edge = brep.half_edges.get(previous_cap_twin)?.edge;
    let old_curve = brep.edges.get(old_edge)?.curve.clone();
    let old_role = brep.edges.get(old_edge)?.role;
    let old_face = brep
        .half_edges
        .get(previous_cap_twin)
        .and_then(|half_edge| brep.loops.get(half_edge.loop_))
        .map(|loop_| loop_.face)?;
    let preserved_edge = brep.edges.insert(popkern_core::topology::Edge {
        half_edge: previous_cap_twin,
        curve: old_curve,
        name: None,
        tolerance: None,
        canonical: Default::default(),
        role: old_role,
    });
    brep.half_edges[previous_cap_twin].edge = preserved_edge;

    let balancing_half_edge = brep.half_edges.insert(HalfEdge {
        origin: cap_start,
        twin: previous_cap_twin,
        next: HalfEdgeHandle::default(),
        prev: HalfEdgeHandle::default(),
        edge: preserved_edge,
        loop_: LoopHandle::default(),
        pcurve: None,
    });
    brep.half_edges[previous_cap_twin].twin = balancing_half_edge;
    let closing = author_support_patch_closing_side(brep, cap_end, cap_start)?;
    let residual_face =
        insert_face_with_outer_loop(brep, old_face, &[balancing_half_edge, closing]).ok()?;
    brep.faces[residual_face].role = TopologyRole::ConstructionResidual;
    Some(())
}

fn author_support_patch_closing_side(
    brep: &mut BRep,
    start: popkern_core::topology::VertexHandle,
    end: popkern_core::topology::VertexHandle,
) -> Option<HalfEdgeHandle> {
    let start_point = brep.vertices.get(start)?.position;
    let end_point = brep.vertices.get(end)?.position;
    let curve = CurveOnEdge::line_segment(start_point, end_point, brep.tolerance.linear())?;
    let (edge, half_edge) = insert_boundary_half_edge(brep, start);
    brep.edges[edge].curve = Some(curve);
    Some(half_edge)
}

fn endpoints_lie_on_support(
    support_surface: &SurfaceGeometry,
    start: popkern_core::geometry::Point3,
    end: popkern_core::geometry::Point3,
    tolerance: f64,
) -> bool {
    match support_surface {
        SurfaceGeometry::Plane { origin, normal, .. } => {
            let tolerance = tolerance.max(1.0e-9);
            (start - *origin).dot(normal).abs() <= tolerance
                && (end - *origin).dot(normal).abs() <= tolerance
        }
        _ => true,
    }
}

fn validate_support_subdivision(
    brep: &BRep,
    cap_half_edge: HalfEdgeHandle,
    support_half_edge: HalfEdgeHandle,
    support_face: FaceHandle,
) -> Option<()> {
    let cap = brep.half_edges.get(cap_half_edge)?;
    let support = brep.half_edges.get(support_half_edge)?;
    if cap.twin != support_half_edge || support.twin != cap_half_edge || cap.edge != support.edge {
        return None;
    }
    let owner = brep
        .loops
        .get(support.loop_)
        .map(|loop_data| loop_data.face)?;
    (owner == support_face).then_some(())
}
