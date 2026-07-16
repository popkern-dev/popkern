//! Sector-face topology materialization for exact-boundary apex caps.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::topology::{
    BRep, FaceHandle, HalfEdgeHandle, TopologyRole, Vertex, VertexHandle,
};

#[path = "topology_outer_seams.rs"]
mod topology_outer_seams;

use super::corner::{PreparedExactBoundaryInnerSpokes, PreparedExactBoundarySectorBoundary};
use super::sector::ExactBoundaryCapSectorInput;
use topology_outer_seams::refresh_sector_outer_seams;

use crate::apex_construction::repeated_apex_partition_topology::{
    insert_face_with_outer_loop, insert_unwired_edge_pair, wire_existing_loop,
};

/// Sector topology ready for exact-boundary BRep commit.
#[derive(Debug, Clone)]
pub(crate) struct ExactBoundarySectorTopology {
    /// Sector inputs remapped to the materialized sector face handles.
    pub sectors: Vec<ExactBoundaryCapSectorInput>,
    /// Shared pole vertex used by all generated inner spokes.
    pub pole_vertex: VertexHandle,
    /// Generated spoke edges from boundary anchors to the pole.
    pub spoke_edges: Vec<popkern_core::topology::EdgeHandle>,
}

/// Stable decline reason for sector topology materialization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ExactBoundaryTopologyDeclineReason {
    /// The aggregate cap face is stale.
    CapFaceMissing,
    /// Sector and generated-spoke inventories disagree.
    SectorInventoryMismatch,
    /// A sector's two cap-side outer half-edges are not consecutive.
    OuterHalfEdgesNotConsecutive,
    /// A generated spoke anchor could not be matched to a live boundary vertex.
    SpokeAnchorVertexMissing,
    /// A generated spoke edge for a sector corner is missing.
    SpokeEdgeMissing,
    /// Low-level topology wiring failed.
    TopologyWiringFailed,
    /// Final BRep validation failed, so the transaction was rolled back.
    ValidationFailed,
}

impl ExactBoundaryTopologyDeclineReason {
    /// Stable diagnostic code for tests and future UI/API consumers.
    pub(crate) fn code(self) -> &'static str {
        match self {
            Self::CapFaceMissing => "exact_boundary_cap.topology.cap_face_missing",
            Self::SectorInventoryMismatch => "exact_boundary_cap.topology.inventory_mismatch",
            Self::OuterHalfEdgesNotConsecutive => {
                "exact_boundary_cap.topology.outer_half_edges_not_consecutive"
            }
            Self::SpokeAnchorVertexMissing => {
                "exact_boundary_cap.topology.spoke_anchor_vertex_missing"
            }
            Self::SpokeEdgeMissing => "exact_boundary_cap.topology.spoke_edge_missing",
            Self::TopologyWiringFailed => "exact_boundary_cap.topology.wiring_failed",
            Self::ValidationFailed => "exact_boundary_cap.topology.validation_failed",
        }
    }
}

/// Structured topology-materialization decline.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ExactBoundaryTopologyDiagnostic {
    /// Stable decline reason.
    pub reason: ExactBoundaryTopologyDeclineReason,
    /// Stable diagnostic code.
    pub code: &'static str,
    /// Stable sector name when the failure is sector-scoped.
    pub sector_name: Option<String>,
    /// Measured residual when available.
    pub residual: Option<f64>,
}

impl ExactBoundaryTopologyDiagnostic {
    fn new(
        reason: ExactBoundaryTopologyDeclineReason,
        sector_name: Option<String>,
        residual: Option<f64>,
    ) -> Self {
        Self {
            reason,
            code: reason.code(),
            sector_name,
            residual,
        }
    }

    pub(crate) fn into_apex_diagnostic(self) -> crate::types::ApexConstructionDiagnostic {
        crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
            crate::types::ExactBoundaryCapDiagnostic {
                stage: "exact_boundary_cap.topology".to_string(),
                code: self.code.to_string(),
                detail_code: None,
                sector_name: self
                    .sector_name
                    .unwrap_or_else(|| "exact_boundary_cap".to_string()),
                edge: None,
                support_half_edge: None,
                support_face: None,
                residual: self.residual,
                outer_side: None,
                sample_parameter: None,
                tolerance: None,
                shared_fact_ids: Vec::new(),
            },
        )
    }
}

/// Materialize one sector face per exact-boundary sector under rollback.
pub(crate) fn materialize_exact_boundary_sector_topology(
    brep: &mut BRep,
    cap_face: FaceHandle,
    sectors: &[ExactBoundaryCapSectorInput],
    spokes: &PreparedExactBoundaryInnerSpokes,
) -> Result<ExactBoundarySectorTopology, ExactBoundaryTopologyDiagnostic> {
    let snapshot = brep.clone();
    let result = materialize_validated_sector_topology(brep, cap_face, sectors, spokes);
    if result.is_err() {
        *brep = snapshot;
    }
    result
}

fn materialize_validated_sector_topology(
    brep: &mut BRep,
    cap_face: FaceHandle,
    sectors: &[ExactBoundaryCapSectorInput],
    spokes: &PreparedExactBoundaryInnerSpokes,
) -> Result<ExactBoundarySectorTopology, ExactBoundaryTopologyDiagnostic> {
    if !brep.faces.contains_key(cap_face) {
        return Err(ExactBoundaryTopologyDiagnostic::new(
            ExactBoundaryTopologyDeclineReason::CapFaceMissing,
            None,
            None,
        ));
    }
    if sectors.len() != spokes.sectors.len() || sectors.is_empty() {
        return Err(ExactBoundaryTopologyDiagnostic::new(
            ExactBoundaryTopologyDeclineReason::SectorInventoryMismatch,
            None,
            None,
        ));
    }

    let pole_vertex = brep.vertices.insert(Vertex {
        position: spokes.pole,
        outgoing: HalfEdgeHandle::default(),
        name: None,
        tolerance: None,
    });
    let mut spoke_edges = Vec::with_capacity(sectors.len());
    let mut spoke_half_edges = Vec::with_capacity(sectors.len());
    for sector in &spokes.sectors {
        let vertex = vertex_for_spoke_anchor(brep, sector)?;
        let (edge, anchor_to_pole, pole_to_anchor) =
            insert_unwired_edge_pair(brep, vertex, pole_vertex);
        brep.edges[edge].role = TopologyRole::Production;
        // Semantic line storage (PK-092): a spoke whose IMAGE is
        // straight is stored as an analytic Line at construction
        // time. Nothing has consumed this edge yet — every pcurve,
        // split, and export downstream derives from the edge curve —
        // so canonicalizing here is contract-neutral, whereas the
        // post-hoc normalization pass must (and does) roll back once
        // pcurves have been fitted to the NURBS parameterization.
        let spoke_geometry = CurveGeometry::Nurbs(sector.face_inner_spoke.curve.clone());
        let spoke_geometry = spoke_geometry
            .canonical_line_in_range((0.0, 1.0), brep.tolerance.linear().max(1.0e-8))
            .unwrap_or(spoke_geometry);
        brep.edges[edge].curve = Some(CurveOnEdge {
            geometry: spoke_geometry,
            t_start: 0.0,
            t_end: 1.0,
        });
        spoke_edges.push(edge);
        spoke_half_edges.push((anchor_to_pole, pole_to_anchor));
    }

    let mut remapped = sectors.to_vec();
    for index in 0..sectors.len() {
        let sector = &sectors[index];
        let previous_index = (index + sectors.len() - 1) % sectors.len();
        let cycle = sector_cycle(
            brep,
            sector,
            spoke_half_edges[previous_index],
            spoke_half_edges[index],
        )?;
        let sector_face = if index == 0 {
            let loop_handle = brep.faces[cap_face].outer_loop;
            wire_existing_loop(brep, loop_handle, cap_face, &cycle).map_err(|_| {
                ExactBoundaryTopologyDiagnostic::new(
                    ExactBoundaryTopologyDeclineReason::TopologyWiringFailed,
                    Some(sector.sector_name.clone()),
                    None,
                )
            })?;
            cap_face
        } else {
            insert_face_with_outer_loop(brep, cap_face, &cycle).map_err(|_| {
                ExactBoundaryTopologyDiagnostic::new(
                    ExactBoundaryTopologyDeclineReason::TopologyWiringFailed,
                    Some(sector.sector_name.clone()),
                    None,
                )
            })?
        };
        remapped[index].cap_face = sector_face;
    }
    for sector in &mut remapped {
        refresh_sector_outer_seams(brep, sector)?;
    }

    if let Err(err) = brep.validate() {
        return Err(ExactBoundaryTopologyDiagnostic::new(
            ExactBoundaryTopologyDeclineReason::ValidationFailed,
            None,
            Some(err.to_string().len() as f64),
        ));
    }
    Ok(ExactBoundarySectorTopology {
        sectors: remapped,
        pole_vertex,
        spoke_edges,
    })
}

fn vertex_for_spoke_anchor(
    brep: &BRep,
    sector: &PreparedExactBoundarySectorBoundary,
) -> Result<VertexHandle, ExactBoundaryTopologyDiagnostic> {
    let tolerance = brep.tolerance.linear().max(1.0e-7);
    brep.vertices
        .iter()
        .filter_map(|(vertex, data)| {
            let distance = data.position.distance_to(&sector.face_inner_spoke.anchor);
            (distance <= tolerance).then_some((vertex, distance))
        })
        .min_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(vertex, _)| vertex)
        .ok_or_else(|| {
            ExactBoundaryTopologyDiagnostic::new(
                ExactBoundaryTopologyDeclineReason::SpokeAnchorVertexMissing,
                Some(sector.sector_name.clone()),
                None,
            )
        })
}

fn sector_cycle(
    brep: &BRep,
    sector: &ExactBoundaryCapSectorInput,
    previous_spoke: (HalfEdgeHandle, HalfEdgeHandle),
    current_spoke: (HalfEdgeHandle, HalfEdgeHandle),
) -> Result<Vec<HalfEdgeHandle>, ExactBoundaryTopologyDiagnostic> {
    let fillet = sector.cap_half_edges[0];
    let face = sector.cap_half_edges[1];
    let fillet_data = brep.half_edges.get(fillet).ok_or_else(|| {
        ExactBoundaryTopologyDiagnostic::new(
            ExactBoundaryTopologyDeclineReason::OuterHalfEdgesNotConsecutive,
            Some(sector.sector_name.clone()),
            None,
        )
    })?;
    let face_data = brep.half_edges.get(face).ok_or_else(|| {
        ExactBoundaryTopologyDiagnostic::new(
            ExactBoundaryTopologyDeclineReason::OuterHalfEdgesNotConsecutive,
            Some(sector.sector_name.clone()),
            None,
        )
    })?;
    let outer_half_edges = if fillet_data.next == face {
        [fillet, face]
    } else if face_data.next == fillet {
        [face, fillet]
    } else {
        return Err(ExactBoundaryTopologyDiagnostic::new(
            ExactBoundaryTopologyDeclineReason::OuterHalfEdgesNotConsecutive,
            Some(sector.sector_name.clone()),
            None,
        ));
    };
    sector_cycle_from_endpoints(brep, outer_half_edges, previous_spoke, current_spoke).ok_or_else(
        || {
            ExactBoundaryTopologyDiagnostic::new(
                ExactBoundaryTopologyDeclineReason::SpokeEdgeMissing,
                Some(sector.sector_name.clone()),
                None,
            )
        },
    )
}

fn sector_cycle_from_endpoints(
    brep: &BRep,
    outer_half_edges: [HalfEdgeHandle; 2],
    previous_spoke: (HalfEdgeHandle, HalfEdgeHandle),
    current_spoke: (HalfEdgeHandle, HalfEdgeHandle),
) -> Option<Vec<HalfEdgeHandle>> {
    let spokes = [
        previous_spoke.0,
        previous_spoke.1,
        current_spoke.0,
        current_spoke.1,
    ];
    for first_spoke in spokes {
        for second_spoke in spokes {
            if first_spoke == second_spoke {
                continue;
            }
            let candidate = vec![
                outer_half_edges[0],
                outer_half_edges[1],
                first_spoke,
                second_spoke,
            ];
            if cycle_is_contiguous(brep, &candidate) {
                return Some(candidate);
            }
        }
    }
    None
}

fn cycle_is_contiguous(brep: &BRep, cycle: &[HalfEdgeHandle]) -> bool {
    !cycle.is_empty()
        && cycle.iter().enumerate().all(|(index, half_edge)| {
            let next = cycle[(index + 1) % cycle.len()];
            half_edges_connect(brep, *half_edge, next)
        })
}

fn half_edges_connect(brep: &BRep, left: HalfEdgeHandle, right: HalfEdgeHandle) -> bool {
    let Some(left_end) = half_edge_destination_vertex(brep, left) else {
        return false;
    };
    let Some(right_start) = brep.half_edges.get(right).map(|half_edge| half_edge.origin) else {
        return false;
    };
    if left_end == right_start {
        return true;
    }
    brep.vertices
        .get(left_end)
        .zip(brep.vertices.get(right_start))
        .is_some_and(|(left_vertex, right_vertex)| {
            left_vertex.position.distance_to(&right_vertex.position)
                <= brep.tolerance.linear().max(1.0e-9)
        })
}

fn half_edge_destination_vertex(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<VertexHandle> {
    let data = brep.half_edges.get(half_edge)?;
    brep.half_edges
        .get(data.twin)
        .map(|twin| twin.origin)
        .or_else(|| brep.half_edges.get(data.next).map(|next| next.origin))
}
