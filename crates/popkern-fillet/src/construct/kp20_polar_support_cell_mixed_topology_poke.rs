//! Topology commit for complete mixed KP20 support-cell patchwork.

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{GregoryPatch, Kp20PieceId};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::construct::kp20_polar_expanded_support_binding::{
    Kp20ExpandedSupportSurfaceAssignment, Kp20ExpandedSupportSurfaceSlot,
    attach_kp20_expanded_support_surfaces, author_kp20_expanded_support_pcurves,
    plan_kp20_expanded_support_surface_bindings_with_sources,
};
use crate::construct::kp20_polar_support_cell_inventory::Kp20SupportCellInventory;

#[path = "kp20_polar_support_cell_topology_pattern.rs"]
#[allow(dead_code)]
mod support_cell_topology_pattern;

use support_cell_topology_pattern::Kp20SupportCellTopologyPattern;

#[derive(Debug, Clone)]
pub(crate) struct Kp20MixedSupportCellCommitResult {
    pub(crate) sub_faces: Vec<FaceHandle>,
    pub(crate) support_surface_indices: Vec<usize>,
}

pub(crate) fn commit_kp20_mixed_support_cell_patchwork(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    sector_count: usize,
    support_surfaces: &[GregoryPatch],
    inventory: &Kp20SupportCellInventory,
    target_piece: Kp20PieceId,
) -> Result<Kp20MixedSupportCellCommitResult, CoreError> {
    let before = brep.clone();
    match commit_kp20_mixed_support_cell_patchwork_inner(
        brep,
        face,
        central_point,
        sector_count,
        support_surfaces,
        inventory,
        target_piece,
    ) {
        Ok(result) => Ok(result),
        Err(error) => {
            *brep = before;
            Err(error)
        }
    }
}

fn commit_kp20_mixed_support_cell_patchwork_inner(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    sector_count: usize,
    support_surfaces: &[GregoryPatch],
    inventory: &Kp20SupportCellInventory,
    target_piece: Kp20PieceId,
) -> Result<Kp20MixedSupportCellCommitResult, CoreError> {
    if inventory.sector_count != sector_count {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 mixed support-cell commit inventory sector count {} does not match ring sector count {sector_count}",
            inventory.sector_count,
        )));
    }
    let expected_surface_count =
        4 * sector_count + inventory.target_sectors.len() * (inventory.cells_per_target_sector - 1);
    if support_surfaces.len() != expected_surface_count {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 mixed support-cell commit expected {expected_surface_count} surfaces, got {}",
            support_surfaces.len(),
        )));
    }
    let Some(pattern) =
        Kp20SupportCellTopologyPattern::for_cells_per_sector(inventory.cells_per_target_sector)
    else {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 mixed support-cell commit does not support {} cells per target sector",
            inventory.cells_per_target_sector,
        )));
    };
    if !matches!(
        target_piece,
        Kp20PieceId::NearMidpointRight | Kp20PieceId::NearMidpointLeft
    ) {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 mixed support-cell commit requires a support-owner target piece, got {target_piece:?}",
        )));
    }
    let sub_faces = poke_mixed_support_cell_faces(
        brep,
        face,
        central_point,
        sector_count,
        inventory,
        pattern,
        target_piece,
    )?;
    let slots = mixed_support_cell_surface_slots(brep, &sub_faces)?;
    let source_map = (0..sub_faces.len()).collect::<Vec<_>>();
    let binding_plan = plan_kp20_expanded_support_surface_bindings_with_sources(
        &slots,
        support_surfaces.len(),
        &source_map,
    )?;
    attach_kp20_expanded_support_surfaces(brep, &binding_plan, support_surfaces)?;
    author_kp20_expanded_support_pcurves(brep, &binding_plan)?;
    let support_surface_indices = binding_plan
        .bindings
        .iter()
        .map(|binding| binding.support_surface_index)
        .collect();
    Ok(Kp20MixedSupportCellCommitResult {
        sub_faces,
        support_surface_indices,
    })
}

fn poke_mixed_support_cell_faces(
    brep: &mut BRep,
    face: FaceHandle,
    central_point: Point3,
    sector_count: usize,
    inventory: &Kp20SupportCellInventory,
    pattern: Kp20SupportCellTopologyPattern,
    target_piece: Kp20PieceId,
) -> Result<Vec<FaceHandle>, CoreError> {
    let sector_poke = crate::construct::poke_hex::poke_n_gon_face_into_expected_n_quads(
        brep,
        face,
        central_point,
        sector_count,
    )?;
    split_radial_edges_for_kp20_semantic_poke(brep, sector_poke.central_vertex, sector_count)?;
    let semantic_faces_by_sector = sector_poke
        .sub_faces
        .into_iter()
        .map(|sector_face| {
            subdivide_sector_into_kp20_semantic_four(brep, sector_face, sector_poke.central_vertex)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut sub_faces = Vec::new();
    for (sector, semantic_faces) in semantic_faces_by_sector.into_iter().enumerate() {
        for piece in Kp20PieceId::all() {
            let semantic_face = semantic_faces[piece.index()];
            if piece == target_piece && inventory.target_sectors.contains(&sector) {
                sub_faces.extend(pattern.split_sector_support_face(brep, semantic_face)?);
            } else {
                sub_faces.push(semantic_face);
            }
        }
    }
    Ok(sub_faces)
}

fn split_radial_edges_for_kp20_semantic_poke(
    brep: &mut BRep,
    central_vertex: VertexHandle,
    sector_count: usize,
) -> Result<(), CoreError> {
    let radial_edges = collect_vertex_outgoing_edges(brep, central_vertex)?;
    if radial_edges.len() != sector_count {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 semantic mixed support-cell poke central vertex has {} outgoing edges, expected {sector_count}",
            radial_edges.len()
        )));
    }
    for edge in radial_edges {
        brep.split_edge(edge, 0.5)?;
    }
    Ok(())
}

fn subdivide_sector_into_kp20_semantic_four(
    brep: &mut BRep,
    face: FaceHandle,
    central_vertex: VertexHandle,
) -> Result<[FaceHandle; 4], CoreError> {
    let hes = walk_loop_from_vertex(brep, face, central_vertex)?;
    if hes.len() != 6 {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 semantic mixed support-cell poke expected 6-edge sector after radial splits, got {}",
            hes.len()
        )));
    }
    let left_mid_vertex = brep.half_edges[hes[1]].origin;
    let outer_mid_vertex = brep.half_edges[hes[3]].origin;
    let right_mid_vertex = brep.half_edges[hes[5]].origin;
    let left_mid = brep.vertices[left_mid_vertex].position;
    let outer_mid = brep.vertices[outer_mid_vertex].position;
    let right_mid = brep.vertices[right_mid_vertex].position;
    let central = brep.vertices[central_vertex].position;
    let midpoint_lane = average_points(&[left_mid, outer_mid, right_mid]);
    let origin_lane = central + (midpoint_lane - central) * 0.5;

    let mid_mev = brep.mev(hes[1], midpoint_lane)?;
    let midpoint_vertex = mid_mev.vertex;
    let outer_mid_half_edge = half_edge_from_origin_on_face(brep, face, outer_mid_vertex)?;
    let near_midpoint_left = brep.mef(mid_mev.half_edge_twin, outer_mid_half_edge)?;

    let midpoint_half_edge = half_edge_from_origin_on_face(brep, face, midpoint_vertex)?;
    let right_mid_half_edge = half_edge_from_origin_on_face(brep, face, right_mid_vertex)?;
    let near_midpoint_right = brep.mef(midpoint_half_edge, right_mid_half_edge)?;

    let central_half_edge = half_edge_from_origin_on_face(brep, face, central_vertex)?;
    let origin_mev = brep.mev(central_half_edge, origin_lane)?;
    let midpoint_half_edge = half_edge_from_origin_on_face(brep, face, midpoint_vertex)?;
    let near_origin_left = brep.mef(origin_mev.half_edge_twin, midpoint_half_edge)?;

    Ok([
        near_origin_left.face,
        face,
        near_midpoint_right.face,
        near_midpoint_left.face,
    ])
}

fn collect_vertex_outgoing_edges(
    brep: &BRep,
    vertex: VertexHandle,
) -> Result<Vec<EdgeHandle>, CoreError> {
    let start = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| {
            CoreError::InvalidHandle("kp20 semantic mixed support-cell poke vertex invalid".into())
        })?
        .outgoing;
    if !brep.half_edges.contains_key(start) {
        return Err(CoreError::PreconditionFailed(
            "kp20 semantic mixed support-cell poke central vertex has no outgoing half-edge".into(),
        ));
    }
    let mut edges = Vec::new();
    let mut cursor = start;
    for _ in 0..brep.half_edges.len() {
        let he = &brep.half_edges[cursor];
        edges.push(he.edge);
        cursor = brep.half_edges[he.twin].next;
        if cursor == start {
            return Ok(edges);
        }
    }
    Err(CoreError::PreconditionFailed(
        "kp20 semantic mixed support-cell poke central vertex fan did not close".into(),
    ))
}

fn walk_loop_from_vertex(
    brep: &BRep,
    face: FaceHandle,
    start_vertex: VertexHandle,
) -> Result<Vec<HalfEdgeHandle>, CoreError> {
    let start = half_edge_from_origin_on_face(brep, face, start_vertex)?;
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return Ok(out);
        }
    }
    Err(CoreError::PreconditionFailed(
        "kp20 semantic mixed support-cell poke face loop did not close".into(),
    ))
}

fn half_edge_from_origin_on_face(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
) -> Result<HalfEdgeHandle, CoreError> {
    face_outer_half_edges(brep, face)?
        .into_iter()
        .find(|&half_edge| brep.half_edges[half_edge].origin == vertex)
        .ok_or_else(|| {
            CoreError::PreconditionFailed(
                "kp20 semantic mixed support-cell poke expected vertex on face loop".into(),
            )
        })
}

fn average_points(points: &[Point3]) -> Point3 {
    let scale = 1.0 / points.len() as f64;
    Point3::new(
        points.iter().map(|point| point.x).sum::<f64>() * scale,
        points.iter().map(|point| point.y).sum::<f64>() * scale,
        points.iter().map(|point| point.z).sum::<f64>() * scale,
    )
}

fn mixed_support_cell_surface_slots(
    brep: &BRep,
    sub_faces: &[FaceHandle],
) -> Result<Vec<Kp20ExpandedSupportSurfaceSlot>, CoreError> {
    let mut slots = Vec::with_capacity(sub_faces.len());
    for &face in sub_faces {
        let boundary_edge_count = face_outer_half_edges(brep, face)?.len();
        let assignment = match boundary_edge_count {
            4 => Kp20ExpandedSupportSurfaceAssignment::FullDomainRectangularPatch,
            5.. => Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace,
            len => {
                return Err(CoreError::PreconditionFailed(format!(
                    "kp20 mixed support-cell slot cannot classify {len}-edge face"
                )));
            }
        };
        slots.push(Kp20ExpandedSupportSurfaceSlot {
            face,
            boundary_edge_count,
            assignment,
        });
    }
    Ok(slots)
}

fn face_outer_half_edges(brep: &BRep, face: FaceHandle) -> Result<Vec<HalfEdgeHandle>, CoreError> {
    let outer_loop = brep
        .faces
        .get(face)
        .ok_or_else(|| CoreError::InvalidHandle("kp20 mixed support-cell face invalid".into()))?
        .outer_loop;
    let start = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| CoreError::InvalidHandle("kp20 mixed support-cell loop invalid".into()))?
        .half_edge;
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return Ok(out);
        }
    }
    Err(CoreError::PreconditionFailed(
        "kp20 mixed support-cell loop did not close".into(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::tolerance::Tolerance;
    use std::collections::BTreeSet;

    #[test]
    fn mixed_support_cell_commit_attaches_complete_inventory() {
        let (mut brep, face) = build_flat_2ngon(6);
        let support_surfaces = (0..33).map(dummy_gregory_surface).collect::<Vec<_>>();
        let inventory = Kp20SupportCellInventory {
            sector_count: 6,
            target_sectors: vec![0, 2, 5],
            cells_per_target_sector: 4,
        };

        let commit = commit_kp20_mixed_support_cell_patchwork(
            &mut brep,
            face,
            Point3::new(0.0, 0.0, 0.0),
            6,
            &support_surfaces,
            &inventory,
            Kp20PieceId::NearMidpointLeft,
        )
        .expect("mixed support-cell commit");

        assert_eq!(commit.sub_faces.len(), 33);
        assert_eq!(commit.support_surface_indices.len(), 33);
        assert!(commit.sub_faces.iter().all(|face| matches!(
            brep.faces[*face].surface.as_ref().map(|surface| &surface.geometry),
            Some(SurfaceGeometry::Gregory(_))
        )));
        brep.validate()
            .expect("mixed support-cell commit leaves valid BRep");
    }

    #[test]
    fn legacy_4n_topology_exposes_three_outer_slots() {
        let (mut brep, face) = build_flat_2ngon(6);
        let poke = crate::construct::poke_4n::poke_n_gon_face_into_expected_4n_quads(
            &mut brep,
            face,
            Point3::new(0.0, 0.0, 0.0),
            6,
        )
        .expect("legacy 4n poke");
        let cap_faces = poke.sub_faces.iter().copied().collect::<BTreeSet<_>>();
        let external_counts = (0..4)
            .map(|piece| external_boundary_edge_count(&brep, &cap_faces, poke.sub_faces[piece]))
            .collect::<Vec<_>>();

        assert_eq!(
            external_counts,
            vec![1, 2, 1, 0],
            "legacy [O_L,O_C,O_R,C_star] topology is not KP20-semantic: three slots touch four support-boundary edges and only C_star is central"
        );
    }

    #[test]
    fn semantic_kp20_topology_keeps_origin_slots_off_support_boundary() {
        let (mut brep, face) = build_flat_2ngon(6);
        let inventory = Kp20SupportCellInventory {
            sector_count: 6,
            target_sectors: Vec::new(),
            cells_per_target_sector: 4,
        };
        let sub_faces = poke_mixed_support_cell_faces(
            &mut brep,
            face,
            Point3::new(0.0, 0.0, 0.0),
            6,
            &inventory,
            Kp20SupportCellTopologyPattern::FourCellSupportGrid,
            Kp20PieceId::NearMidpointLeft,
        )
        .expect("semantic mixed support-cell poke");
        let cap_faces = sub_faces.iter().copied().collect::<BTreeSet<_>>();
        let external_counts = (0..4)
            .map(|piece| external_boundary_edge_count(&brep, &cap_faces, sub_faces[piece]))
            .collect::<Vec<_>>();

        assert_eq!(
            external_counts,
            vec![0, 0, 1, 1],
            "semantic KP20 topology keeps near-origin pieces internal and puts support-boundary ownership on near-midpoint pieces"
        );
        brep.validate()
            .expect("semantic KP20 mixed topology leaves valid BRep");
    }

    fn build_flat_2ngon(n: usize) -> (BRep, FaceHandle) {
        let mut brep = BRep::new(Tolerance::default());
        let count = 2 * n;
        let pts = (0..count)
            .map(|k| {
                let theta = std::f64::consts::TAU * (k as f64) / (count as f64);
                Point3::new(theta.cos(), theta.sin(), 0.0)
            })
            .collect::<Vec<_>>();
        let mvfs = brep.mvfs(pts[0]).expect("mvfs");
        let first = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, pts[1])
            .expect("first edge");
        let mut previous = first;
        for point in pts.iter().take(count).skip(2) {
            previous = brep.mev(previous.half_edge_twin, *point).expect("edge");
        }
        let close = brep
            .mef(previous.half_edge_twin, first.half_edge)
            .expect("close");
        (brep, close.face)
    }

    fn dummy_gregory_surface(index: usize) -> GregoryPatch {
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
        for (u, row) in grid.iter_mut().enumerate() {
            for (v, point) in row.iter_mut().enumerate() {
                *point = Point3::new(u as f64 / 9.0, v as f64 / 9.0, index as f64 * 1.0e-6);
            }
        }
        GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid }))
    }

    fn external_boundary_edge_count(
        brep: &BRep,
        cap_faces: &BTreeSet<FaceHandle>,
        face: FaceHandle,
    ) -> usize {
        face_outer_half_edges(brep, face)
            .expect("face loop")
            .into_iter()
            .filter(|half_edge| {
                let twin = brep.half_edges[*half_edge].twin;
                let Some(twin_face) = brep
                    .loops
                    .get(brep.half_edges[twin].loop_)
                    .map(|loop_| loop_.face)
                else {
                    return true;
                };
                !cap_faces.contains(&twin_face)
            })
            .count()
    }
}
