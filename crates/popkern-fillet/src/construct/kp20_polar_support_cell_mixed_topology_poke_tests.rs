use super::*;
use crate::construct::kp20_polar_support_cell_inventory::{
    Kp20SupportCellInventory, Kp20SupportCellOwnerInventory, Kp20SupportCellOwnerTarget,
};
use popkern_core::geometry::gregory::{GregoryPatch, Kp20PieceId, RectangularGregoryBiNonic};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
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
fn owner_mixed_support_cell_commit_splits_multiple_support_owners() {
    let (mut brep, face) = build_flat_2ngon(6);
    let support_surfaces = (0..30).map(dummy_gregory_surface).collect::<Vec<_>>();
    let inventory = Kp20SupportCellOwnerInventory {
        sector_count: 6,
        targets: vec![
            Kp20SupportCellOwnerTarget {
                sector: 0,
                target_piece: Kp20PieceId::NearMidpointRight,
            },
            Kp20SupportCellOwnerTarget {
                sector: 2,
                target_piece: Kp20PieceId::NearMidpointLeft,
            },
        ],
        cells_per_target_sector: 4,
    };

    let commit = commit_kp20_owner_mixed_support_cell_patchwork(
        &mut brep,
        face,
        Point3::new(0.0, 0.0, 0.0),
        6,
        &support_surfaces,
        &inventory,
    )
    .expect("owner mixed support-cell commit");

    assert_eq!(commit.sub_faces.len(), 30);
    assert_eq!(commit.support_surface_indices.len(), 30);
    assert!(commit.sub_faces.iter().all(|face| matches!(
        brep.faces[*face].surface.as_ref().map(|surface| &surface.geometry),
        Some(SurfaceGeometry::Gregory(_))
    )));
    brep.validate()
        .expect("owner mixed support-cell commit leaves valid BRep");
}

#[test]
fn owner_mixed_support_cell_commit_splits_adjacent_support_owners_in_same_sector() {
    let (mut brep, face) = build_flat_2ngon(6);
    let support_surfaces = (0..30).map(dummy_gregory_surface).collect::<Vec<_>>();
    let inventory = Kp20SupportCellOwnerInventory {
        sector_count: 6,
        targets: vec![
            Kp20SupportCellOwnerTarget {
                sector: 0,
                target_piece: Kp20PieceId::NearMidpointRight,
            },
            Kp20SupportCellOwnerTarget {
                sector: 0,
                target_piece: Kp20PieceId::NearMidpointLeft,
            },
        ],
        cells_per_target_sector: 4,
    };

    let commit = commit_kp20_owner_mixed_support_cell_patchwork(
        &mut brep,
        face,
        Point3::new(0.0, 0.0, 0.0),
        6,
        &support_surfaces,
        &inventory,
    )
    .expect("same-sector owner mixed support-cell commit");

    assert_eq!(commit.sub_faces.len(), 30);
    assert_eq!(commit.support_surface_indices.len(), 30);
    assert!(commit.sub_faces.iter().all(|face| matches!(
        brep.faces[*face].surface.as_ref().map(|surface| &surface.geometry),
        Some(SurfaceGeometry::Gregory(_))
    )));
    brep.validate()
        .expect("same-sector owner mixed support-cell commit leaves valid BRep");
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
