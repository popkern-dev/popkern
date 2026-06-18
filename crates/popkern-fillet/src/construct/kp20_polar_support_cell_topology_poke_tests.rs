use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryBiNonic};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::support_cell_topology_pattern::Kp20SupportCellTopologyPattern;
use super::{
    commit_kp20_support_cell_patchwork, commit_kp20_support_cell_patchwork_for_inventory,
    poke_kp20_support_cell_patchwork,
};
use super::{support_cell_poke_shape, support_cell_surface_slots};
use crate::construct::kp20_polar_expanded_support_binding::{
    attach_kp20_expanded_support_surfaces, author_kp20_expanded_support_pcurves,
    plan_kp20_expanded_support_surface_bindings,
};
use crate::construct::kp20_polar_support_cell_inventory::support_cell_inventory_from_provenance;
use crate::types::CornerSubPatchProvenance;

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

#[test]
fn support_cell_poke_emits_bindable_rectangular_support_cells() {
    let (mut brep, face) = build_flat_2ngon(6);
    let result = poke_kp20_support_cell_patchwork(
        &mut brep,
        face,
        Point3::new(0.0, 0.0, 0.0),
        &[0, 1, 2, 3, 4, 5],
    )
    .expect("support-cell poke");

    assert_eq!(result.support_cell_faces.len(), 24);
    assert_eq!(result.shape.face_count, 24);
    assert_eq!(
        result.shape.quad_count + result.shape.pentagon_count,
        result.shape.face_count
    );
    assert_eq!(result.surface_slots.len(), 24);

    let binding_plan = plan_kp20_expanded_support_surface_bindings(&result.surface_slots, 24)
        .expect("one-to-one cell binding plan");
    assert_eq!(binding_plan.bindings.len(), 24);
    assert_eq!(
        binding_plan.full_domain_binding_count,
        result.shape.quad_count
    );
    assert_eq!(
        binding_plan.trimmed_macro_binding_count,
        result.shape.pentagon_count
    );
    assert_eq!(binding_plan.reused_source_count, 0);

    let support_surfaces = (0..24).map(dummy_gregory_surface).collect::<Vec<_>>();
    attach_kp20_expanded_support_surfaces(&mut brep, &binding_plan, &support_surfaces)
        .expect("surface attachment");
    assert!(result.support_cell_faces.iter().all(|face| matches!(
        brep.faces[*face].surface.as_ref().map(|surface| &surface.geometry),
        Some(SurfaceGeometry::Gregory(_))
    )));
    author_kp20_expanded_support_pcurves(&mut brep, &binding_plan).expect("pcurve authoring");
    assert!(
        result
            .support_cell_faces
            .iter()
            .all(|face| outer_loop_has_pcurves(&brep, *face))
    );
    brep.validate()
        .expect("support-cell poke leaves valid BRep");
}

#[test]
fn batched_eight_cell_poke_emits_bindable_rectangular_support_cells() {
    let (mut brep, face) = build_flat_2ngon(6);
    let pattern = Kp20SupportCellTopologyPattern::eight_cell_support_grid();
    let macro_poke = crate::construct::poke_4n::poke_n_gon_face_into_expected_4n_quads(
        &mut brep,
        face,
        Point3::new(0.0, 0.0, 0.0),
        6,
    )
    .expect("macro support-cell poke");
    let support_faces = macro_poke
        .sub_faces
        .chunks_exact(4)
        .map(|sector| sector[3])
        .collect::<Vec<_>>();
    let support_faces = pattern
        .split_sector_support_faces(&mut brep, &support_faces)
        .expect("batched eight-cell support split");

    assert_eq!(macro_poke.sub_faces.len(), 24);
    assert_eq!(support_faces.len(), 48);
    let shape = support_cell_poke_shape(&brep, &support_faces).expect("support-cell shape");
    assert_eq!(shape.face_count, 48);
    assert_eq!(shape.quad_count + shape.pentagon_count, shape.face_count);
    let surface_slots = support_cell_surface_slots(&brep, &support_faces)
        .expect("eight-cell support surface slots");
    let binding_plan = plan_kp20_expanded_support_surface_bindings(&surface_slots, 48)
        .expect("one-to-one eight-cell binding plan");
    assert_eq!(binding_plan.bindings.len(), 48);
    assert_eq!(binding_plan.reused_source_count, 0);
    brep.validate()
        .expect("batched eight-cell support split leaves valid BRep");
}

#[test]
fn eight_cell_commit_attaches_surfaces_and_pcurves() {
    let (mut brep, face) = build_flat_2ngon(6);
    let support_surfaces = (0..48).map(dummy_gregory_surface).collect::<Vec<_>>();

    let commit = commit_kp20_support_cell_patchwork(
        &mut brep,
        face,
        Point3::new(0.0, 0.0, 0.0),
        6,
        &support_surfaces,
    )
    .expect("eight-cell support-cell commit");

    assert_eq!(commit.support_cell_faces.len(), 48);
    assert_eq!(commit.support_surface_indices.len(), 48);
    assert!(commit.support_cell_faces.iter().all(|face| matches!(
        brep.faces[*face].surface.as_ref().map(|surface| &surface.geometry),
        Some(SurfaceGeometry::Gregory(_))
    )));
    assert!(
        commit
            .support_cell_faces
            .iter()
            .all(|face| outer_loop_has_pcurves(&brep, *face))
    );
    assert_shared_edge_curves_match_topology_vertices(&brep, &commit.support_cell_faces);
    brep.validate()
        .expect("eight-cell support-cell commit leaves valid BRep");
}

#[test]
fn sixteen_cell_commit_attaches_surfaces_and_pcurves() {
    let (mut brep, face) = build_flat_2ngon(6);
    let support_surfaces = (0..96).map(dummy_gregory_surface).collect::<Vec<_>>();

    let commit = commit_kp20_support_cell_patchwork(
        &mut brep,
        face,
        Point3::new(0.0, 0.0, 0.0),
        6,
        &support_surfaces,
    )
    .expect("sixteen-cell support-cell commit");

    assert_eq!(commit.support_cell_faces.len(), 96);
    assert_eq!(commit.support_surface_indices.len(), 96);
    assert!(commit.support_cell_faces.iter().all(|face| matches!(
        brep.faces[*face].surface.as_ref().map(|surface| &surface.geometry),
        Some(SurfaceGeometry::Gregory(_))
    )));
    assert!(
        commit
            .support_cell_faces
            .iter()
            .all(|face| outer_loop_has_pcurves(&brep, *face))
    );
    assert_shared_edge_curves_match_topology_vertices(&brep, &commit.support_cell_faces);
    brep.validate()
        .expect("sixteen-cell support-cell commit leaves valid BRep");
}

#[test]
fn sparse_sector_commit_uses_original_sector_provenance() {
    let (mut brep, face) = build_flat_2ngon(6);
    let support_surfaces = (0..12).map(dummy_gregory_surface).collect::<Vec<_>>();
    let provenance = [0, 0, 0, 0, 2, 2, 2, 2, 5, 5, 5, 5]
        .into_iter()
        .enumerate()
        .map(|(index, sector)| CornerSubPatchProvenance::kp20_support_cell(sector, 6, index, 12))
        .collect::<Vec<_>>();
    let inventory = support_cell_inventory_from_provenance(6, support_surfaces.len(), &provenance)
        .expect("sparse support-cell inventory");

    let commit = commit_kp20_support_cell_patchwork_for_inventory(
        &mut brep,
        face,
        Point3::new(0.0, 0.0, 0.0),
        6,
        &support_surfaces,
        &inventory,
    )
    .expect("sparse support-cell commit");

    assert_eq!(commit.support_cell_faces.len(), 12);
    assert_eq!(commit.support_surface_indices.len(), 12);
    assert_eq!(inventory.target_sectors, vec![0, 2, 5]);
    assert!(commit.support_cell_faces.iter().all(|face| matches!(
        brep.faces[*face].surface.as_ref().map(|surface| &surface.geometry),
        Some(SurfaceGeometry::Gregory(_))
    )));
    assert!(
        commit
            .support_cell_faces
            .iter()
            .all(|face| outer_loop_has_pcurves(&brep, *face))
    );
    assert_shared_edge_curves_match_topology_vertices(&brep, &commit.support_cell_faces);
    brep.validate()
        .expect("sparse support-cell commit leaves valid BRep");
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

fn outer_loop_has_pcurves(brep: &BRep, face: FaceHandle) -> bool {
    let start = brep.loops[brep.faces[face].outer_loop].half_edge;
    let mut cursor = start;
    for _ in 0..brep.half_edges.len().max(1) {
        if brep.half_edges[cursor].pcurve.is_none() {
            return false;
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return true;
        }
    }
    false
}

fn assert_shared_edge_curves_match_topology_vertices(brep: &BRep, faces: &[FaceHandle]) {
    let mut checked = 0;
    for &face in faces {
        let start = brep.loops[brep.faces[face].outer_loop].half_edge;
        let mut cursor = start;
        for _ in 0..brep.half_edges.len().max(1) {
            let twin = brep.half_edges[cursor].twin;
            let other_face = brep.loops[brep.half_edges[twin].loop_].face;
            let edge = brep.half_edges[cursor].edge;
            if faces.contains(&other_face) && brep.edges[edge].half_edge == cursor {
                let canonical_twin = brep.half_edges[cursor].twin;
                let start_vertex = brep.half_edges[cursor].origin;
                let end_vertex = brep.half_edges[canonical_twin].origin;
                let curve = brep.edges[edge].curve.as_ref().expect("shared edge curve");
                assert!(
                    curve
                        .geometry
                        .point_at(curve.t_start)
                        .distance_to(&brep.vertices[start_vertex].position)
                        <= 1.0e-9
                );
                assert!(
                    curve
                        .geometry
                        .point_at(curve.t_end)
                        .distance_to(&brep.vertices[end_vertex].position)
                        <= 1.0e-9
                );
                checked += 1;
            }
            cursor = brep.half_edges[cursor].next;
            if cursor == start {
                break;
            }
        }
    }
    assert!(checked > 0, "expected support-cell shared edges");
}
