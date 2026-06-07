use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryBiNonic};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::poke_kp20_expanded_support_patchwork;
use crate::construct::kp20_polar_expanded_support_binding::{
    Kp20ExpandedSupportSurfaceAssignment, attach_kp20_expanded_support_surfaces,
    author_kp20_expanded_support_pcurves, plan_kp20_expanded_support_shared_trim_bindings,
    plan_kp20_expanded_support_surface_bindings,
};

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

fn loop_he_count(brep: &BRep, face: FaceHandle) -> usize {
    let start = brep.loops[brep.faces[face].outer_loop].half_edge;
    let mut cursor = start;
    let mut count = 0;
    for _ in 0..=brep.half_edges.len() {
        count += 1;
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return count;
        }
    }
    panic!("loop did not close");
}

#[test]
fn expanded_support_poke_splits_selected_support_faces() {
    let (mut brep, face) = build_flat_2ngon(6);
    let result = poke_kp20_expanded_support_patchwork(
        &mut brep,
        face,
        Point3::new(0.0, 0.0, 0.0),
        &[0, 1, 2, 3, 4, 5],
    )
    .expect("expanded support poke");

    assert_eq!(result.sub_faces.len(), 45);
    assert_eq!(result.shape.face_count, 45);
    assert_eq!(result.shape.quad_count, 39);
    assert_eq!(result.shape.pentagon_count, 6);
    assert!(result.shape.requires_trimmed_surface_assignment());
    assert_eq!(result.surface_slots.len(), result.sub_faces.len());
    for (slot, &sub_face) in result.surface_slots.iter().zip(result.sub_faces.iter()) {
        assert_eq!(slot.face, sub_face);
    }
    assert_eq!(
        result
            .surface_slots
            .iter()
            .filter(|slot| matches!(
                slot.assignment,
                Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace
            ))
            .count(),
        6
    );
    let current_patchwork = plan_kp20_expanded_support_surface_bindings(&result.surface_slots, 42);
    assert!(
        matches!(
            current_patchwork,
            Err(popkern_core::error::CoreError::PreconditionFailed(_))
        ),
        "42 rectangular support surfaces must not silently bind 45 topology slots",
    );
    let binding_plan = plan_kp20_expanded_support_surface_bindings(&result.surface_slots, 45)
        .expect("complete binding plan");
    assert_eq!(binding_plan.bindings.len(), result.sub_faces.len());
    assert_eq!(binding_plan.full_domain_binding_count, 39);
    assert_eq!(binding_plan.trimmed_macro_binding_count, 6);
    assert_eq!(binding_plan.reused_source_count, 0);
    let shared_binding_plan =
        plan_kp20_expanded_support_shared_trim_bindings(&result.surface_slots, 42)
            .expect("shared trimmed macro binding plan");
    assert_eq!(shared_binding_plan.bindings.len(), result.sub_faces.len());
    assert_eq!(shared_binding_plan.full_domain_binding_count, 39);
    assert_eq!(shared_binding_plan.trimmed_macro_binding_count, 6);
    assert_eq!(shared_binding_plan.reused_source_count, 3);
    assert!(
        shared_binding_plan
            .bindings
            .iter()
            .all(|binding| binding.support_surface_index < 42)
    );
    let support_surfaces = (0..42).map(dummy_gregory_surface).collect::<Vec<_>>();
    attach_kp20_expanded_support_surfaces(&mut brep, &shared_binding_plan, &support_surfaces)
        .expect("surface attachment");
    assert!(result.sub_faces.iter().all(|face| matches!(
        brep.faces[*face].surface.as_ref().map(|surface| &surface.geometry),
        Some(SurfaceGeometry::Gregory(_))
    )));
    author_kp20_expanded_support_pcurves(&mut brep, &shared_binding_plan)
        .expect("pcurve authoring");
    assert!(
        result
            .sub_faces
            .iter()
            .all(|face| outer_loop_has_pcurves(&brep, *face))
    );
    let loop_counts = result
        .sub_faces
        .iter()
        .map(|&sub_face| loop_he_count(&brep, sub_face))
        .collect::<Vec<_>>();
    assert!(loop_counts.iter().all(|count| (4..=5).contains(count)));
    assert_eq!(loop_counts.iter().filter(|count| **count == 5).count(), 6);
    brep.validate()
        .expect("expanded support poke leaves valid BRep");
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
