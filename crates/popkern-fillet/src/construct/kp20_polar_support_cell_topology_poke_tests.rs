use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregoryBiNonic};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::poke_kp20_support_cell_patchwork;
use crate::construct::kp20_polar_expanded_support_binding::{
    attach_kp20_expanded_support_surfaces, author_kp20_expanded_support_pcurves,
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
