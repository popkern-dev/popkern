//! Phase 6c-C.2 plan-builder tests.
//!
//! Verifies `plan_apex_region_deletion` against the
//! cube_uniform fixture's known structure: 6 cap-loop vertices
//! split 3 sphere-fit + 3 rolling-ball footprints, 6 generator
//! edges to delete, 3 cap-tangent arcs to insert.

use super::plan::plan_apex_region_deletion;
use crate::FilletOptions;
use crate::types::{FilletRadius, FilletSpec};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

fn corner_edges(brep: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    use std::collections::BTreeSet;
    let edges: BTreeSet<EdgeHandle> = brep
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    edges.into_iter().collect()
}

fn cube_apex(brep: &BRep) -> VertexHandle {
    brep.vertices
        .keys()
        .find(|&v| brep.half_edges.values().filter(|he| he.origin == v).count() == 3)
        .expect("cube has at least one valence-3 vertex")
}

fn fillet_cube_uniform() -> (BRep, popkern_core::topology::FaceHandle) {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let apex_pre = cube_apex(&cube);
    let edges = corner_edges(&cube, apex_pre);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.10),
        })
        .collect();
    let result = crate::fillet_edges(
        &cube,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("fillet uniform cube");
    let cap_face = result.corner_patches[0].face;
    (result.brep, cap_face)
}

#[test]
fn plan_cube_uniform_keeps_three_sphere_fit_vertices() {
    let (brep, cap_face) = fillet_cube_uniform();
    let plan = plan_apex_region_deletion(&brep, cap_face);
    assert_eq!(plan.cap_face, Some(cap_face));
    assert_eq!(
        plan.keep_vertices.len(),
        3,
        "uniform cube cap has exactly 3 sphere-fit endpoints"
    );
    // Each kept vertex must lie ON the cap sphere (within 1%
    // of cap radius).
    let (centre, radius) = plan.cap_sphere.expect("plan has cap sphere");
    for &v in &plan.keep_vertices {
        let p = brep.vertices[v].position;
        let off = (p.distance_to(&centre) - radius).abs();
        assert!(
            off < radius * 0.01,
            "kept vertex {v:?} at {p:?} is {off} off the cap sphere (centre={centre:?}, r={radius})"
        );
    }
}

#[test]
#[ignore = "route B corner reconstruction (contact collapse + exact boundary arcs) removed the legacy anchor-hexagon cap anatomy this machinery was built to clean up; the apex-region deletion plan finds no redundant lens structures on the reconstructed corner. Retire or repurpose with PK-101."]
fn plan_cube_uniform_deletes_three_footprints() {
    let (brep, cap_face) = fillet_cube_uniform();
    let plan = plan_apex_region_deletion(&brep, cap_face);
    assert_eq!(
        plan.delete_vertices.len(),
        3,
        "uniform cube cap has exactly 3 rolling-ball-footprint vertices to delete"
    );
    // Each deleted vertex must NOT lie on the cap sphere
    // (it's the off-sphere case).
    let (centre, radius) = plan.cap_sphere.unwrap();
    for &v in &plan.delete_vertices {
        let p = brep.vertices[v].position;
        let off = (p.distance_to(&centre) - radius).abs();
        assert!(
            off >= radius * 0.01,
            "delete-marked vertex {v:?} at {p:?} is unexpectedly *on* the cap sphere"
        );
    }
}

#[test]
#[ignore = "route B corner reconstruction (contact collapse + exact boundary arcs) removed the legacy anchor-hexagon cap anatomy this machinery was built to clean up; the apex-region deletion plan finds no redundant lens structures on the reconstructed corner. Retire or repurpose with PK-101."]
fn plan_cube_uniform_deletes_six_cap_cylinder_generators() {
    let (brep, cap_face) = fillet_cube_uniform();
    let plan = plan_apex_region_deletion(&brep, cap_face);
    assert_eq!(
        plan.delete_cap_cylinder_generators.len(),
        6,
        "cube_uniform cap loop has 6 generator edges (2 per fillet cylinder)"
    );
}

#[test]
fn plan_cube_uniform_inserts_three_cap_tangent_arcs() {
    let (brep, cap_face) = fillet_cube_uniform();
    let plan = plan_apex_region_deletion(&brep, cap_face);
    assert_eq!(
        plan.insert_cap_tangent_arcs.len(),
        3,
        "uniform cube post-rebuild cap has 3 cap-tangent arc edges"
    );
    // Each insertion connects two distinct sphere-fit
    // (kept) vertices.
    for &(a, b, _) in &plan.insert_cap_tangent_arcs {
        assert_ne!(a, b, "arc endpoints must be distinct");
        assert!(plan.keep_vertices.contains(&a));
        assert!(plan.keep_vertices.contains(&b));
    }
}

#[test]
#[ignore = "route B corner reconstruction (contact collapse + exact boundary arcs) removed the legacy anchor-hexagon cap anatomy this machinery was built to clean up; the apex-region deletion plan finds no redundant lens structures on the reconstructed corner. Retire or repurpose with PK-101."]
fn plan_cube_uniform_is_actionable() {
    let (brep, cap_face) = fillet_cube_uniform();
    let plan = plan_apex_region_deletion(&brep, cap_face);
    assert!(
        plan.is_actionable(),
        "uniform cube plan should be actionable: {plan:#?}"
    );
}

#[test]
fn plan_unfilleted_cube_is_empty() {
    // A face that *isn't* a cap face shouldn't yield an
    // actionable plan — in particular, a pre-fillet cube face
    // has no cap sphere assigned.
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let any_face = cube.faces.keys().next().unwrap();
    let plan = plan_apex_region_deletion(&cube, any_face);
    assert!(plan.cap_face.is_none(), "non-cap face yields no plan");
    assert!(!plan.is_actionable());
}
