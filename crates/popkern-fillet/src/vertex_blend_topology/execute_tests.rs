//! Phase 6c-C.2b executor tests.
//!
//! Drives the executor against cube_uniform's plan and verifies
//! the partial-execution invariants:
//! * 3 cap-tangent arc edges inserted
//! * 3 lens sub-faces created
//! * 3 footprint vertices remain (each still incident to its
//!   2 obsolete generator edges)
//! * BRep `validate()` still passes after the insertions

use super::execute::{ExecuteError, execute_deletion_plan};
use super::plan::plan_apex_region_deletion;
use crate::FilletOptions;
use crate::types::{FilletRadius, FilletSpec};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

fn corner_edges(brep: &BRep, vertex: VertexHandle) -> Vec<EdgeHandle> {
    use std::collections::BTreeSet;
    brep.half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn cube_apex(brep: &BRep) -> VertexHandle {
    brep.vertices
        .keys()
        .find(|&v| brep.half_edges.values().filter(|he| he.origin == v).count() == 3)
        .expect("cube has at least one valence-3 vertex")
}

fn fillet_cube_uniform() -> (BRep, FaceHandle) {
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
fn execute_inserts_three_cap_tangent_arcs() {
    let (mut brep, cap_face) = fillet_cube_uniform();
    let plan = plan_apex_region_deletion(&brep, cap_face);
    let report = execute_deletion_plan(&mut brep, &plan).expect("execution succeeds");
    assert_eq!(report.inserted_arcs.len(), 3);
    assert_eq!(report.lens_faces.len(), 3);
}

#[test]
fn execute_post_brep_passes_validate() {
    let (mut brep, cap_face) = fillet_cube_uniform();
    let plan = plan_apex_region_deletion(&brep, cap_face);
    let _ = execute_deletion_plan(&mut brep, &plan).expect("execution succeeds");
    // Post-execution BRep should still validate (the lens sub-
    // faces are valid topology even though they're geometrically
    // redundant — they share the cap surface with the real cap).
    brep.validate()
        .expect("BRep validates after partial execution");
}

#[test]
fn execute_rejects_unactionable_plan() {
    use super::plan::DeletionPlan;
    let (mut brep, _) = fillet_cube_uniform();
    let empty = DeletionPlan::default();
    let err = execute_deletion_plan(&mut brep, &empty).expect_err("empty plan rejected");
    assert_eq!(err, ExecuteError::PlanNotActionable);
}

#[test]
fn execute_inserted_arcs_lie_on_cap_sphere() {
    // Each inserted cap-tangent arc must trace a curve on the
    // cap sphere — at samples, the points must satisfy
    // |p - centre| ≈ radius within sphere_tol.
    use popkern_core::geometry::CurveEval;
    let (mut brep, cap_face) = fillet_cube_uniform();
    let plan = plan_apex_region_deletion(&brep, cap_face);
    let (centre, radius) = plan.cap_sphere.unwrap();
    let report = execute_deletion_plan(&mut brep, &plan).expect("execution succeeds");
    let sphere_tol = (radius * 0.01).max(1e-6);
    for &edge in &report.inserted_arcs {
        let coe = brep.edges[edge]
            .curve
            .as_ref()
            .expect("inserted edge has curve");
        let (t0, t1) = coe.geometry.parameter_range();
        for i in 0..=8 {
            let s = i as f64 / 8.0;
            let t = t0 + s * (t1 - t0);
            let p = coe.geometry.point_at(t);
            let off = (p.distance_to(&centre) - radius).abs();
            assert!(
                off < sphere_tol,
                "arc point at t={t:.3} is {off:.3e} off cap sphere (centre={centre:?}, r={radius})"
            );
        }
    }
}
