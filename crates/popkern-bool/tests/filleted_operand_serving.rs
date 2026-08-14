//! Serving-quality gate for booleans on filleted operands (PK-272).
//!
//! Promoted from a PK-272 diagnostic probe: the pre-postprocess
//! `validate_all` diagnostic reported a winding-dishonest cylinder
//! wall on this exact union, and the probe established that the
//! finding was transient — `reconcile_orientation_flags` in
//! `run_output_postprocess` repairs it before the result is served.
//! This test pins that contract: a boolean whose operand carries a
//! periodic torus fillet strip must serve a `validate_all`-clean
//! result (winding honesty included), not merely pass `post_check`.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[test]
fn served_torus_union_winding_honesty() {
    let tol = Tolerance::default();
    let cyl = BRep::cylinder(1.0, 2.0, tol).expect("cylinder builds");
    let base_edge = cyl
        .edges
        .keys()
        .find(|&e| {
            let he = cyl.edges[e].half_edge;
            let twin = cyl.half_edges[he].twin;
            let f1 = cyl.loops[cyl.half_edges[he].loop_].face;
            let f2 = cyl.loops[cyl.half_edges[twin].loop_].face;
            if f1 == f2 {
                return false;
            }
            let a = cyl.vertices[cyl.half_edges[he].origin].position;
            let b = cyl.vertices[cyl.half_edges[twin].origin].position;
            a.z.abs() < 1e-9 && b.z.abs() < 1e-9
        })
        .expect("cylinder has a base edge");
    let specs = vec![FilletSpec {
        edge: base_edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let filleted = fillet_edges(&cyl, &specs, &FilletOptions::default(), tol, None)
        .expect("cylinder base fillet succeeds")
        .brep;
    let profile = [
        Point3::new(0.5, -0.25, 0.75),
        Point3::new(1.5, -0.25, 0.75),
        Point3::new(1.5, 0.25, 0.75),
        Point3::new(0.5, 0.25, 0.75),
    ];
    let boss = BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), 0.5, tol).expect("boss builds");
    let served = popkern_bool::union_ref(&filleted, &boss, tol, None).expect("union succeeds");
    let errs = served.validate_all();
    for e in &errs {
        eprintln!("[serving] validate_all issue: {e}");
    }
    assert!(
        errs.is_empty(),
        "served torus-carrying union result fails validate_all with {} issues",
        errs.len()
    );
}
