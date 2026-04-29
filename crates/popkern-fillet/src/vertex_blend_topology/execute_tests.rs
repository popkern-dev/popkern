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
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

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

/// Phase 6c-C.2e probe — measure how close
/// `cap_surf(pcurve_a(t))` is to `arc_3d(t)` at samples along
/// each inserted arc. If pcurves are perfectly fit, the
/// distance should be at machine precision (1e-12). Larger
/// values reveal fit-residual issues.
#[test]
#[ignore]
fn pcurve_fit_residual_diagnostic() {
    use super::collapse::collapse_apex_lens;
    use popkern_core::geometry::CurveEval;
    use popkern_core::geometry::pcurve::PcurveGeometry;
    use popkern_core::geometry::surface::SurfaceEval;
    let (mut brep, cap_face) = fillet_cube_uniform();
    let plan = plan_apex_region_deletion(&brep, cap_face);
    let report = execute_deletion_plan(&mut brep, &plan).expect("execution succeeds");
    // Collapse each lens so pb's loop-face is the cylinder, not the
    // lens-copy of the sphere.
    for &lens in &report.lens_faces {
        let lens_loop = brep.faces[lens].outer_loop;
        let start_he = brep.loops[lens_loop].half_edge;
        let mut gens: Vec<HalfEdgeHandle> = Vec::new();
        let mut cursor = start_he;
        for _ in 0..3 {
            let he_edge = brep.half_edges[cursor].edge;
            let he_next = brep.half_edges[cursor].next;
            if !report.inserted_arcs.contains(&he_edge) {
                gens.push(cursor);
            }
            cursor = he_next;
            if cursor == start_he {
                break;
            }
        }
        let _ = collapse_apex_lens(&mut brep, lens, gens[0], gens[1]);
    }
    let cap_surf = brep.faces[cap_face]
        .surface
        .as_ref()
        .map(|s| s.geometry.clone())
        .expect("cap has surface");
    eprintln!("\n=== pcurve_fit_residual_diagnostic ===");
    eprintln!(
        "cap_surf={:?}",
        match &cap_surf {
            popkern_core::geometry::surface::SurfaceGeometry::Sphere { center, radius } => {
                format!("Sphere(centre={center:?}, r={radius})")
            }
            other => format!("{other:?}"),
        }
    );
    for (arc_idx, &edge) in report.inserted_arcs.iter().enumerate() {
        let coe = brep.edges[edge].curve.as_ref().unwrap().clone();
        let (t0, t1) = (coe.t_start, coe.t_end);
        // Find the cap-side half-edge for this edge.
        let cap_he = brep
            .half_edges
            .iter()
            .find(|(_, he)| he.edge == edge && brep.loops[he.loop_].face == cap_face)
            .map(|(h, _)| h)
            .expect("cap-side half-edge");
        let twin_he = brep.half_edges[cap_he].twin;
        let pa = brep.half_edges[cap_he].pcurve.clone().expect("cap pcurve");
        let pb = brep.half_edges[twin_he].pcurve.clone().expect("cyl pcurve");
        let twin_face = brep.loops[brep.half_edges[twin_he].loop_].face;
        let cyl_surf = brep.faces[twin_face]
            .surface
            .as_ref()
            .map(|s| s.geometry.clone())
            .expect("cyl surface");
        eprintln!("\narc[{arc_idx}] edge={edge:?}");
        match &pa.geometry {
            PcurveGeometry::Nurbs(nc) => {
                let cps: Vec<String> = nc
                    .control_points
                    .iter()
                    .map(|p| format!("({:.4},{:.4})", p.x, p.y))
                    .collect();
                eprintln!("  pa(cap) Nurbs deg={} cps=[{}]", nc.degree, cps.join(", "));
            }
            other => eprintln!("  pa(cap) {other:?}"),
        }
        match &pb.geometry {
            PcurveGeometry::Nurbs(nc) => {
                let cps: Vec<String> = nc
                    .control_points
                    .iter()
                    .map(|p| format!("({:.4},{:.4})", p.x, p.y))
                    .collect();
                eprintln!("  pb(cyl) Nurbs deg={} cps=[{}]", nc.degree, cps.join(", "));
            }
            other => eprintln!("  pb(cyl) {other:?}"),
        }
        let mut max_pa_resid = 0.0_f64;
        let mut max_pb_resid = 0.0_f64;
        for i in 0..=10 {
            let s = i as f64 / 10.0;
            let t = t0 + s * (t1 - t0);
            let p3 = coe.geometry.point_at(t);
            let (pa_lo, pa_hi) = pa.parameter_range();
            let (pb_lo, pb_hi) = pb.parameter_range();
            let uv_a = pa.point_at(pa_lo + s * (pa_hi - pa_lo));
            let uv_b = pb.point_at(pb_lo + s * (pb_hi - pb_lo));
            let p_a = cap_surf.point_at(uv_a.x, uv_a.y);
            let p_b = cyl_surf.point_at(uv_b.x, uv_b.y);
            let r_a = p3.distance_to(&p_a);
            let r_b = p3.distance_to(&p_b);
            if r_a > max_pa_resid {
                max_pa_resid = r_a;
            }
            if r_b > max_pb_resid {
                max_pb_resid = r_b;
            }
        }
        eprintln!(
            "  max |arc(t) - cap(pa(t))| = {max_pa_resid:.3e}\n  \
             max |arc(t) - cyl(pb(t))| = {max_pb_resid:.3e}"
        );
    }
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
