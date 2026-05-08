use super::tests_support::*;

#[test]
#[ignore]
fn cube_uniform_apex_no_rebuild_baseline() {
    let (brep, cap_face) = fillet_uniform_with_rebuild(false);
    brep.validate().expect("baseline BRep validates");
    let n_edges = cap_loop_edge_count(&brep, cap_face);
    let (g0, g1) = cap_vs_cylinder_g0_g1(&brep, cap_face);
    let surf_kind = brep.faces[cap_face]
        .surface
        .as_ref()
        .map(|s| match &s.geometry {
            SurfaceGeometry::Sphere { .. } => "Sphere",
            SurfaceGeometry::Gregory(_) => "Gregory",
            SurfaceGeometry::Plane { .. } => "Plane",
            SurfaceGeometry::Cylinder { .. } => "Cyl",
            SurfaceGeometry::Nurbs(_) => "Nurbs",
            _ => "Other",
        })
        .unwrap_or("?");
    eprintln!(
        "BASELINE (no rebuild): cap_face={cap_face:?} surf={surf_kind} \
         cap_loop_edges={n_edges} \
         BRep V={} E={} F={}, max_g0={g0:.3e} max_g1={g1:.3e}",
        brep.vertices.len(),
        brep.edges.len(),
        brep.faces.len(),
    );
}

#[test]
#[ignore]
fn cube_uniform_apex_rebuild_e2e() {
    let (brep, cap_face) = fillet_uniform_with_rebuild(true);
    brep.validate().expect("post-rebuild BRep validates");
    let n_edges = cap_loop_edge_count(&brep, cap_face);
    let (g0, g1) = cap_vs_cylinder_g0_g1(&brep, cap_face);
    eprintln!(
        "REBUILD ON: cap_loop_edges={n_edges} \
         max_g0={g0:.3e} max_g1={g1:.3e}"
    );
    // Phase 6c-C topology gate: cap loop reduces 6 → 3 edges
    // (the sphere-fit triangle).  This is the hard guarantee.
    assert_eq!(n_edges, 3, "Phase 6c-C: cap loop should reduce 6 → 3 edges");
    // Phase 6c-C.2e G1 commercial gate: 5e-2. cube_uniform's 3
    // cap-tangent arcs hit ≤ 1e-3 G1 deviation after the
    // degree-9 NURBS pcurve fit + atan2 unwrap + Sphere u-period
    // clamp-skip. The residual is dominated by the pole-crossing
    // F_x arc; without that arc, all G1 values are at machine
    // precision.
    assert!(
        g1 < 5e-2,
        "Phase 6c-C: cap-vs-cylinder G1 deviation should pass commercial gate (5e-2); got {g1:.3e}"
    );
    // G0 commercial gate: 1e-2 (1 cm on a 1 m unit cube). Today's
    // residual sits ~4 mm — driven by the pole-crossing arc's
    // degree-9 fit. Tightening to machine precision needs an
    // analytic small-circle pcurve construction (deferred).
    assert!(
        g0 < 1e-2,
        "Phase 6c-C: cap-vs-cylinder G0 should pass commercial gate (1e-2); got {g0:.3e}"
    );
}
