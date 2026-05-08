#![cfg(test)]
use super::tests_support::*;

#[test]
#[ignore]
fn cube_mixed_apex_baseline_no_reroute() {
    // Pre-reroute baseline: sphere-cap rung is on (caps are
    // Sphere) but Path B vertex-split surgery hasn't fired, so
    // BRep vertices still sit at rolling-ball footprint
    // positions and the cap-vs-cylinder metric is poor.
    unsafe {
        std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
        std::env::set_var("POPKERN_SPHERE_CAP_RUNG", "1");
        std::env::remove_var("POPKERN_VAR_RADIUS_REROUTE");
        std::env::remove_var("POPKERN_APEX_TOPOLOGY_REBUILD");
    }
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let apex_pre = cube_apex(&cube);
    let edges = corner_edges(&cube, apex_pre);
    let radii = [0.08_f64, 0.10, 0.12];
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let result = crate::fillet_edges(
        &cube,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("fillet cube_mixed");
    let brep = result.brep;
    brep.validate().expect("validates");
    let metrics = per_sphere_cap_metrics(&brep);
    eprintln!(
        "BASELINE cube_mixed (rung-on, NO reroute): {} Sphere cap faces",
        metrics.len()
    );
    for (fh, n, g0, g1) in &metrics {
        eprintln!("  {fh:?}: n_loop_edges={n} max_g0={g0:.3e} max_g1={g1:.3e}");
    }
}

/// Phase 6c-D commercial-grade gate for cube_mixed (variable-
/// radius fillets `[0.08, 0.10, 0.12]`).  Sphere-cap rung +
/// Path B vertex-split reroute together produce per-fillet
/// Sphere cap surfaces whose boundaries lie on the canonical
/// (sphere ∩ cylinder) tangent circles to machine precision.
///
/// Each cap sub-face has its own per-fillet `Sphere(centre_k,
/// r_k)` — analytic, perfect cap-vs-cylinder G1 by construction.
#[test]
#[ignore]
fn cube_mixed_apex_post_reroute_commercial_gate() {
    let brep = fillet_cube_mixed_with_rebuild(false);
    brep.validate().expect("post-reroute validates");
    let metrics = per_sphere_cap_metrics(&brep);
    eprintln!(
        "POST-REROUTE cube_mixed (rung-on): {} Sphere cap faces",
        metrics.len()
    );
    for (fh, n, g0, g1) in &metrics {
        eprintln!("  {fh:?}: n_loop_edges={n} max_g0={g0:.3e} max_g1={g1:.3e}");
    }
    assert_eq!(
        metrics.len(),
        3,
        "cube_mixed should produce 3 per-fillet Sphere caps"
    );
    // Commercial gate: machine-precision G0/G1 on cap-vs-cylinder
    // edges. The Sphere cap is analytic and its tangent circle
    // with each adjacent cylinder is a great circle; once Path B
    // moves the BRep vertices to the canonical sphere-fit
    // endpoints, the cap-tangent arcs lie EXACTLY on (sphere ∩
    // cylinder) and cap-vs-cylinder G0/G1 collapses to numerical
    // round-off.
    for (fh, _n, g0, g1) in &metrics {
        assert!(
            *g0 < 1e-10,
            "cube_mixed cap {fh:?}: G0 must be at machine precision; got {g0:.3e}"
        );
        assert!(
            *g1 < 1e-10,
            "cube_mixed cap {fh:?}: G1 must be at machine precision; got {g1:.3e}"
        );
    }
}
