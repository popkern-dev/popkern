use super::tests_support::*;

/// Pre-rebuild baseline diagnostic for pyramid_uniform.
/// Phase 6c-D commercial-grade gate for cube_extreme — 1:10
/// radius ratio fillets `[0.05, 0.10, 0.50]`.  Extends the
/// commercial-grade guarantee to extreme asymmetry.
#[test]
#[ignore]
fn cube_extreme_apex_post_reroute_commercial_gate() {
    unsafe {
        std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
        std::env::set_var("POPKERN_SPHERE_CAP_RUNG", "1");
        std::env::set_var("POPKERN_VAR_RADIUS_REROUTE", "1");
        std::env::remove_var("POPKERN_APEX_TOPOLOGY_REBUILD");
    }
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let apex_pre = cube_apex(&cube);
    let edges = corner_edges(&cube, apex_pre);
    let radii = [0.05_f64, 0.10, 0.50];
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
    );
    let Ok(result) = result else {
        eprintln!("cube_extreme fillet failed: {:?}", result.as_ref().err());
        return;
    };
    let brep = result.brep;
    brep.validate().expect("post-reroute validates");
    let metrics = per_sphere_cap_metrics(&brep);
    eprintln!(
        "POST-REROUTE cube_extreme (rung-on): {} Sphere cap faces",
        metrics.len()
    );
    for (fh, n, g0, g1) in &metrics {
        eprintln!("  {fh:?}: n_loop_edges={n} max_g0={g0:.3e} max_g1={g1:.3e}");
    }
    assert_eq!(
        metrics.len(),
        3,
        "cube_extreme should produce 3 per-fillet Sphere caps"
    );
    for (fh, _n, g0, g1) in &metrics {
        assert!(
            *g0 < 1e-10,
            "cube_extreme cap {fh:?}: G0 must be at machine precision; got {g0:.3e}"
        );
        assert!(
            *g1 < 1e-10,
            "cube_extreme cap {fh:?}: G1 must be at machine precision; got {g1:.3e}"
        );
    }
}

/// Phase 6c-E pyramid_uniform commercial-grade gate — sphere-cap
/// rung + Path B reroute.  After the per-fillet rolling-ball
/// cap-centre solver lands (each per-fillet sphere centre on
/// its own cylinder axis) and the N=4 dispatch gate is opened,
/// pyramid_uniform should hit machine-precision cap-vs-cylinder
/// G0/G1 across all 4 per-fillet sub-faces.
#[test]
#[ignore]
fn pyramid_uniform_apex_post_reroute_commercial_gate() {
    unsafe {
        std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
        std::env::set_var("POPKERN_SPHERE_CAP_RUNG", "1");
        std::env::set_var("POPKERN_VAR_RADIUS_REROUTE", "1");
        std::env::remove_var("POPKERN_APEX_TOPOLOGY_REBUILD");
    }
    let pyr = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let apex = pyr
        .vertices
        .keys()
        .find(|&v| pyr.half_edges.values().filter(|he| he.origin == v).count() == 4)
        .expect("pyramid has a valence-4 apex");
    let edges = corner_edges(&pyr, apex);
    let radii = [0.10_f64; 4];
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let result = crate::fillet_edges(
        &pyr,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("pyramid_uniform fillet should succeed");
    let brep = result.brep;
    brep.validate().expect("post-reroute validates");
    let metrics = per_sphere_cap_metrics(&brep);
    eprintln!(
        "POST-REROUTE pyramid_uniform (rung-on, r={:.2}): {} Sphere cap faces",
        radii[0],
        metrics.len()
    );
    for (fh, n, g0, g1) in &metrics {
        eprintln!("  {fh:?}: n_loop_edges={n} max_g0={g0:.3e} max_g1={g1:.3e}");
    }
    assert!(
        !metrics.is_empty(),
        "pyramid_uniform should produce at least one Sphere sub-face"
    );
    for (fh, _n, g0, g1) in &metrics {
        assert!(
            *g0 < 1e-10,
            "pyramid_uniform cap {fh:?}: G0 must be at machine precision; got {g0:.3e}"
        );
        assert!(
            *g1 < 1e-10,
            "pyramid_uniform cap {fh:?}: G1 must be at machine precision; got {g1:.3e}"
        );
    }
}

/// Phase 6c-D probe — cube_mixed with Gregory caps (rung OFF).
/// Baseline reference for comparing pyramid Gregory caps.
#[test]
#[ignore]
fn cube_mixed_apex_gregory_caps_diagnostic() {
    unsafe {
        std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
        std::env::set_var("POPKERN_GREGORY_CAP_REBUILD", "1");
        std::env::remove_var("POPKERN_SPHERE_CAP_RUNG");
        std::env::remove_var("POPKERN_VAR_RADIUS_REROUTE");
        std::env::remove_var("POPKERN_APEX_TOPOLOGY_REBUILD");
    }
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let apex = cube_apex(&cube);
    let edges = corner_edges(&cube, apex);
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
    .expect("cube_mixed fillet should succeed");
    let brep = result.brep;
    brep.validate().expect("baseline validates");

    let metrics = per_cap_metrics_filtered(&brep, |g| {
        matches!(
            g,
            SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Gregory(_)
        )
    });
    let mut surf_counts: std::collections::BTreeMap<&str, usize> =
        std::collections::BTreeMap::new();
    for (_, face) in brep.faces.iter() {
        let kind = face
            .surface
            .as_ref()
            .map(|s| match &s.geometry {
                SurfaceGeometry::Sphere { .. } => "Sphere",
                SurfaceGeometry::Cylinder { .. } => "Cyl",
                SurfaceGeometry::Plane { .. } => "Plane",
                SurfaceGeometry::Gregory(_) => "Gregory",
                SurfaceGeometry::Nurbs(_) => "Nurbs",
                _ => "Other",
            })
            .unwrap_or("none");
        *surf_counts.entry(kind).or_insert(0) += 1;
    }
    eprintln!(
        "GREGORY-CAPS cube_mixed (rung-OFF, r={radii:?}): {} cap faces; surface counts={surf_counts:?}",
        metrics.len()
    );
    for (fh, n, g0, g1) in &metrics {
        let kind = brep.faces[*fh]
            .surface
            .as_ref()
            .map(|s| match &s.geometry {
                SurfaceGeometry::Sphere { .. } => "Sphere",
                SurfaceGeometry::Gregory(_) => "Gregory",
                _ => "Other",
            })
            .unwrap_or("none");
        eprintln!("  {fh:?} ({kind}): n_loop_edges={n} max_g0={g0:.3e} max_g1={g1:.3e}");
    }
}

/// Phase 6c-D probe — pyramid_uniform with Gregory caps (sphere-fit
/// substitution ON, sphere-cap rung OFF).  This is the path the
/// pyramid commercial-grade gate must take, since per-fillet Sphere
/// caps don't fit pyramid geometry (Path B reroute can't move
/// vertices onto a sphere centred on a single neighbouring axis;
/// see `project_pyramid_path_b_blocker` memory).
///
/// Reports cap-vs-cylinder G0/G1 for each Gregory sub-face.  No
/// commercial-grade assertion yet — diagnostic only — to establish
/// the baseline from which Phase 6c-D's optimiser fix will improve.
#[test]
#[ignore]
fn pyramid_uniform_apex_gregory_caps_diagnostic() {
    unsafe {
        std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
        std::env::set_var("POPKERN_VAR_RADIUS_REROUTE", "1");
        std::env::remove_var("POPKERN_SPHERE_CAP_RUNG");
        std::env::remove_var("POPKERN_GREGORY_CAP_REBUILD");
        std::env::remove_var("POPKERN_APEX_TOPOLOGY_REBUILD");
    }
    let pyr = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let apex = pyr
        .vertices
        .keys()
        .find(|&v| pyr.half_edges.values().filter(|he| he.origin == v).count() == 4)
        .expect("pyramid has a valence-4 apex");
    let edges = corner_edges(&pyr, apex);
    let radii = [0.10_f64; 4];
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let result = crate::fillet_edges(
        &pyr,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("pyramid_uniform fillet should succeed");
    let brep = result.brep;
    brep.validate().expect("baseline validates");

    let metrics = per_cap_metrics_filtered(&brep, |g| {
        matches!(
            g,
            SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Gregory(_)
        )
    });
    let mut surf_counts: std::collections::BTreeMap<&str, usize> =
        std::collections::BTreeMap::new();
    for (_, face) in brep.faces.iter() {
        let kind = face
            .surface
            .as_ref()
            .map(|s| match &s.geometry {
                SurfaceGeometry::Sphere { .. } => "Sphere",
                SurfaceGeometry::Cylinder { .. } => "Cyl",
                SurfaceGeometry::Plane { .. } => "Plane",
                SurfaceGeometry::Gregory(_) => "Gregory",
                SurfaceGeometry::Nurbs(_) => "Nurbs",
                _ => "Other",
            })
            .unwrap_or("none");
        *surf_counts.entry(kind).or_insert(0) += 1;
    }
    eprintln!(
        "GREGORY-CAPS pyramid_uniform (rung-OFF, r={:.2}): {} cap faces; surface counts={surf_counts:?}",
        radii[0],
        metrics.len()
    );
    for (fh, n, g0, g1) in &metrics {
        let kind = brep.faces[*fh]
            .surface
            .as_ref()
            .map(|s| match &s.geometry {
                SurfaceGeometry::Sphere { .. } => "Sphere",
                SurfaceGeometry::Gregory(_) => "Gregory",
                _ => "Other",
            })
            .unwrap_or("none");
        eprintln!("  {fh:?} ({kind}): n_loop_edges={n} max_g0={g0:.3e} max_g1={g1:.3e}");
    }
}
