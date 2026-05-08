#![cfg(test)]
use super::tests_support::*;

/// Phase 6c-E pyramid baseline (default pipeline, no opt-in
/// flags).  Documents the gap metrics — pyramid's default path
/// today produces a single Sphere cap with a wrong centre.
#[test]
#[ignore]
fn pyramid_uniform_apex_baseline_metrics() {
    unsafe {
        std::env::remove_var("POPKERN_APEX_TOPOLOGY_REBUILD");
        std::env::remove_var("POPKERN_SPHERE_FIT_BOUNDARY");
        std::env::remove_var("POPKERN_SPHERE_CAP_RUNG");
        std::env::remove_var("POPKERN_VAR_RADIUS_REROUTE");
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
    let metrics = per_sphere_cap_metrics(&brep);
    eprintln!(
        "BASELINE pyramid_uniform (r={:.2}): {} Sphere cap faces",
        radii[0],
        metrics.len()
    );
    for (fh, n, g0, g1) in &metrics {
        eprintln!("  {fh:?}: n_loop_edges={n} max_g0={g0:.3e} max_g1={g1:.3e}");
    }
    // Dump all face surface kinds to see the post-dispatch state.
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
    eprintln!("Surface counts: {surf_counts:?}");
}

/// Phase 6c-E pyramid_uniform diagnostic — documents the gap.
///
/// **Status: known-broken upstream.** Pyramid_uniform's cap face
/// is assigned a single `Sphere` surface, but:
///
/// * The Sphere centre computed by `try_sphere_for_corner` uses
///   a `corner ± r` heuristic that's correct only for
///   axis-aligned trihedral apices (cube). For pyramid the
///   symmetric mean defeats the directional tiebreak and the
///   centre lands at `(-r, -r, h-r)` instead of on the symmetry
///   axis.
/// * Even with a correct symmetric-axis sphere centre `(0, 0,
///   z*)`, the cap-tangent geometry doesn't admit a single
///   sphere — the four slant cylinders have axes that don't
///   meet at one point with sphere centre (sphere centre would
///   be at `(0, 0, h+r)` outside the solid). The cap is
///   geometrically over-constrained for a single Sphere
///   primitive.
/// * Boundary vertices on the cap face aren't on the assigned
///   sphere AND aren't on the adjacent fillet cylinders — the
///   pre-cap pipeline placed them at intermediate positions.
///
/// **Path forward (deferred):** route pyramid_uniform through
/// the rational-Gregory N=4 sub-face split with sphere-cap rung,
/// analogous to cube_mixed.  Today the rational-Gregory dispatch
/// doesn't fire for pyramid_uniform (root cause TBD).
#[test]
#[ignore]
fn pyramid_uniform_planner_diagnostic() {
    use super::super::plan::plan_apex_region_deletion;
    unsafe {
        std::env::remove_var("POPKERN_APEX_TOPOLOGY_REBUILD");
        std::env::remove_var("POPKERN_SPHERE_FIT_BOUNDARY");
        std::env::remove_var("POPKERN_SPHERE_CAP_RUNG");
        std::env::remove_var("POPKERN_VAR_RADIUS_REROUTE");
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
    let cap_face = result.corner_patches[0].face;
    let plan = plan_apex_region_deletion(&brep, cap_face);
    eprintln!(
        "pyramid_uniform plan: cap_sphere={:?} keep={} delete={} arcs={} actionable={}",
        plan.cap_sphere,
        plan.keep_vertices.len(),
        plan.delete_vertices.len(),
        plan.insert_cap_tangent_arcs.len(),
        plan.is_actionable(),
    );
    // Phase 6c-C planner generalisation diagnostic — dump the cap
    // face's outer loop structure so we can see what the planner
    // is working against. For pyramid_uniform with sphere-cap-rung:
    // each cap face is a per-fillet Sphere; the outer loop has
    // 6 vertices (5 outer + 1 central) and 6 edges. Some neighbours
    // are cylinders (cap-tangent edges), some are spheres (cap-cap
    // spokes), some are planes (slant-bridge fragments). The
    // planner needs to identify which are cap-vs-cylinder and
    // construct fresh sphere∩cylinder tangent points there.
    if let Some(face) = brep.faces.get(cap_face) {
        let outer = face.outer_loop;
        let start = brep.loops[outer].half_edge;
        eprintln!(
            "  cap face outer loop (cap surface = {:?}):",
            face.surface.as_ref().map(|s| match &s.geometry {
                SurfaceGeometry::Sphere { center, radius } =>
                    format!("Sphere(c={center:?}, r={radius:.3})"),
                SurfaceGeometry::Gregory(_) => "Gregory".to_string(),
                _ => "Other".to_string(),
            })
        );
        let mut cur = start;
        for step in 0..16 {
            let he = &brep.half_edges[cur];
            let v_origin = brep.vertices[he.origin].position;
            let twin_face = brep.loops[brep.half_edges[he.twin].loop_].face;
            let twin_kind = brep
                .faces
                .get(twin_face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| match &s.geometry {
                    SurfaceGeometry::Sphere { .. } => "Sphere",
                    SurfaceGeometry::Cylinder { .. } => "Cyl",
                    SurfaceGeometry::Plane { .. } => "Plane",
                    SurfaceGeometry::Gregory(_) => "Gregory",
                    _ => "?",
                })
                .unwrap_or("none");
            // Distance from vertex to cap sphere (signed).
            let on_sphere_signed = if let Some((c, r)) = plan.cap_sphere {
                v_origin.distance_to(&c) - r
            } else {
                f64::NAN
            };
            eprintln!(
                "    step{step}: V({:.4},{:.4},{:.4})  twin={twin_kind} (face {twin_face:?})  dist_to_sphere={on_sphere_signed:+.3e}",
                v_origin.x, v_origin.y, v_origin.z
            );
            cur = he.next;
            if cur == start {
                break;
            }
        }
    }
}

/// Phase 6c-C #172 — dump the per-fillet sphere cap loop
/// structure under `POPKERN_SPHERE_CAP_RUNG=1`, which splits
/// the apex into 4 separate cap faces (one per fillet).
/// Each cap is a Sphere primitive; this diagnostic walks each
/// cap's outer loop and prints vertex positions, twin surface
/// kinds, and dist-to-this-cap-sphere signed values.
///
/// Used to design the planner generalisation that constructs
/// fresh on-sphere vertices for pyramid (where no existing
/// vertex lies on the per-fillet sphere).
#[test]
#[ignore]
fn pyramid_uniform_per_fillet_cap_loop_dump() {
    use super::super::plan::plan_apex_region_deletion;
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
    let metrics = per_sphere_cap_metrics(&brep);
    eprintln!(
        "\n=== Per-fillet sphere cap loop dump (rung ON) ===\n{} cap faces total",
        metrics.len()
    );
    for (cap_face_idx, (fh, n_edges, g0, g1)) in metrics.iter().enumerate() {
        eprintln!(
            "\nCap face #{cap_face_idx} ({fh:?}): n_loop_edges={n_edges} g0={g0:.3e} g1={g1:.3e}"
        );
        let plan = plan_apex_region_deletion(&brep, *fh);
        eprintln!(
            "  plan: cap_sphere={:?}  keep={} delete={} arcs={} actionable={}",
            plan.cap_sphere,
            plan.keep_vertices.len(),
            plan.delete_vertices.len(),
            plan.insert_cap_tangent_arcs.len(),
            plan.is_actionable(),
        );
        let face = &brep.faces[*fh];
        let outer = face.outer_loop;
        let start = brep.loops[outer].half_edge;
        let mut cur = start;
        for step in 0..16 {
            let he = &brep.half_edges[cur];
            let v_origin = brep.vertices[he.origin].position;
            let twin_face = brep.loops[brep.half_edges[he.twin].loop_].face;
            let twin_surf = brep.faces.get(twin_face).and_then(|f| f.surface.as_ref());
            let twin_kind = twin_surf
                .map(|s| match &s.geometry {
                    SurfaceGeometry::Sphere { .. } => "Sphere",
                    SurfaceGeometry::Cylinder { .. } => "Cyl",
                    SurfaceGeometry::Plane { .. } => "Plane",
                    SurfaceGeometry::Gregory(_) => "Gregory",
                    _ => "?",
                })
                .unwrap_or("none");
            let on_sphere_signed = if let Some((c, r)) = plan.cap_sphere {
                v_origin.distance_to(&c) - r
            } else {
                f64::NAN
            };
            // Also: distance from edge midpoint to cap sphere.
            let edge = &brep.edges[he.edge];
            let edge_mid_dist =
                if let (Some(coe), Some((c, r))) = (edge.curve.as_ref(), plan.cap_sphere) {
                    use popkern_core::geometry::CurveEval;
                    let t_mid = 0.5 * (coe.t_start + coe.t_end);
                    let p_mid = coe.geometry.point_at(t_mid);
                    p_mid.distance_to(&c) - r
                } else {
                    f64::NAN
                };
            eprintln!(
                "    step{step}: V({:.4},{:.4},{:.4})  twin={twin_kind:7} (face {twin_face:?})  v_to_sphere={on_sphere_signed:+.3e}  edge_mid_to_sphere={edge_mid_dist:+.3e}",
                v_origin.x, v_origin.y, v_origin.z
            );
            cur = he.next;
            if cur == start {
                break;
            }
        }
    }
}
