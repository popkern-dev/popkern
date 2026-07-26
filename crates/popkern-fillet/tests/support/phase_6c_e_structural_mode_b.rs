use super::*;

/// Phase 6c-X structural dump: for each fixture, dump the cap-loop
/// topology that the rebuild needs to consume.  For each cap face,
/// reports vertex count, vertex positions, vertex-on-cap-curve
/// classification, edge curves, and twin-face surface kinds.
#[ignore]
#[test]
fn phase_6c_x_structural_dump() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceGeometry;
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");

    let cases: &[(&str, &[f64], bool, f64)] = &[
        ("cube_uniform", &[0.10, 0.10, 0.10], false, 0.0),
        ("cube_mixed", &[0.08, 0.10, 0.12], false, 0.0),
        ("cube_extreme", &[0.02, 0.10, 0.20], false, 0.0),
        ("pyramid4_mixed", &[0.04, 0.06, 0.05, 0.07], true, 0.8),
        ("pyramid4_uniform", &[0.05, 0.05, 0.05, 0.05], true, 0.8),
    ];
    let kind = |g: &SurfaceGeometry| -> &'static str {
        match g {
            SurfaceGeometry::Plane { .. } => "Plane",
            SurfaceGeometry::Cylinder { .. } => "Cyl",
            SurfaceGeometry::Sphere { .. } => "Sphere",
            SurfaceGeometry::Torus { .. } => "Torus",
            SurfaceGeometry::Cone { .. } => "Cone",
            SurfaceGeometry::Gregory(_) => "Gregory",
            SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline",
            SurfaceGeometry::Nurbs(_) => "Nurbs",
        }
    };
    eprintln!("\n=== Phase 6c-X structural dump ===");
    for (name, radii, is_pyramid, h) in cases {
        let (brep, faces) = if *is_pyramid {
            let pyr = BRep::n_gonal_pyramid(radii.len(), 1.0, *h, Tolerance::default()).unwrap();
            let vertex = pick_vertex_with_valence(&pyr, radii.len()).unwrap();
            let edges = incident_edges(&pyr, vertex);
            fillet_corner(&pyr, &edges, radii, true)
        } else {
            let cube = BRep::unit_cube(Tolerance::default()).unwrap();
            let vertex = cube.vertices.keys().next().unwrap();
            let edges = incident_edges(&cube, vertex);
            fillet_corner(&cube, &edges, radii, true)
        };
        eprintln!("\n--- {name} ---");
        eprintln!(
            "  {} cap faces, {} verts/edges/faces total",
            faces.len(),
            brep.vertices.len()
        );
        for (k, &cap) in faces.iter().enumerate() {
            let face_data = &brep.faces[cap];
            let surface_kind = face_data
                .surface
                .as_ref()
                .map(|s| kind(&s.geometry))
                .unwrap_or("none");
            eprintln!("  cap[{k}] face={cap:?} surface={surface_kind}");
            let outer = face_data.outer_loop;
            let start = brep.loops[outer].half_edge;
            let mut cur = start;
            let mut step = 0;
            loop {
                let he = &brep.half_edges[cur];
                let v_origin = brep.vertices[he.origin].position;
                let edge_h = he.edge;
                let edge = &brep.edges[edge_h];
                let twin_face = brep.loops[brep.half_edges[he.twin].loop_].face;
                let twin_kind = brep
                    .faces
                    .get(twin_face)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| kind(&s.geometry))
                    .unwrap_or("?");
                let curve_info = if let Some(coe) = edge.curve.as_ref() {
                    let p_start = coe.geometry.point_at(coe.t_start);
                    let p_end = coe.geometry.point_at(coe.t_end);
                    let gap_origin = v_origin.distance_to(&p_start);
                    let gap_dest = {
                        let dest_he = he.next;
                        let v_dest = brep.vertices[brep.half_edges[dest_he].origin].position;
                        v_dest.distance_to(&p_end)
                    };
                    format!(
                        "span={:.3e} v↔curve gaps: origin={:.3e} dest={:.3e}",
                        coe.t_end - coe.t_start,
                        gap_origin,
                        gap_dest
                    )
                } else {
                    "<no curve>".to_string()
                };
                eprintln!(
                    "    step{step:2} v={:?} pos=({:.4},{:.4},{:.4}) edge={edge_h:?} twin={twin_kind:6} {curve_info}",
                    he.origin, v_origin.x, v_origin.y, v_origin.z
                );
                cur = he.next;
                step += 1;
                if cur == start || step > 32 {
                    break;
                }
            }
        }
    }
}

/// Mode B caveat diagnostic: run BRep validation on the post-fix
/// fixtures.  The Mode B fix sets edge.curve.t_start/t_end to span
/// the full segment range, but the BRep edge's vertex endpoints
/// don't match `curve(t_start) / curve(t_end)`.  This test reports
/// whether `validate_watertight` (Check 10) actually flags the
/// mismatch and at what magnitude.
#[ignore]
#[test]
fn phase_6c_e_mode_b_validation_diagnostic() {
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");

    let cases: &[(&str, &[f64], bool)] = &[
        ("cube_uniform", &[0.10, 0.10, 0.10], false),
        ("cube_mixed", &[0.08, 0.10, 0.12], false),
        ("cube_extreme", &[0.02, 0.10, 0.20], false),
        ("pyramid4_mixed", &[0.04, 0.06, 0.05, 0.07], true),
        ("pyramid4_uniform", &[0.05, 0.05, 0.05, 0.05], true),
    ];
    eprintln!("\n=== Mode B validation diagnostic ===");
    eprintln!("tolerance.linear() = {:.3e}", Tolerance::default().linear());
    for (name, radii, is_pyramid) in cases {
        let (brep, _faces) = if *is_pyramid {
            let pyr = BRep::n_gonal_pyramid(radii.len(), 1.0, 0.8, Tolerance::default()).unwrap();
            let vertex = pick_vertex_with_valence(&pyr, radii.len()).unwrap();
            let edges = incident_edges(&pyr, vertex);
            fillet_corner(&pyr, &edges, radii, true)
        } else {
            let cube = BRep::unit_cube(Tolerance::default()).unwrap();
            let vertex = cube.vertices.keys().next().unwrap();
            let edges = incident_edges(&cube, vertex);
            fillet_corner(&cube, &edges, radii, true)
        };
        let errors = brep.validate_all();
        eprintln!("  {name:18}: {} validation errors", errors.len());
        for (i, err) in errors.iter().enumerate().take(8) {
            eprintln!("    [{i}] {err:?}");
        }
        if errors.len() > 8 {
            eprintln!("    ... ({} more)", errors.len() - 8);
        }
    }
}

/// Mode B root-cause diagnostic: dump the curve3d parameter range
/// plus endpoint positions for each cap-vs-Plane bridge edge on
/// `cube_mixed`.  Confirms whether `t_start == t_end` (collapsed
/// parameter range) is the actual cause of the lift failure.
#[ignore]
#[test]
fn phase_6c_e_cube_mixed_bridge_edge_param_ranges() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceGeometry;
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    // Note: APEX_FACE_LIFT *off* so we observe the pre-lift edge state.
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let radii = [0.08_f64, 0.10, 0.12];
    let (brep, faces) = fillet_corner(&cube, &edges, &radii, true);
    let cap_set: std::collections::BTreeSet<_> = faces.iter().copied().collect();

    eprintln!("\n=== cube_mixed bridge-edge param ranges (apex-face-lift OFF) ===");
    for (k, &cap) in faces.iter().enumerate() {
        let outer = brep.faces[cap].outer_loop;
        let start = brep.loops[outer].half_edge;
        let mut cur = start;
        let mut step = 0;
        loop {
            let he = &brep.half_edges[cur];
            let twin_face = brep.loops[brep.half_edges[he.twin].loop_].face;
            if !cap_set.contains(&twin_face)
                && let Some(sof) = brep.faces.get(twin_face).and_then(|f| f.surface.as_ref())
                && matches!(sof.geometry, SurfaceGeometry::Plane { .. })
            {
                let edge_h = he.edge;
                let edge = &brep.edges[edge_h];
                if let Some(coe) = edge.curve.as_ref() {
                    let p_start = coe.geometry.point_at(coe.t_start);
                    let p_end = coe.geometry.point_at(coe.t_end);
                    let chord = p_start.distance_to(&p_end);
                    eprintln!(
                        "  cap[{k}] step{step} edge={edge_h:?} t=[{:.6},{:.6}] span={:.3e} chord={chord:.3e}",
                        coe.t_start,
                        coe.t_end,
                        coe.t_end - coe.t_start,
                    );
                    eprintln!(
                        "      p_start=({:.4},{:.4},{:.4}) p_end=({:.4},{:.4},{:.4})",
                        p_start.x, p_start.y, p_start.z, p_end.x, p_end.y, p_end.z
                    );
                    let v_origin = brep.vertices[he.origin].position;
                    let dest_he = he.next;
                    let v_dest = brep.vertices[brep.half_edges[dest_he].origin].position;
                    eprintln!(
                        "      v_origin=({:.4},{:.4},{:.4}) v_dest=({:.4},{:.4},{:.4})",
                        v_origin.x, v_origin.y, v_origin.z, v_dest.x, v_dest.y, v_dest.z
                    );
                }
            }
            cur = he.next;
            step += 1;
            if cur == start || step > 32 {
                break;
            }
        }
    }
}

/// Mode B diagnostic: walk each cube_mixed cap face's outer loop
/// and print every half-edge with its twin-face surface kind and
/// twin-face handle.  Used to verify how many distinct plane
/// neighbours each cap has.
#[ignore]
#[test]
fn phase_6c_e_cube_mixed_cap_loop_walk() {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let radii = [0.08_f64, 0.10, 0.12];
    let (brep, faces) = fillet_corner(&cube, &edges, &radii, true);
    let cap_set: std::collections::BTreeSet<_> = faces.iter().copied().collect();
    let kind = |g: &SurfaceGeometry| -> &'static str {
        match g {
            SurfaceGeometry::Plane { .. } => "Plane",
            SurfaceGeometry::Cylinder { .. } => "Cyl",
            SurfaceGeometry::Sphere { .. } => "Sphere",
            SurfaceGeometry::Torus { .. } => "Torus",
            SurfaceGeometry::Cone { .. } => "Cone",
            SurfaceGeometry::Gregory(_) => "Gregory",
            SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline",
            SurfaceGeometry::Nurbs(_) => "Nurbs",
        }
    };
    eprintln!("\n=== cube_mixed cap-loop walk ===");
    for (k, &cap) in faces.iter().enumerate() {
        let face_data = &brep.faces[cap];
        let outer = face_data.outer_loop;
        let start = brep.loops[outer].half_edge;
        eprintln!("  cap[{k}] face={cap:?} outer_loop={outer:?}");
        let mut cur = start;
        let mut step = 0;
        loop {
            let he = &brep.half_edges[cur];
            let edge_h = he.edge;
            let twin_he = he.twin;
            let twin_face = brep.loops[brep.half_edges[twin_he].loop_].face;
            let twin_kind = brep
                .faces
                .get(twin_face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| kind(&s.geometry))
                .unwrap_or("?");
            let is_cap = cap_set.contains(&twin_face);
            eprintln!(
                "    step{step:2} he={cur:?} edge={edge_h:?} twin_face={twin_face:?} ({twin_kind}) cap_neighbor={is_cap}"
            );
            cur = he.next;
            step += 1;
            if cur == start || step > 32 {
                break;
            }
        }
    }
}
