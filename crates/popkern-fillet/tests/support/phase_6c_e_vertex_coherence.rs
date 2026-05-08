use super::*;

/// Phase 6c-X.1: apex vertex coherence diagnostic.  For each
/// apex vertex (vertex on a cap face outer loop), enumerate ALL
/// incident edges and report the curve endpoint each one demands
/// at this vertex.  Cluster by tolerance.  If all targets agree
/// (consensus), report the consensus and current vertex offset.
/// If they disagree, report the cluster count and each cluster's
/// representative.
///
/// The output drives the X.2 vs X.4 decision:
/// * **All apex vertices have consensus** → Stage X.2 alone
///   suffices: walk cap loops, mutate vertex.position to the
///   consensus.  ~1 day.
/// * **Some apex vertices have disagreement** → Stage X.4 (mev
///   to split shared vertices) is required.  ~3-5 days.
#[ignore]
#[test]
fn phase_6c_x1_apex_vertex_coherence() {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::topology::VertexHandle;

    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");

    let cases: &[(&str, &[f64], bool, f64)] = &[
        ("cube_uniform", &[0.10, 0.10, 0.10], false, 0.0),
        ("cube_mixed", &[0.08, 0.10, 0.12], false, 0.0),
        ("cube_extreme", &[0.02, 0.10, 0.20], false, 0.0),
        ("pyramid4_mixed", &[0.04, 0.06, 0.05, 0.07], true, 0.8),
        ("pyramid4_uniform", &[0.05, 0.05, 0.05, 0.05], true, 0.8),
    ];

    let cluster_tol = 1e-6_f64;

    eprintln!("\n=== Phase 6c-X.1 apex vertex coherence ===");
    eprintln!("(cluster tolerance: {:.1e} m)", cluster_tol);
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

        // Collect apex vertices = those on any cap face's outer loop.
        let mut apex_vertices: std::collections::BTreeSet<VertexHandle> =
            std::collections::BTreeSet::new();
        for &cap in faces.iter() {
            let outer = brep.faces[cap].outer_loop;
            let start = brep.loops[outer].half_edge;
            let mut cur = start;
            for _ in 0..32 {
                let he = &brep.half_edges[cur];
                apex_vertices.insert(he.origin);
                cur = he.next;
                if cur == start {
                    break;
                }
            }
        }

        eprintln!("\n--- {name}: {} apex vertices ---", apex_vertices.len());

        let mut summary_consensus = 0usize;
        let mut summary_disagree = 0usize;
        for &v in &apex_vertices {
            let v_pos = brep.vertices[v].position;
            // Find ALL incident edges across the BRep (any HE with origin == v).
            let mut targets: Vec<(Point3, &'static str, popkern_core::topology::EdgeHandle)> =
                Vec::new();
            for (he_h, he) in brep.half_edges.iter() {
                if he.origin != v {
                    continue;
                }
                let edge_h = he.edge;
                let edge = &brep.edges[edge_h];
                let Some(coe) = edge.curve.as_ref() else {
                    continue;
                };
                // Decide whether this HE represents the curve in
                // its native direction (origin == curve(t_start)) or
                // reversed (origin == curve(t_end)) by comparing
                // distances at the (potentially stale) vertex pos.
                let p_t_start = coe.geometry.point_at(coe.t_start);
                let p_t_end = coe.geometry.point_at(coe.t_end);
                let d_to_start = v_pos.distance_to(&p_t_start);
                let d_to_end = v_pos.distance_to(&p_t_end);
                let target = if d_to_start <= d_to_end {
                    p_t_start
                } else {
                    p_t_end
                };
                let twin_face = brep.loops[brep.half_edges[he.twin].loop_].face;
                let face_kind = brep
                    .faces
                    .get(twin_face)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| match &s.geometry {
                        popkern_core::geometry::surface::SurfaceGeometry::Plane { .. } => "Plane",
                        popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. } => "Cyl",
                        popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. } => "Sphere",
                        popkern_core::geometry::surface::SurfaceGeometry::Gregory(_) => "Gregory",
                        popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_) => "Nurbs",
                        _ => "?",
                    })
                    .unwrap_or("?");
                let _ = he_h; // Unused, but kept for symmetry.
                targets.push((target, face_kind, edge_h));
            }

            // Cluster targets by tolerance.
            let mut clusters: Vec<Vec<(Point3, &'static str, popkern_core::topology::EdgeHandle)>> =
                Vec::new();
            for t in targets {
                let mut found = false;
                for c in clusters.iter_mut() {
                    if t.0.distance_to(&c[0].0) < cluster_tol {
                        c.push(t);
                        found = true;
                        break;
                    }
                }
                if !found {
                    clusters.push(vec![t]);
                }
            }

            if clusters.len() == 1 {
                summary_consensus += 1;
                let consensus = clusters[0][0].0;
                let move_dist = v_pos.distance_to(&consensus);
                eprintln!(
                    "  {v:?} pos=({:.4},{:.4},{:.4}) → CONSENSUS dist={:.3e} from {} edges",
                    v_pos.x,
                    v_pos.y,
                    v_pos.z,
                    move_dist,
                    clusters[0].len()
                );
            } else if !clusters.is_empty() {
                summary_disagree += 1;
                eprintln!(
                    "  {v:?} pos=({:.4},{:.4},{:.4}) → DISAGREEMENT {} clusters",
                    v_pos.x,
                    v_pos.y,
                    v_pos.z,
                    clusters.len()
                );
                for (i, c) in clusters.iter().enumerate() {
                    let p = c[0].0;
                    let kinds: Vec<&str> = c.iter().map(|t| t.1).collect();
                    eprintln!(
                        "    cluster[{i}] ({:.4},{:.4},{:.4}) {} edges via {:?}",
                        p.x,
                        p.y,
                        p.z,
                        c.len(),
                        kinds
                    );
                }
            } else {
                eprintln!(
                    "  {v:?} pos=({:.4},{:.4},{:.4}) → NO INCIDENT CURVES (inner spoke only)",
                    v_pos.x, v_pos.y, v_pos.z
                );
            }
        }
        eprintln!(
            "  summary: {} consensus, {} disagreement, {} no-curve",
            summary_consensus,
            summary_disagree,
            apex_vertices.len() - summary_consensus - summary_disagree,
        );
    }
}
