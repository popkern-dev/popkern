use super::*;

/// Phase 6c-X.A.3 verification: with PATH_A on, every cap-loop edge
/// should have full-range curve (span = curve.parameter_range()
/// difference, no projection sub-ranging).  Verify on
/// pyramid4_mixed where cap-vs-Plane edges previously had span=0.6.
#[ignore]
#[test]
fn phase_6c_a3_full_range_check_pyramid4_mixed() {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let pyr = BRep::n_gonal_pyramid(4, 1.0, 0.8, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, 4).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let radii = [0.04_f64, 0.06, 0.05, 0.07];
    let (brep, faces) = fillet_corner(&pyr, &edges, &radii, true);
    let cap_set: std::collections::BTreeSet<_> = faces.iter().copied().collect();
    eprintln!("\n=== Phase 6c-X.A.3 full-range check pyramid4_mixed ===");
    let mut total = 0;
    let mut full_range = 0;
    for (eh, edge) in brep.edges.iter() {
        let Some(coe) = edge.curve.as_ref() else {
            continue;
        };
        let he_a = edge.half_edge;
        let he_b = brep.half_edges[he_a].twin;
        let face_a = brep.loops[brep.half_edges[he_a].loop_].face;
        let face_b = brep.loops[brep.half_edges[he_b].loop_].face;
        let a_is_cap = cap_set.contains(&face_a);
        let b_is_cap = cap_set.contains(&face_b);
        if !a_is_cap && !b_is_cap {
            continue;
        }
        total += 1;
        let span = coe.t_end - coe.t_start;
        let kind = brep
            .faces
            .get(if a_is_cap { face_b } else { face_a })
            .and_then(|f| f.surface.as_ref())
            .map(|s| match &s.geometry {
                SurfaceGeometry::Plane { .. } => "Plane",
                SurfaceGeometry::Cylinder { .. } => "Cyl",
                SurfaceGeometry::Nurbs(_) => "Nurbs",
                _ => "?",
            })
            .unwrap_or("?");
        let is_full = (span - 1.0).abs() < 1e-9;
        if is_full {
            full_range += 1;
        }
        eprintln!("  {eh:?} cap-vs-{kind}  span={span:.6} full_range={is_full}");
    }
    eprintln!("  {full_range}/{total} cap-loop edges at full segment range");
}

/// Phase 6c-X.A.2: candidate T discovery + wire-match verification.
///
/// Goal: determine whether apex-vertex disagreement is caused by
/// (a) wire matching the same vertex's incident edges to mismatched
/// boundary segments, or (b) genuine geometric disagreement.
///
/// Strategy: collect every cap-loop edge's curve endpoints across
/// the BRep.  Cluster all those 3D points by 1 mm tolerance.  Each
/// cluster represents one logical "T" (FilletEnd termination).
/// For each apex vertex V, the nearest cluster is V's "logical T".
/// Report whether V's cap-loop edges' endpoints agree with this T
/// or with DIFFERENT clusters (wire mismatch).
#[ignore]
#[test]
fn phase_6c_a2_t_discovery_and_wire_match() {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::topology::VertexHandle;

    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");

    let cases: &[(&str, &[f64], bool, f64)] = &[
        ("cube_mixed", &[0.08, 0.10, 0.12], false, 0.0),
        ("pyramid4_mixed", &[0.04, 0.06, 0.05, 0.07], true, 0.8),
    ];
    let cluster_tol = 1e-3_f64; // 1 mm

    eprintln!("\n=== Phase 6c-X.A.2 T discovery + wire match ===");
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
        let cap_set: std::collections::BTreeSet<_> = faces.iter().copied().collect();

        // Apex vertices.
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

        // Collect all cap-loop edge curve endpoints, cluster.
        let mut all_endpoints: Vec<Point3> = Vec::new();
        for (_eh, edge) in brep.edges.iter() {
            let Some(coe) = edge.curve.as_ref() else {
                continue;
            };
            let he_a = edge.half_edge;
            let he_b = brep.half_edges[he_a].twin;
            let face_a = brep.loops[brep.half_edges[he_a].loop_].face;
            let face_b = brep.loops[brep.half_edges[he_b].loop_].face;
            let a_is_cap = cap_set.contains(&face_a);
            let b_is_cap = cap_set.contains(&face_b);
            if !a_is_cap && !b_is_cap {
                continue;
            }
            all_endpoints.push(coe.geometry.point_at(coe.t_start));
            all_endpoints.push(coe.geometry.point_at(coe.t_end));
        }
        // Cluster.
        let mut clusters: Vec<Vec<Point3>> = Vec::new();
        for p in &all_endpoints {
            let mut found = false;
            for c in clusters.iter_mut() {
                if p.distance_to(&c[0]) < cluster_tol {
                    c.push(*p);
                    found = true;
                    break;
                }
            }
            if !found {
                clusters.push(vec![*p]);
            }
        }
        let cluster_centroids: Vec<Point3> = clusters
            .iter()
            .map(|c| {
                let n = c.len() as f64;
                let mut sum = Point3::new(0.0, 0.0, 0.0);
                for p in c {
                    sum.x += p.x;
                    sum.y += p.y;
                    sum.z += p.z;
                }
                Point3::new(sum.x / n, sum.y / n, sum.z / n)
            })
            .collect();
        eprintln!(
            "\n--- {name}: {} apex vertices, {} cap-loop endpoint clusters ---",
            apex_vertices.len(),
            clusters.len()
        );
        for (ci, (cs, c)) in cluster_centroids.iter().zip(clusters.iter()).enumerate() {
            eprintln!(
                "  T[{ci}] centroid=({:.4},{:.4},{:.4}) {} samples",
                cs.x,
                cs.y,
                cs.z,
                c.len()
            );
        }

        // For each apex vertex, find nearest cluster (its "logical T").
        // Then check whether all incident cap-loop edges' nearest endpoints
        // map to the SAME cluster.
        for &v in &apex_vertices {
            let v_pos = brep.vertices[v].position;
            // Find nearest cluster to v.
            let (best_cluster, dist_to_t) = cluster_centroids
                .iter()
                .enumerate()
                .map(|(i, c)| (i, v_pos.distance_to(c)))
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                .unwrap_or((usize::MAX, f64::INFINITY));
            // Walk incident cap-loop edges' nearest endpoints, map each to its cluster.
            let mut clusters_seen: std::collections::BTreeSet<usize> =
                std::collections::BTreeSet::new();
            let mut edges_with_curve = 0;
            for (_he_h, he) in brep.half_edges.iter() {
                if he.origin != v {
                    continue;
                }
                let edge = &brep.edges[he.edge];
                let Some(coe) = edge.curve.as_ref() else {
                    continue;
                };
                let face_a = brep.loops[brep.half_edges[edge.half_edge].loop_].face;
                let face_b =
                    brep.loops[brep.half_edges[brep.half_edges[edge.half_edge].twin].loop_].face;
                let a_is_cap = cap_set.contains(&face_a);
                let b_is_cap = cap_set.contains(&face_b);
                if !a_is_cap && !b_is_cap {
                    continue;
                }
                edges_with_curve += 1;
                let p_start = coe.geometry.point_at(coe.t_start);
                let p_end = coe.geometry.point_at(coe.t_end);
                let near = if v_pos.distance_to(&p_start) <= v_pos.distance_to(&p_end) {
                    p_start
                } else {
                    p_end
                };
                let nearest_cluster = cluster_centroids
                    .iter()
                    .enumerate()
                    .map(|(i, c)| (i, near.distance_to(c)))
                    .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                    .map(|x| x.0)
                    .unwrap_or(usize::MAX);
                clusters_seen.insert(nearest_cluster);
            }
            let agreement = if clusters_seen.len() <= 1 {
                "AGREE"
            } else {
                "DISAGREE"
            };
            eprintln!(
                "  {v:?} pos=({:.4},{:.4},{:.4}) nearest=T[{best_cluster}] dist={dist_to_t:.3e} \
                 {agreement}: cap-loop edges={edges_with_curve}, clusters_hit={:?}",
                v_pos.x, v_pos.y, v_pos.z, clusters_seen
            );
        }
    }
}

/// Phase 6c-X.A.1: apex region detection diagnostic.
///
/// For each apex vertex on every fixture, walk ALL incident
/// half-edges and classify each by edge role:
/// * **Cap-loop edge** (one HE on a cap face): cap-vs-Cyl or
///   cap-vs-Plane.  Has substituted curve3d ending at FilletEnd
///   termination T.
/// * **Apex generator** (BOTH HEs on non-cap faces, both faces
///   share the apex region — neither is the cap): cyl-vs-plane
///   apex generator.  Original (pre-substitution) curve.
/// * **Base edge** (incident to apex vertex but going AWAY from
///   the apex region — toward the cube base or another non-apex
///   region): unchanged geometry, no Path A action needed.
///
/// Distinguishes "apex generator" from "base edge" by checking
/// whether the OTHER endpoint of the edge is also an apex
/// vertex (both endpoints in the apex set ⇒ apex generator;
/// only one endpoint ⇒ base edge).
///
/// This drives A.2 — the rebuild plan needs to know which edges
/// to insert/remove/keep.
#[ignore]
#[test]
fn phase_6c_a1_apex_region_detection() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::topology::VertexHandle;

    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");

    let cases: &[(&str, &[f64], bool, f64)] = &[
        ("cube_mixed", &[0.08, 0.10, 0.12], false, 0.0),
        ("pyramid4_mixed", &[0.04, 0.06, 0.05, 0.07], true, 0.8),
    ];
    let kind = |g: &SurfaceGeometry| -> &'static str {
        match g {
            SurfaceGeometry::Plane { .. } => "Plane",
            SurfaceGeometry::Cylinder { .. } => "Cyl",
            SurfaceGeometry::Sphere { .. } => "Sphere",
            SurfaceGeometry::Torus { .. } => "Torus",
            SurfaceGeometry::Cone { .. } => "Cone",
            SurfaceGeometry::Gregory(_) => "Gregory",
            SurfaceGeometry::Nurbs(_) => "Nurbs",
        }
    };

    eprintln!("\n=== Phase 6c-X.A.1 apex region detection ===");
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
        let cap_set: std::collections::BTreeSet<_> = faces.iter().copied().collect();

        // Apex set = vertices on any cap face's outer loop.
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

        for &v in &apex_vertices {
            let v_pos = brep.vertices[v].position;
            eprintln!("  {v:?} pos=({:.4},{:.4},{:.4})", v_pos.x, v_pos.y, v_pos.z);
            // Find every edge incident to v (HE.origin == v OR HE.dest == v).
            let mut edges_seen: std::collections::BTreeSet<popkern_core::topology::EdgeHandle> =
                std::collections::BTreeSet::new();
            for (_he_h, he) in brep.half_edges.iter() {
                if he.origin != v {
                    continue;
                }
                edges_seen.insert(he.edge);
            }
            for edge_h in edges_seen {
                let edge = &brep.edges[edge_h];
                let he_a = edge.half_edge;
                let he_b = brep.half_edges[he_a].twin;
                let face_a = brep.loops[brep.half_edges[he_a].loop_].face;
                let face_b = brep.loops[brep.half_edges[he_b].loop_].face;
                let a_is_cap = cap_set.contains(&face_a);
                let b_is_cap = cap_set.contains(&face_b);
                let kind_a = brep
                    .faces
                    .get(face_a)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| kind(&s.geometry))
                    .unwrap_or("?");
                let kind_b = brep
                    .faces
                    .get(face_b)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| kind(&s.geometry))
                    .unwrap_or("?");
                // Find the other endpoint of this edge.
                let other_v = if brep.half_edges[he_a].origin == v {
                    brep.half_edges[brep.half_edges[he_a].next].origin
                } else {
                    brep.half_edges[he_a].origin
                };
                let other_in_apex = apex_vertices.contains(&other_v);
                let role = if a_is_cap || b_is_cap {
                    "CAP-LOOP "
                } else if other_in_apex {
                    "APEX-GEN "
                } else {
                    "BASE-EDGE"
                };
                let curve_status = if let Some(coe) = edge.curve.as_ref() {
                    let p_start = coe.geometry.point_at(coe.t_start);
                    let p_end = coe.geometry.point_at(coe.t_end);
                    let near = if v_pos.distance_to(&p_start) <= v_pos.distance_to(&p_end) {
                        p_start
                    } else {
                        p_end
                    };
                    let gap = v_pos.distance_to(&near);
                    format!(
                        "curve span={:.3e} v→curve_near={:.3e} near=({:.4},{:.4},{:.4})",
                        coe.t_end - coe.t_start,
                        gap,
                        near.x,
                        near.y,
                        near.z
                    )
                } else {
                    "<no curve>".to_string()
                };
                eprintln!(
                    "    {role} {edge_h:?} {kind_a:6}↔{kind_b:6} other={other_v:?} apex={other_in_apex}  {curve_status}"
                );
            }
        }
    }
}
