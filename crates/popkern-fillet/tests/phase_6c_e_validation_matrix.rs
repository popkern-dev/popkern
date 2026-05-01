//! Phase 6c-E validation matrix.
//!
//! Runs the apex-fillet fixture matrix under Stage 4.7e
//! (`POPKERN_SPHERE_FIT_BOUNDARY=1` + `POPKERN_STAGE4_WIDE_WIRE_TOL=1`)
//! and captures `cap_outer_g0` / `cap_outer_g1` for each fixture.
//!
//! Goal: confirm the Stage 4.7e win on `pyramid4_mixed`
//! (4.85e-10 m worst case) generalises to the rest of the matrix
//! before flipping the flag default in Phase 6c-F.
//!
//! Commercial gate: `cap_outer_g0` ≤ 5e-3 m.
//!
//! Each fixture runs as a separate test so failures localise.
//! `#[ignore]` keeps them out of the default test run while we
//! shake out which ones still need attention.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{
    CapOuterMetrics, N_SAMPLES, fillet_corner, incident_edges, measure_cap_outer_boundary,
    pick_vertex_with_valence,
};

/// Phase 6c-E commercial-grade gate. Stage 4.7e on
/// `pyramid4_mixed` reports 4.85e-10 m; 5e-3 m is the user-
/// visible quality bar.
const PHASE_6C_E_G0_GATE: f64 = 5e-3;

/// RAII env-var guard.
struct EnvGuard {
    key: &'static str,
    prior: Option<String>,
}

impl EnvGuard {
    fn set(key: &'static str, value: &str) -> Self {
        let prior = std::env::var(key).ok();
        unsafe {
            std::env::set_var(key, value);
        }
        Self { key, prior }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        unsafe {
            match &self.prior {
                Some(v) => std::env::set_var(self.key, v),
                None => std::env::remove_var(self.key),
            }
        }
    }
}

fn run_cube_with_stage4(radii: &[f64]) -> CapOuterMetrics {
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let (brep, faces) = fillet_corner(&cube, &edges, radii, true);
    measure_cap_outer_boundary(&brep, &faces, N_SAMPLES)
}

fn run_n_pyramid_with_stage4(n: usize, h: f64, radii: &[f64]) -> CapOuterMetrics {
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let pyr = BRep::n_gonal_pyramid(n, 1.0, h, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, n).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let (brep, faces) = fillet_corner(&pyr, &edges, radii, true);
    measure_cap_outer_boundary(&brep, &faces, N_SAMPLES)
}

fn report(fixture: &str, m: &CapOuterMetrics) {
    eprintln!(
        "  {fixture:24}: cap_outer_g0={:.3e} g1={:.3e} edges={}",
        m.worst_g0_gap, m.worst_g1_coplanarity, m.edges_measured
    );
}

fn assert_g0(fixture: &str, m: &CapOuterMetrics) {
    assert!(
        m.edges_measured >= 1,
        "[{fixture}] no cap-vs-non-cap edges measured"
    );
    assert!(
        m.worst_g0_gap < PHASE_6C_E_G0_GATE,
        "[{fixture}] cap_outer_g0 {:.3e} exceeds {:.1e} commercial gate",
        m.worst_g0_gap,
        PHASE_6C_E_G0_GATE
    );
}

#[ignore]
#[test]
fn phase_6c_e_matrix_probe() {
    eprintln!(
        "\n=== Phase 6c-E validation matrix \
         (POPKERN_SPHERE_FIT_BOUNDARY=1, POPKERN_STAGE4_WIDE_WIRE_TOL=1) ==="
    );
    let cu = run_cube_with_stage4(&[0.10, 0.10, 0.10]);
    report("cube_uniform", &cu);
    let cm = run_cube_with_stage4(&[0.08, 0.10, 0.12]);
    report("cube_mixed", &cm);
    let ce = run_cube_with_stage4(&[0.02, 0.10, 0.20]);
    report("cube_extreme", &ce);
    let p3 = run_n_pyramid_with_stage4(3, 1.2, &[0.10, 0.10, 0.10]);
    report("pyramid3_uniform", &p3);
    let p4u = run_n_pyramid_with_stage4(4, 0.8, &[0.05, 0.05, 0.05, 0.05]);
    report("pyramid4_uniform", &p4u);
    let p4m = run_n_pyramid_with_stage4(4, 0.8, &[0.04, 0.06, 0.05, 0.07]);
    report("pyramid4_mixed", &p4m);
}

#[test]
fn phase_6c_e_pyramid4_mixed_under_5mm() {
    let m = run_n_pyramid_with_stage4(4, 0.8, &[0.04, 0.06, 0.05, 0.07]);
    report("pyramid4_mixed", &m);
    assert_g0("pyramid4_mixed", &m);
}

#[ignore]
#[test]
fn phase_6c_e_cube_uniform_under_5mm() {
    let m = run_cube_with_stage4(&[0.10, 0.10, 0.10]);
    report("cube_uniform", &m);
    assert_g0("cube_uniform", &m);
}

#[test]
fn phase_6c_e_cube_mixed_under_5mm() {
    let m = run_cube_with_stage4(&[0.08, 0.10, 0.12]);
    report("cube_mixed", &m);
    assert_g0("cube_mixed", &m);
}

#[test]
fn phase_6c_e_cube_extreme_under_5mm() {
    let m = run_cube_with_stage4(&[0.02, 0.10, 0.20]);
    report("cube_extreme", &m);
    assert_g0("cube_extreme", &m);
}

#[ignore]
#[test]
fn phase_6c_e_pyramid3_uniform_under_5mm() {
    let m = run_n_pyramid_with_stage4(3, 1.2, &[0.10, 0.10, 0.10]);
    report("pyramid3_uniform", &m);
    assert_g0("pyramid3_uniform", &m);
}

#[ignore]
#[test]
fn phase_6c_e_pyramid4_uniform_under_5mm() {
    let m = run_n_pyramid_with_stage4(4, 0.8, &[0.05, 0.05, 0.05, 0.05]);
    report("pyramid4_uniform", &m);
    assert_g0("pyramid4_uniform", &m);
}

/// Generic per-edge breakdown — used by per-fixture diagnostics.
fn breakdown(brep: &popkern_core::topology::BRep, faces: &[popkern_core::topology::FaceHandle]) {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::pcurve::PcurveGeometry;
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
    let cap_set: std::collections::BTreeSet<_> = faces.iter().copied().collect();
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
    let pkind = |g: &PcurveGeometry| -> &'static str {
        match g {
            PcurveGeometry::Linear { .. } => "Linear",
            PcurveGeometry::Elliptical { .. } => "Elliptical",
            PcurveGeometry::Sinusoidal { .. } => "Sinusoidal",
            PcurveGeometry::Nurbs(_) => "Nurbs",
            PcurveGeometry::CylinderArcViaRationalQuad { .. } => "CylArcRatQuad",
        }
    };
    let mut idx = 0;
    for (_eh, edge) in brep.edges.iter() {
        let Some(coe) = edge.curve.as_ref() else {
            continue;
        };
        let he_a = edge.half_edge;
        let Some(ha) = brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = brep.half_edges.get(ha.twin) else {
            continue;
        };
        let face_a = brep.loops[ha.loop_].face;
        let face_b = brep.loops[hb.loop_].face;
        if face_a == face_b {
            continue;
        }
        let a_is_cap = cap_set.contains(&face_a);
        let b_is_cap = cap_set.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        let Some(pa) = ha.pcurve.as_ref() else {
            continue;
        };
        let Some(pb) = hb.pcurve.as_ref() else {
            continue;
        };
        let Some(sa) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let Some(sb) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let cap_kind = if a_is_cap {
            kind(&sa.geometry)
        } else {
            kind(&sb.geometry)
        };
        let twin_kind = if a_is_cap {
            kind(&sb.geometry)
        } else {
            kind(&sa.geometry)
        };
        let (t_lo, t_hi) = (coe.t_start, coe.t_end);
        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
        let t_a = pa_lo + 0.5 * (pa_hi - pa_lo);
        let t_b = pb_lo + 0.5 * (pb_hi - pb_lo);
        let p_a = {
            let uv = pa.point_at(t_a);
            sa.geometry.point_at(uv.x, uv.y)
        };
        let p_b = {
            let uv = pb.point_at(t_b);
            sb.geometry.point_at(uv.x, uv.y)
        };
        let p3_mid = coe.geometry.point_at(0.5 * (t_lo + t_hi));
        let g0 = p_a.distance_to(&p_b);
        let dev_a = p3_mid.distance_to(&p_a);
        let dev_b = p3_mid.distance_to(&p_b);
        idx += 1;
        eprintln!(
            "  edge#{idx} cap={cap_kind} twin={twin_kind} cap_is_a={a_is_cap}  \
             pa={} pb={}  g0={g0:.3e}  dev_cap={dev_a:.3e}  dev_twin={dev_b:.3e}",
            pkind(&pa.geometry),
            pkind(&pb.geometry),
        );
    }
}

/// Diagnostic: per-edge breakdown for a non-passing fixture.
/// Reports cap+twin surface kind, pcurve kind on each side, and
/// G0 gap.  Used to localise why Stage 4.7e isn't activating
/// on a given fixture.
#[ignore]
#[test]
fn phase_6c_e_cube_mixed_per_edge_breakdown() {
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let radii = [0.08_f64, 0.10, 0.12];
    let (brep, faces) = fillet_corner(&cube, &edges, &radii, true);
    eprintln!("\n=== cube_mixed cap-outer per-edge breakdown ===");
    eprintln!(
        "{} cap faces; {} BRep edges total",
        faces.len(),
        brep.edges.len()
    );
    breakdown(&brep, &faces);
}

#[ignore]
#[test]
fn phase_6c_e_pyramid4_uniform_per_edge_breakdown() {
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let pyr = BRep::n_gonal_pyramid(4, 1.0, 0.8, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, 4).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let radii = [0.05_f64, 0.05, 0.05, 0.05];
    let (brep, faces) = fillet_corner(&pyr, &edges, &radii, true);
    eprintln!("\n=== pyramid4_uniform cap-outer per-edge breakdown ===");
    eprintln!(
        "{} cap faces; {} BRep edges total",
        faces.len(),
        brep.edges.len()
    );
    breakdown(&brep, &faces);
}

#[ignore]
#[test]
fn phase_6c_e_cube_uniform_per_edge_breakdown() {
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let radii = [0.10_f64, 0.10, 0.10];
    let (brep, faces) = fillet_corner(&cube, &edges, &radii, true);
    eprintln!("\n=== cube_uniform cap-outer per-edge breakdown ===");
    eprintln!(
        "{} cap faces; {} BRep edges total",
        faces.len(),
        brep.edges.len()
    );
    breakdown(&brep, &faces);
}

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

#[ignore]
#[test]
fn phase_6c_e_pyramid3_uniform_per_edge_breakdown() {
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let pyr = BRep::n_gonal_pyramid(3, 1.0, 1.2, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, 3).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let radii = [0.10_f64, 0.10, 0.10];
    let (brep, faces) = fillet_corner(&pyr, &edges, &radii, true);
    eprintln!("\n=== pyramid3_uniform cap-outer per-edge breakdown ===");
    eprintln!(
        "{} cap faces; {} BRep edges total",
        faces.len(),
        brep.edges.len()
    );
    breakdown(&brep, &faces);
}
