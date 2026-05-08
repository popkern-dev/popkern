use super::*;

#[path = "kp_binonic_boundary_injection_pcurve_dump.rs"]
mod kp_binonic_boundary_injection_pcurve_dump;

use kp_binonic_boundary_injection_pcurve_dump::dump_pcurve_trajectory;

#[ignore]
#[test]
fn cap_surface_at_row0_tracks_analytic_boundary() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::gregory::GregoryPatch as GP;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let (_guard, brep, faces) = pyramid4_mixed_injection_fixture();

    eprintln!("\n=== pyramid4_mixed cap-surface direct eval at row 0 ===");
    for (i, &fh) in faces.iter().enumerate() {
        let face = &brep.faces[fh];
        let Some(sof) = face.surface.as_ref() else {
            eprintln!("  cap[{i}]: no surface");
            continue;
        };
        let SurfaceGeometry::Gregory(GP::RectangularBiNonic(g)) = &sof.geometry else {
            eprintln!("  cap[{i}]: not RectangularBiNonic");
            continue;
        };
        for (label, samples) in [
            ("row0", [(0.0, 0.0), (0.0, 0.5), (0.0, 1.0)]),
            ("col9", [(0.0, 1.0), (0.5, 1.0), (1.0, 1.0)]),
        ] {
            let points = samples
                .into_iter()
                .map(|(u, v)| g.point_at(u, v))
                .collect::<Vec<_>>();
            eprintln!(
                "  cap[{i}] {label}: ({:.4},{:.4},{:.4}) ({:.4},{:.4},{:.4}) ({:.4},{:.4},{:.4})",
                points[0].x,
                points[0].y,
                points[0].z,
                points[1].x,
                points[1].y,
                points[1].z,
                points[2].x,
                points[2].y,
                points[2].z,
            );
        }
    }

    eprintln!("\n=== checking each cap-vs-cylinder edge ===");
    for &fh in &faces {
        let face = &brep.faces[fh];
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        let SurfaceGeometry::Gregory(GP::RectangularBiNonic(g)) = &sof.geometry else {
            continue;
        };
        let outer = face.outer_loop;
        let start = brep.loops[outer].half_edge;
        let mut cur = start;
        let mut step = 0;
        loop {
            let he = &brep.half_edges[cur];
            let twin_face = brep.loops[brep.half_edges[he.twin].loop_].face;
            let twin = brep.faces.get(twin_face).and_then(|f| f.surface.as_ref());
            let is_cyl = matches!(
                twin.map(|s| &s.geometry),
                Some(SurfaceGeometry::Cylinder { .. })
            );
            if is_cyl && let Some(coe) = brep.edges[he.edge].curve.as_ref() {
                let t_lo = coe.t_start;
                let t_hi = coe.t_end;
                let t_mid = 0.5 * (t_lo + t_hi);
                let p3 = coe.geometry.point_at(t_mid);
                let mut best = (f64::INFINITY, 0.0_f64, 0u8);
                for k in 0..=20 {
                    let s = (k as f64) / 20.0;
                    for (boundary, p) in [(0, g.point_at(0.0, s)), (1, g.point_at(s, 1.0))] {
                        let d = p3.distance_to(&p);
                        if d < best.0 {
                            best = (d, s, boundary);
                        }
                    }
                }
                let bnd = if best.2 == 0 { "row0" } else { "col9" };
                eprintln!(
                    "  cap face step{step}: edge_chord={:.3e} t_mid={t_mid:.3e} cap_to_edge_min={:.3e} on {bnd}@{:.2}",
                    coe.geometry
                        .point_at(t_lo)
                        .distance_to(&coe.geometry.point_at(t_hi)),
                    best.0,
                    best.1,
                );
            }
            cur = he.next;
            step += 1;
            if cur == start || step > 16 {
                break;
            }
        }
    }
}

#[ignore]
#[test]
fn cap_loop_vs_cap_outer_edge_coverage() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let (_guard, brep, faces) = pyramid4_mixed_injection_fixture();
    let cap_set: std::collections::BTreeSet<_> = faces.iter().copied().collect();

    eprintln!(
        "\n=== cap face outer-loop edges (from wire's perspective) ===\n{} cap faces",
        faces.len()
    );
    let mut wire_edges: std::collections::BTreeSet<popkern_core::topology::EdgeHandle> =
        std::collections::BTreeSet::new();
    for (face_idx, &fh) in faces.iter().enumerate() {
        let face = &brep.faces[fh];
        let outer = face.outer_loop;
        let start = brep.loops[outer].half_edge;
        eprintln!("\nCap face #{face_idx} ({fh:?}):");
        let mut cur = start;
        for step in 0..16 {
            let he = &brep.half_edges[cur];
            let edge = he.edge;
            wire_edges.insert(edge);
            let twin_face = brep.loops[brep.half_edges[he.twin].loop_].face;
            let twin_kind = brep
                .faces
                .get(twin_face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| match &s.geometry {
                    SurfaceGeometry::Plane { .. } => "Plane",
                    SurfaceGeometry::Cylinder { .. } => "Cyl",
                    SurfaceGeometry::Sphere { .. } => "Sphere",
                    SurfaceGeometry::Gregory(_) => "Gregory",
                    SurfaceGeometry::Nurbs(_) => "Nurbs",
                    _ => "Other",
                })
                .unwrap_or("none");
            let v_origin = brep.vertices[he.origin].position;
            let curve_mid = brep.edges[edge].curve.as_ref().map(|coe| {
                let t_mid = 0.5 * (coe.t_start + coe.t_end);
                coe.geometry.point_at(t_mid)
            });
            eprintln!(
                "  step{step}: he={cur:?} edge={edge:?} twin={twin_kind:7} (face {twin_face:?}) v_origin=({:.3},{:.3},{:.3}) curve_mid={:?}",
                v_origin.x,
                v_origin.y,
                v_origin.z,
                curve_mid.map(|p| (p.x, p.y, p.z))
            );
            cur = he.next;
            if cur == start {
                break;
            }
        }
    }

    let cap_outer_edges = cap_outer_edge_set(&brep, &cap_set);
    eprintln!("cap_outer measures {} unique edges", cap_outer_edges.len());
    eprintln!(
        "wire-only edges (not measured by cap_outer): {:?}",
        wire_edges.difference(&cap_outer_edges).collect::<Vec<_>>()
    );
    eprintln!(
        "cap_outer-only edges (not visited by wire): {:?}",
        cap_outer_edges.difference(&wire_edges).collect::<Vec<_>>()
    );
}

#[ignore]
#[test]
fn cap_outer_per_edge_breakdown_pyramid4_mixed() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

    let (_guard, brep, faces) = pyramid4_mixed_injection_fixture();
    let cap_set: std::collections::BTreeSet<_> = faces.iter().copied().collect();

    eprintln!("\n=== pyramid4_mixed cap-outer per-edge breakdown ===");
    eprintln!(
        "{} cap faces; {} BRep edges total",
        faces.len(),
        brep.edges.len()
    );

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
        let kind_a = kind(&sa.geometry);
        let kind_b = kind(&sb.geometry);
        let cap_kind = if a_is_cap { kind_a } else { kind_b };
        let twin_kind = if a_is_cap { kind_b } else { kind_a };

        let (t_lo, t_hi) = (coe.t_start, coe.t_end);
        let t_mid = 0.5 * (t_lo + t_hi);
        let p3_mid = coe.geometry.point_at(t_mid);
        let p3_lo = coe.geometry.point_at(t_lo);
        let p3_hi = coe.geometry.point_at(t_hi);
        let edge_len_chord = p3_lo.distance_to(&p3_hi);

        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
        let t_a = pa_lo + 0.5 * (pa_hi - pa_lo);
        let t_b = pb_lo + 0.5 * (pb_hi - pb_lo);
        let uv_a = pa.point_at(t_a);
        let uv_b = pb.point_at(t_b);
        let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
        let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
        let g0 = p_a.distance_to(&p_b);
        let dev_a = p3_mid.distance_to(&p_a);
        let dev_b = p3_mid.distance_to(&p_b);

        idx += 1;
        let (cap_uv, twin_uv) = if a_is_cap { (uv_a, uv_b) } else { (uv_b, uv_a) };
        eprintln!(
            "  edge#{idx} cap={cap_kind} twin={twin_kind} cap_is_a={a_is_cap} chord={edge_len_chord:.3e} g0={g0:.3e} dev_cap={dev_a:.3e} dev_twin={dev_b:.3e}",
        );
        eprintln!(
            "         cap_uv=({:.3},{:.3}) twin_uv=({:.3},{:.3}) ranges pa=[{pa_lo:.3},{pa_hi:.3}] pb=[{pb_lo:.3},{pb_hi:.3}] coe=[{t_lo:.3},{t_hi:.3}]",
            cap_uv.x, cap_uv.y, twin_uv.x, twin_uv.y
        );
        eprintln!(
            "         curve3d: lo=({:.3},{:.3},{:.3}) mid=({:.3},{:.3},{:.3}) hi=({:.3},{:.3},{:.3})",
            p3_lo.x, p3_lo.y, p3_lo.z, p3_mid.x, p3_mid.y, p3_mid.z, p3_hi.x, p3_hi.y, p3_hi.z,
        );
        eprintln!(
            "         cap_eval=({:.3},{:.3},{:.3}) twin_eval=({:.3},{:.3},{:.3})",
            p_a.x, p_a.y, p_a.z, p_b.x, p_b.y, p_b.z,
        );
    }
}

#[test]
#[ignore = "diagnostic - run on demand"]
fn cap_surface_boundary_and_pcurve_trajectory_pyramid4_mixed() {
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

    let (_guard, brep, faces) = pyramid4_mixed_injection_fixture();
    eprintln!("\n=== pyramid4_mixed cap surface boundary + pcurve trajectory ===");
    for (cap_k, cap_face) in faces.iter().enumerate() {
        let Some(face_data) = brep.faces.get(*cap_face) else {
            continue;
        };
        let Some(sof) = face_data.surface.as_ref() else {
            continue;
        };
        let SurfaceGeometry::Gregory(_) = &sof.geometry else {
            continue;
        };

        eprintln!("\n  cap[{cap_k}] {cap_face:?}:");
        for (name, samples) in [
            (
                "row 0",
                [(0.0, 0.0), (0.0, 0.25), (0.0, 0.5), (0.0, 0.75), (0.0, 1.0)],
            ),
            (
                "col 9",
                [(0.0, 1.0), (0.25, 1.0), (0.5, 1.0), (0.75, 1.0), (1.0, 1.0)],
            ),
            (
                "row 9",
                [(1.0, 0.0), (1.0, 0.25), (1.0, 0.5), (1.0, 0.75), (1.0, 1.0)],
            ),
            (
                "col 0",
                [(0.0, 0.0), (0.25, 0.0), (0.5, 0.0), (0.75, 0.0), (1.0, 0.0)],
            ),
        ] {
            eprintln!("    {name}:");
            for (u, v) in samples {
                let p = sof.geometry.point_at(u, v);
                eprintln!(
                    "      ({u:.2},{v:.2}): ({:.4}, {:.4}, {:.4})",
                    p.x, p.y, p.z
                );
            }
        }

        let cap_set: std::collections::BTreeSet<_> = faces.iter().copied().collect();
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let start_he = loop_data.half_edge;
        let mut cur = start_he;
        let mut step = 0;
        loop {
            let Some(he) = brep.half_edges.get(cur) else {
                break;
            };
            let twin_face = brep.loops[brep.half_edges[he.twin].loop_].face;
            if !cap_set.contains(&twin_face)
                && let Some(pa) = he.pcurve.as_ref()
            {
                dump_pcurve_trajectory(&brep, sof, he, twin_face, step, pa);
            }
            cur = he.next;
            step += 1;
            if cur == start_he || step > 32 {
                break;
            }
        }
    }
}

#[test]
#[ignore = "diagnostic - run on demand"]
fn faceedge_bridge_offplane_distance_pyramid4_mixed() {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let (_guard, brep, faces) = pyramid4_mixed_injection_fixture();
    let cap_set: std::collections::BTreeSet<_> = faces.iter().copied().collect();

    eprintln!("\n=== pyramid4_mixed FaceEdge bridge off-plane distance ===");
    let mut idx = 0;
    let mut max_off = 0.0f64;
    for (eh, edge) in brep.edges.iter() {
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
        let Some(coe) = edge.curve.as_ref() else {
            continue;
        };
        let twin_face = if a_is_cap { face_b } else { face_a };
        let Some(twin_sof) = brep.faces.get(twin_face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let SurfaceGeometry::Plane { origin, normal, .. } = &twin_sof.geometry else {
            continue;
        };
        let Some(n_unit) = normal.normalized() else {
            continue;
        };
        let plane_dist = |p: popkern_core::geometry::Point3| -> f64 {
            let d = p - *origin;
            d.dot(&n_unit).abs()
        };
        let p_lo = coe.geometry.point_at(coe.t_start);
        let p_mid = coe.geometry.point_at(0.5 * (coe.t_start + coe.t_end));
        let p_hi = coe.geometry.point_at(coe.t_end);
        let d_lo = plane_dist(p_lo);
        let d_mid = plane_dist(p_mid);
        let d_hi = plane_dist(p_hi);
        max_off = max_off.max(d_lo.max(d_mid).max(d_hi));
        idx += 1;
        eprintln!(
            "  bridge#{idx} edge={eh:?} twin_face={twin_face:?} d_lo={d_lo:.3e} d_mid={d_mid:.3e} d_hi={d_hi:.3e}",
        );
    }
    eprintln!("\n  max bridge-to-plane orthogonal distance: {max_off:.3e}");
}

fn pyramid4_mixed_injection_fixture() -> (EnvGuard, BRep, Vec<popkern_core::topology::FaceHandle>) {
    let guard = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let pyr = BRep::square_pyramid(1.0, 0.8, Tolerance::default()).unwrap();
    let vertex = pick_vertex_with_valence(&pyr, 4).unwrap();
    let edges = incident_edges(&pyr, vertex);
    let radii: Vec<f64> = vec![0.04, 0.06, 0.05, 0.07];
    let (brep, faces) = fillet_corner(&pyr, &edges, &radii, true);
    (guard, brep, faces)
}

fn cap_outer_edge_set(
    brep: &BRep,
    cap_set: &std::collections::BTreeSet<popkern_core::topology::FaceHandle>,
) -> std::collections::BTreeSet<popkern_core::topology::EdgeHandle> {
    let mut cap_outer_edges = std::collections::BTreeSet::new();
    for (eh, edge) in brep.edges.iter() {
        let Some(_coe) = edge.curve.as_ref() else {
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
        if a_is_cap != b_is_cap && ha.pcurve.is_some() && hb.pcurve.is_some() {
            cap_outer_edges.insert(eh);
        }
    }
    cap_outer_edges
}
