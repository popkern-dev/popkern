use super::*;

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
            SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline",
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
