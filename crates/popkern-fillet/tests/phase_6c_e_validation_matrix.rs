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
    let _g2 = EnvGuard::set("POPKERN_STAGE4_WIDE_WIRE_TOL", "1");
    let _g3 = EnvGuard::set("POPKERN_APEX_FACE_LIFT", "1");
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let edges = incident_edges(&cube, vertex);
    let (brep, faces) = fillet_corner(&cube, &edges, radii, true);
    measure_cap_outer_boundary(&brep, &faces, N_SAMPLES)
}

fn run_n_pyramid_with_stage4(n: usize, h: f64, radii: &[f64]) -> CapOuterMetrics {
    let _g1 = EnvGuard::set("POPKERN_SPHERE_FIT_BOUNDARY", "1");
    let _g2 = EnvGuard::set("POPKERN_STAGE4_WIDE_WIRE_TOL", "1");
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
    let _g2 = EnvGuard::set("POPKERN_STAGE4_WIDE_WIRE_TOL", "1");
    let _g3 = EnvGuard::set("POPKERN_APEX_FACE_LIFT", "1");
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
    let _g2 = EnvGuard::set("POPKERN_STAGE4_WIDE_WIRE_TOL", "1");
    let _g3 = EnvGuard::set("POPKERN_APEX_FACE_LIFT", "1");
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
    let _g2 = EnvGuard::set("POPKERN_STAGE4_WIDE_WIRE_TOL", "1");
    let _g3 = EnvGuard::set("POPKERN_APEX_FACE_LIFT", "1");
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
    let _g2 = EnvGuard::set("POPKERN_STAGE4_WIDE_WIRE_TOL", "1");
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
    let _g2 = EnvGuard::set("POPKERN_STAGE4_WIDE_WIRE_TOL", "1");
    let _g3 = EnvGuard::set("POPKERN_APEX_FACE_LIFT", "1");
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
    let _g2 = EnvGuard::set("POPKERN_STAGE4_WIDE_WIRE_TOL", "1");
    let _g3 = EnvGuard::set("POPKERN_APEX_FACE_LIFT", "1");
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
