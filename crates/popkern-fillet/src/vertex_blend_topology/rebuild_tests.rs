//! End-to-end Phase 6c-C orchestration tests.
//!
//! `cube_uniform_apex_rebuild_e2e` exercises the full
//! `fillet_edges` pipeline with `POPKERN_APEX_TOPOLOGY_REBUILD=1`,
//! verifies the cap face's outer loop drops from 6 edges to 3
//! edges (the canonical sphere-fit triangle), and measures the
//! cap-vs-cylinder G0 / G1 deviation post-rebuild.
//!
//! Counterpart `cube_uniform_apex_no_rebuild_baseline` measures
//! the same metrics WITHOUT the env-var, so the deltas can be
//! reported as a regression-guard.
//!
//! Marked `#[ignore]` because they set process-level env vars;
//! run with `--ignored`.

use crate::FilletOptions;
use crate::types::{FilletRadius, FilletSpec};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

fn cube_apex(brep: &BRep) -> VertexHandle {
    brep.vertices
        .keys()
        .find(|&v| brep.half_edges.values().filter(|he| he.origin == v).count() == 3)
        .expect("cube has at least one valence-3 vertex")
}

fn corner_edges(brep: &BRep, vertex: VertexHandle) -> Vec<popkern_core::topology::EdgeHandle> {
    use std::collections::BTreeSet;
    brep.half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn fillet_uniform_with_rebuild(rebuild: bool) -> (BRep, FaceHandle) {
    unsafe {
        if rebuild {
            std::env::set_var("POPKERN_APEX_TOPOLOGY_REBUILD", "1");
        } else {
            std::env::remove_var("POPKERN_APEX_TOPOLOGY_REBUILD");
        }
    }
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let apex_pre = cube_apex(&cube);
    let edges = corner_edges(&cube, apex_pre);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.10),
        })
        .collect();
    let result = crate::fillet_edges(
        &cube,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("fillet uniform cube");
    let cap_face = result.corner_patches[0].face;
    (result.brep, cap_face)
}

fn cap_loop_edge_count(brep: &BRep, cap_face: FaceHandle) -> usize {
    let outer_loop = brep.faces[cap_face].outer_loop;
    let start = brep.loops[outer_loop].half_edge;
    let mut count = 0;
    let mut cursor = start;
    loop {
        count += 1;
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            break;
        }
        if count > 64 {
            break;
        }
    }
    count
}

/// Measure the worst cap-vs-cylinder G0 / G1 along the cap face's
/// outer loop. Samples `n_samples=7` per edge. Determines pcurve
/// direction by comparing endpoint distances to vertex positions.
fn cap_vs_cylinder_g0_g1(brep: &BRep, cap_face: FaceHandle) -> (f64, f64) {
    use popkern_core::geometry::pcurve::Pcurve;
    let outer_loop = brep.faces[cap_face].outer_loop;
    let start = brep.loops[outer_loop].half_edge;
    let cap_surf = brep.faces[cap_face]
        .surface
        .as_ref()
        .map(|s| s.geometry.clone())
        .expect("cap face has surface");

    let pcurve_forward = |pc: &Pcurve,
                          surf: &SurfaceGeometry,
                          p_origin: popkern_core::geometry::Point3,
                          p_dest: popkern_core::geometry::Point3|
     -> bool {
        let (lo, hi) = pc.parameter_range();
        let uv_lo = pc.point_at(lo);
        let uv_hi = pc.point_at(hi);
        let p_lo = surf.point_at(uv_lo.x, uv_lo.y);
        let p_hi = surf.point_at(uv_hi.x, uv_hi.y);
        p_origin.distance_to(&p_lo) + p_dest.distance_to(&p_hi)
            <= p_origin.distance_to(&p_hi) + p_dest.distance_to(&p_lo)
    };

    let n_samples = 7usize;
    let mut max_g0 = 0.0_f64;
    let mut max_g1 = 0.0_f64;
    let mut cursor = start;
    let mut safety = 64;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        let he_data = brep.half_edges[cursor].clone();
        let twin_data = brep.half_edges[he_data.twin].clone();
        let twin_face = brep.loops[twin_data.loop_].face;
        let other_surf = brep.faces[twin_face]
            .surface
            .as_ref()
            .map(|s| s.geometry.clone());
        let pa = he_data.pcurve.clone();
        let pb = twin_data.pcurve.clone();
        let v_origin = brep.vertices[he_data.origin].position;
        let v_dest = brep.vertices[twin_data.origin].position;
        if let (Some(other), Some(pa), Some(pb)) = (other_surf, pa, pb) {
            // Only measure cap-vs-cylinder edges.
            if matches!(other, SurfaceGeometry::Cylinder { .. }) {
                let pa_fwd = pcurve_forward(&pa, &cap_surf, v_origin, v_dest);
                // The twin's half-edge runs from v_dest to v_origin,
                // so swap the endpoints when checking direction.
                let pb_fwd = pcurve_forward(&pb, &other, v_dest, v_origin);
                let (pa_lo, pa_hi) = pa.parameter_range();
                let (pb_lo, pb_hi) = pb.parameter_range();
                for i in 0..n_samples {
                    let s = 0.05 + 0.90 * (i as f64) / ((n_samples - 1) as f64);
                    let s_a = if pa_fwd { s } else { 1.0 - s };
                    // Twin runs the other way along the edge, so its
                    // forward direction at sample s corresponds to the
                    // cap-side's `1 - s` already; combined with its
                    // own pcurve direction:
                    let s_b = if pb_fwd { 1.0 - s } else { s };
                    let t_a = pa_lo + s_a * (pa_hi - pa_lo);
                    let t_b = pb_lo + s_b * (pb_hi - pb_lo);
                    let uv_a = pa.point_at(t_a);
                    let uv_b = pb.point_at(t_b);
                    let p_a = cap_surf.point_at(uv_a.x, uv_a.y);
                    let p_b = other.point_at(uv_b.x, uv_b.y);
                    let g0 = p_a.distance_to(&p_b);
                    if g0 > max_g0 {
                        max_g0 = g0;
                    }
                    let su_a = cap_surf.partial_u(uv_a.x, uv_a.y);
                    let sv_a = cap_surf.partial_v(uv_a.x, uv_a.y);
                    let su_b = other.partial_u(uv_b.x, uv_b.y);
                    let sv_b = other.partial_v(uv_b.x, uv_b.y);
                    let n_a = su_a.cross(&sv_a);
                    let n_b = su_b.cross(&sv_b);
                    let nn_a = (n_a.dot(&n_a)).sqrt();
                    let nn_b = (n_b.dot(&n_b)).sqrt();
                    if nn_a > 1e-300 && nn_b > 1e-300 {
                        let dot_unit = (n_a.dot(&n_b) / (nn_a * nn_b)).abs();
                        let deficit = 1.0 - dot_unit;
                        if deficit > max_g1 {
                            max_g1 = deficit;
                        }
                    }
                }
            }
        }
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            break;
        }
    }
    (max_g0, max_g1)
}

#[test]
#[ignore]
fn cube_uniform_apex_no_rebuild_baseline() {
    let (brep, cap_face) = fillet_uniform_with_rebuild(false);
    brep.validate().expect("baseline BRep validates");
    let n_edges = cap_loop_edge_count(&brep, cap_face);
    let (g0, g1) = cap_vs_cylinder_g0_g1(&brep, cap_face);
    let surf_kind = brep.faces[cap_face]
        .surface
        .as_ref()
        .map(|s| match &s.geometry {
            SurfaceGeometry::Sphere { .. } => "Sphere",
            SurfaceGeometry::Gregory(_) => "Gregory",
            SurfaceGeometry::Plane { .. } => "Plane",
            SurfaceGeometry::Cylinder { .. } => "Cyl",
            SurfaceGeometry::Nurbs(_) => "Nurbs",
            _ => "Other",
        })
        .unwrap_or("?");
    eprintln!(
        "BASELINE (no rebuild): cap_face={cap_face:?} surf={surf_kind} \
         cap_loop_edges={n_edges} \
         BRep V={} E={} F={}, max_g0={g0:.3e} max_g1={g1:.3e}",
        brep.vertices.len(),
        brep.edges.len(),
        brep.faces.len(),
    );
}

#[test]
#[ignore]
fn cube_uniform_apex_rebuild_e2e() {
    let (brep, cap_face) = fillet_uniform_with_rebuild(true);
    brep.validate().expect("post-rebuild BRep validates");
    let n_edges = cap_loop_edge_count(&brep, cap_face);
    let (g0, g1) = cap_vs_cylinder_g0_g1(&brep, cap_face);
    // Diagnostic: walk cap loop and report per-edge state.
    let outer_loop = brep.faces[cap_face].outer_loop;
    let start = brep.loops[outer_loop].half_edge;
    let mut cursor = start;
    let mut idx = 0;
    loop {
        let he = &brep.half_edges[cursor];
        let twin = &brep.half_edges[he.twin];
        let twin_face = brep.loops[twin.loop_].face;
        let twin_kind = brep.faces[twin_face]
            .surface
            .as_ref()
            .map(|s| match &s.geometry {
                SurfaceGeometry::Sphere { .. } => "Sph",
                SurfaceGeometry::Cylinder { .. } => "Cyl",
                SurfaceGeometry::Plane { .. } => "Pln",
                _ => "?",
            })
            .unwrap_or("?");
        let curve_kind = brep.edges[he.edge]
            .curve
            .as_ref()
            .map(|c| {
                use popkern_core::geometry::CurveGeometry;
                match &c.geometry {
                    CurveGeometry::Nurbs(n) => format!("Nurbs(deg={})", n.degree),
                    CurveGeometry::Line { .. } => "Line".to_string(),
                    CurveGeometry::Arc { .. } => "Arc".to_string(),
                    CurveGeometry::Circle { .. } => "Circle".to_string(),
                    CurveGeometry::Ellipse { .. } => "Ellipse".to_string(),
                    _ => "Other".to_string(),
                }
            })
            .unwrap_or_else(|| "None".to_string());
        let pa_kind = if he.pcurve.is_some() { "Some" } else { "None" };
        let pb_kind = if twin.pcurve.is_some() {
            "Some"
        } else {
            "None"
        };
        eprintln!(
            "  edge[{idx}]: {:?} twin_face={twin_face:?}({twin_kind}) curve={curve_kind} pa={pa_kind} pb={pb_kind}",
            he.edge
        );
        cursor = he.next;
        idx += 1;
        if cursor == start || idx >= 20 {
            break;
        }
    }
    eprintln!(
        "REBUILD ON: cap_loop_edges={n_edges} \
         max_g0={g0:.3e} max_g1={g1:.3e}"
    );
    // Phase 6c-C topology gate: cap loop reduces 6 → 3 edges
    // (the sphere-fit triangle).  This is the hard guarantee.
    assert_eq!(n_edges, 3, "Phase 6c-C: cap loop should reduce 6 → 3 edges");
    // G0/G1 metric is descriptive only at this stage — the
    // pcurves on the inserted arcs are degree-5 NURBS fits
    // through 24 samples of (project(arc(s)) onto surface). Fit
    // residuals can leave a few-cm gap on the sphere chart's
    // non-trivial parameterisation; commercial-grade close-out
    // requires either an analytic small-circle pcurve (Phase
    // 6c-C.2e follow-up) or a higher-fidelity fit. Tracked
    // separately. Empirical baseline: ~9 cm on cube_uniform
    // (down from 13 cm pre-pcurve-fit).
}
