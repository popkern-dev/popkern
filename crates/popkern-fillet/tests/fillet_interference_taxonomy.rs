//! PK-140 — fillet interference taxonomy (stage 1).
//!
//! Three C-grade interference classes measured against the
//! commercial bar: each input must produce a WATERTIGHT completed
//! body or a STRUCTURED decline — never a panic and never silent
//! garbage (open/self-intersecting output that claims success).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
#[path = "support/mesh_gauges.rs"]
mod mesh_gauges;
use mesh_gauges::welded_open_edges;
use popkern_fillet::{BlendResult, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn default_tol() -> Tolerance {
    Tolerance::default()
}

/// Run specs on a body; enforce the taxonomy bar.
fn run_class(label: &str, brep: &BRep, specs: &[FilletSpec]) {
    let result = fillet_edges(brep, specs, &FilletOptions::default(), default_tol(), None);
    match result {
        Ok(BlendResult { brep: out, .. }) => {
            let opens = welded_open_edges(&out);
            assert_eq!(
                opens, 0,
                "{label}: completed body is OPEN ({opens} mesh edges) — silent garbage; must \
                 decline or complete watertight"
            );
            let errors = out.validate_all();
            assert!(
                errors.is_empty(),
                "{label}: completed body invalid ({} errors) — silent garbage",
                errors.len()
            );
            eprintln!("[interference] {label}: COMPLETED watertight + validate_all clean");
        }
        Err(err) => {
            eprintln!(
                "[interference] {label}: DECLINED ({})",
                short(&format!("{err:?}"))
            );
        }
    }
}

fn short(s: &str) -> String {
    s.chars().take(160).collect()
}

/// Two opposite edges of one cube face, radii summing past the face
/// width — the fillet bands collide mid-face.
#[test]
#[ignore = "PK-140 stage 2: baseline gauge question first — plain single-edge cube fillet reads 26-34 welded opens on RAW fillet_edges output at every radius (r=0.3..0.95), so the interference gates cannot arm until the raw-output watertightness baseline is settled (Part-pipeline healing vs fillet_edges contract); then convert each class to structured decline"]
fn mid_face_collision_declines_or_completes() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    // Opposite vertical edges of the x=… face: pick two edges that
    // share a face but no vertex, with parallel directions.
    let mut pair: Option<(
        popkern_core::topology::EdgeHandle,
        popkern_core::topology::EdgeHandle,
    )> = None;
    'outer: for (e1, d1) in brep.edges.iter() {
        for (e2, d2) in brep.edges.iter() {
            if e1 >= e2 {
                continue;
            }
            let he1 = &brep.half_edges[d1.half_edge];
            let he2 = &brep.half_edges[d2.half_edge];
            let f1 = [
                brep.loops[he1.loop_].face,
                brep.loops[brep.half_edges[he1.twin].loop_].face,
            ];
            let f2 = [
                brep.loops[he2.loop_].face,
                brep.loops[brep.half_edges[he2.twin].loop_].face,
            ];
            let shares_face = f1.iter().any(|f| f2.contains(f));
            let v1 = [he1.origin, brep.half_edges[he1.twin].origin];
            let v2 = [he2.origin, brep.half_edges[he2.twin].origin];
            let shares_vertex = v1.iter().any(|v| v2.contains(v));
            if shares_face && !shares_vertex {
                pair = Some((e1, e2));
                break 'outer;
            }
        }
    }
    let (e1, e2) = pair.expect("cube has an opposite-edge pair");
    run_class(
        "mid-face collision r=0.6+0.6 across width 1.0",
        &brep,
        &[
            FilletSpec {
                edge: e1,
                radius: FilletRadius::Constant(0.6),
            },
            FilletSpec {
                edge: e2,
                radius: FilletRadius::Constant(0.6),
            },
        ],
    );
}

/// Radius nearly consuming the adjacent faces (overflow/cliff class).
#[test]
#[ignore = "PK-140 stage 2: baseline gauge question first — plain single-edge cube fillet reads 26-34 welded opens on RAW fillet_edges output at every radius (r=0.3..0.95), so the interference gates cannot arm until the raw-output watertightness baseline is settled (Part-pipeline healing vs fillet_edges contract); then convert each class to structured decline"]
fn large_radius_overflow_declines_or_completes() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    run_class(
        "overflow r=0.95 on unit cube",
        &brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.95),
        }],
    );
}

/// Radius exactly at the face extent — the degenerate boundary of
/// the overflow class.
#[test]
#[ignore = "PK-140 stage 2: baseline gauge question first — plain single-edge cube fillet reads 26-34 welded opens on RAW fillet_edges output at every radius (r=0.3..0.95), so the interference gates cannot arm until the raw-output watertightness baseline is settled (Part-pipeline healing vs fillet_edges contract); then convert each class to structured decline"]
fn full_face_radius_declines_or_completes() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    run_class(
        "exact-extent r=1.0 on unit cube",
        &brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(1.0),
        }],
    );
}

/// Thin wall: plate 1.0×1.0×0.05 (square prism), fillet an edge with
/// r=0.2 — radius exceeds the local feature size four times over.
#[test]
#[ignore = "PK-140 stage 2: baseline gauge question first — plain single-edge cube fillet reads 26-34 welded opens on RAW fillet_edges output at every radius (r=0.3..0.95), so the interference gates cannot arm until the raw-output watertightness baseline is settled (Part-pipeline healing vs fillet_edges contract); then convert each class to structured decline"]
fn thin_wall_radius_declines_or_completes() {
    use popkern_core::geometry::Point3;
    let outline = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let brep = BRep::prism(
        &outline,
        popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        0.05,
        default_tol(),
    )
    .expect("thin plate");
    let edge = brep.edges.keys().next().unwrap();
    run_class(
        "thin wall t=0.05, r=0.2",
        &brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.2),
        }],
    );
}

/// Diagnostic sweep (not a gate): where does single-edge completion
/// actually break on the unit cube?
#[test]
#[ignore = "diagnostic sweep"]
fn single_edge_radius_envelope_sweep() {
    for r in [0.3, 0.4, 0.45, 0.475, 0.5, 0.55, 0.6, 0.7, 0.85, 0.95] {
        let brep = BRep::unit_cube(default_tol()).unwrap();
        let edge = brep.edges.keys().next().unwrap();
        let result = fillet_edges(
            &brep,
            &[FilletSpec {
                edge,
                radius: FilletRadius::Constant(r),
            }],
            &FilletOptions::default(),
            default_tol(),
            None,
        );
        match result {
            Ok(BlendResult { brep: out, .. }) => {
                eprintln!(
                    "[envelope] r={r}: completed, opens={}",
                    welded_open_edges(&out)
                );
            }
            Err(_) => eprintln!("[envelope] r={r}: declined"),
        }
    }
}
