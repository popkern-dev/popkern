//! PK-148 — terminal arc trim at open fillet chain ends: hard gates.
//!
//! A single-edge fillet on a cube terminates at two plain body
//! vertices. The strip's terminal boundary must be the exact
//! cross-section arc: corner vertex ON the strip surface, stub 3-D
//! curves tracing the arc, cap-side pcurves agreeing — one geometric
//! authority. Before PK-148 the stubs carried straight line curves
//! through the original corner ((√2−1)·r off the arc) and STEP
//! export forked face-local edges (8–12 OCC free closed edges).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, SurfaceEval};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "support/mesh_gauges.rs"]
mod mesh_gauges;

fn default_tol() -> Tolerance {
    Tolerance::default()
}

fn single_edge_fillet(radius: f64) -> popkern_fillet::BlendResult {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    fillet_edges(
        &brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        }],
        &FilletOptions::default(),
        default_tol(),
        None,
    )
    .expect("single-edge fillet")
}

/// Both chain-end corners are settled, with no declines, across the
/// feasible radius range.
#[test]
fn single_edge_fillet_settles_both_corners() {
    for radius in [0.05, 0.1, 0.3, 0.5, 0.7, 0.9] {
        let result = single_edge_fillet(radius);
        let report = &result.repair_report.terminal_arc_trim;
        assert_eq!(
            report.repaired_corners, 2,
            "r={radius}: both chain-end corners must be arc-settled (got {report:?})"
        );
        assert!(
            report.declines.is_empty(),
            "r={radius}: no declines expected (got {report:?})"
        );
    }
}

/// Every vertex on the strip's outer loop lies ON the strip surface,
/// and every strip-loop edge's 3-D curve midpoint lies on the strip
/// surface too (single geometric authority — no chord-vs-arc split).
#[test]
fn strip_boundary_geometry_is_on_surface() {
    let result = single_edge_fillet(0.3);
    let brep = &result.brep;
    let (face, f) = brep
        .faces
        .iter()
        .find(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Cylinder { .. }))
        })
        .expect("fillet strip");
    let surface = &f.surface.as_ref().expect("surface").geometry;
    let l = &brep.loops[f.outer_loop];
    let mut heh = l.half_edge;
    let mut checked_vertices = 0usize;
    let mut checked_curves = 0usize;
    loop {
        let he = &brep.half_edges[heh];
        let pos = brep.vertices[he.origin].position;
        // Vertex-on-surface via the strip pcurve seed when available.
        if let Some(pc) = brep.face_use_pcurve(face, heh).or(he.pcurve.as_ref()) {
            let (t0, _) = pc.parameter_range();
            let uv = pc.point_at(t0);
            let on = surface.point_at(uv.x, uv.y);
            assert!(
                (on - pos).length() < 1e-9,
                "loop vertex {pos:?} is {} off the strip surface",
                (on - pos).length()
            );
            checked_vertices += 1;
        }
        if let Some(c) = brep.edges[he.edge].curve.as_ref() {
            let mid = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
            // Every strip boundary curve must trace ON the cylinder:
            // |dist(mid, axis) − r| ≈ 0.
            if let SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } = surface
            {
                let d = mid - *origin;
                let radial = d - *axis * d.dot(axis);
                assert!(
                    (radial.length() - radius).abs() < 1e-9,
                    "edge {:?} curve midpoint {mid:?} off-cylinder by {}",
                    he.edge,
                    (radial.length() - radius).abs()
                );
                checked_curves += 1;
            }
        }
        heh = he.next;
        if heh == l.half_edge {
            break;
        }
    }
    assert!(checked_vertices >= 6, "loop walk covered the boundary");
    assert!(checked_curves >= 6, "every strip edge carries a 3-D curve");
}

/// The completed body is BRep-clean AND its tessellation is fully
/// conforming under the welded gauge (the pre-PK-148 state read
/// 26–34 welded opens from the chord-vs-arc split).
#[test]
fn single_edge_fillet_mesh_is_conforming() {
    for radius in [0.1, 0.3, 0.6] {
        let result = single_edge_fillet(radius);
        let (structural, open, other) = mesh_gauges::brep_open_defects(&result.brep);
        assert_eq!(
            (structural, open, other),
            (0, 0, 0),
            "r={radius}: completed body must validate clean"
        );
        let opens = mesh_gauges::welded_open_edges(&result.brep);
        assert_eq!(
            opens, 0,
            "r={radius}: tessellation must be conforming (welded opens)"
        );
    }
}

/// Idempotence: re-running identify/settlement machinery through a
/// second fillet on the same body must not disturb settled corners
/// (the recognizer skips on-surface vertices).
#[test]
fn settlement_is_idempotent_across_second_fillet() {
    let first = single_edge_fillet(0.2);
    // Fillet a second, disjoint edge of the result.
    let edge = first
        .brep
        .edges
        .iter()
        .find(|(_, d)| {
            let he = &first.brep.half_edges[d.half_edge];
            let a = first.brep.vertices[he.origin].position;
            let b = first.brep.vertices[first.brep.half_edges[he.twin].origin].position;
            // Top far edge — untouched by the first fillet.
            (a.z - 1.0).abs() < 1e-9
                && (b.z - 1.0).abs() < 1e-9
                && (a.y - 1.0).abs() < 1e-9
                && (b.y - 1.0).abs() < 1e-9
        })
        .map(|(h, _)| h)
        .expect("disjoint top edge");
    let second = fillet_edges(
        &first.brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.2),
        }],
        &FilletOptions::default(),
        default_tol(),
        None,
    )
    .expect("second fillet");
    let report = &second.repair_report.terminal_arc_trim;
    assert_eq!(
        report.repaired_corners, 2,
        "second fillet settles exactly its own two corners: {report:?}"
    );
    let (structural, open, other) = mesh_gauges::brep_open_defects(&second.brep);
    assert_eq!((structural, open, other), (0, 0, 0));
}
