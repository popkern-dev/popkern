//! PK-375 — terminal ISO-CURVE settle for curvature-continuous
//! (G2) swept strips: hard gates.
//!
//! A G2 (degree-5 section) single-edge fillet on the unit cube has a
//! terminal section that is neither circular nor conic, so both
//! analytic settle arms decline; before PK-375 the chain-end corner
//! stayed at the cube corner and check 9 reported a 3.6e-2
//! curve-vs-surface lie on every stub edge (`validate()` passed,
//! `validate_all()` did not — `Part::fillet_with(G2)` never served on
//! a cube). The strip's exact terminal iso curve is the section; the
//! corner and both stubs settle onto it.

use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, SurfaceEval};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn g2_single_edge_fillet(radius: f64) -> popkern_blend::BlendResult {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube");
    let edge = brep.edges.keys().next().expect("edge");
    fillet_edges(
        &brep,
        &[FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        }],
        &FilletOptions {
            continuity: ContinuityLevel::G2,
            ..FilletOptions::default()
        },
        Tolerance::default(),
        None,
    )
    .expect("G2 single-edge fillet")
}

/// The ends close exactly with no declines across the radius range,
/// and the body passes the full validator (check 9 included). Since
/// PK-531 the G2 strip reaches past the edge's end stations, so the
/// face trim cuts the exact end section directly: no stale corner is
/// left for the iso settle, which remains the settlement for fitted
/// strips that still end at their stations.
#[test]
fn g2_single_edge_fillet_settles_both_corners_and_validates_all() {
    for radius in [0.05, 0.1, 0.3, 0.5] {
        let result = g2_single_edge_fillet(radius);
        let report = &result.repair_report.terminal_arc_trim;
        assert_eq!(
            report.repaired_corners, 0,
            "r={radius}: the extended strip leaves no stale corner (got {report:?})"
        );
        assert!(
            report.declines.is_empty(),
            "r={radius}: no declines expected (got {report:?})"
        );
        let errors = result.brep.validate_all();
        assert!(
            errors.is_empty(),
            "r={radius}: G2 fillet must pass validate_all: {errors:?}"
        );
    }
}

/// Every vertex on the G2 strip's loop lies ON the strip surface and
/// every stub edge's 3-D curve traces the strip at interior stations
/// and endpoints — one geometric authority, exact to the repair band.
#[test]
fn g2_strip_boundary_geometry_is_on_surface() {
    let result = g2_single_edge_fillet(0.1);
    let brep = &result.brep;
    let tol = Tolerance::default().linear() * 10.0;
    let (strip, surface) = brep
        .faces
        .iter()
        .find_map(|(h, f)| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(s @ SurfaceGeometry::Nurbs(_)) => Some((h, s.clone())),
            _ => None,
        })
        .expect("G2 strip is a NURBS face");
    let mut checked = 0;
    for (_, edge) in &brep.edges {
        let he = edge.half_edge;
        let twin = brep.half_edges[he].twin;
        let faces = [
            brep.loops[brep.half_edges[he].loop_].face,
            brep.loops[brep.half_edges[twin].loop_].face,
        ];
        if !faces.contains(&strip) {
            continue;
        }
        let curve = edge.curve.as_ref().expect("strip edge carries a curve");
        let (t0, t1) = if curve.t_start.is_finite() && curve.t_end.is_finite() {
            (curve.t_start, curve.t_end)
        } else {
            continue;
        };
        for k in 0..=8 {
            let t = t0 + (t1 - t0) * k as f64 / 8.0;
            let p = curve.geometry.point_at(t);
            let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(&surface, &p)
                .expect("projects onto the strip");
            let d = (surface.point_at(u, v) - p).length();
            assert!(
                d <= tol,
                "strip-loop curve sample {p:?} is {d:.3e} off the G2 strip (band {tol:.1e})"
            );
            checked += 1;
        }
        for vh in [brep.half_edges[he].origin, brep.half_edges[twin].origin] {
            let p = brep.vertices[vh].position;
            let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(&surface, &p)
                .expect("projects onto the strip");
            let d = (surface.point_at(u, v) - p).length();
            assert!(
                d <= tol,
                "strip-loop vertex {p:?} is {d:.3e} off the G2 strip (band {tol:.1e})"
            );
        }
    }
    assert!(
        checked >= 4 * 9,
        "expected the four stub edges to be sampled"
    );
}
