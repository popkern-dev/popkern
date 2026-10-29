//! PK-997: chamfers at a convex vertex of valence four or more that
//! leave an edge sharp.
//!
//! A planar bevel is parallel to its own edge, so at a convex vertex it
//! crosses every other edge, not only its two neighbours: the vertex
//! is cut off and the bevel reaches the far faces. Before PK-997 every
//! such corner declined (check 9) unless symmetry made the planes meet
//! where the per-edge pipeline left the vertex. The truncated apex
//! stage settles the vertex at its face planes' point or ridge tree.
//!
//! Every gate is numeric: validation (structure, check 9, open loops,
//! self-intersection) plus the exact volume from the independent
//! half-space clip (a planar chamfer of a convex body is its
//! intersection with the bevel half-spaces). Irregular cones (apex off
//! centre, unequal face slopes) keep symmetry from hiding a wrong
//! position; a scale and rotation sweep keeps the construction from
//! being fixture-specific.

#[path = "support/apex_cones.rs"]
mod apex_cones;
#[path = "support/chamfer_fixtures.rs"]
#[allow(dead_code)] // not every fixture is used here
mod chamfer_fixtures;
#[path = "support/halfspace_oracle.rs"]
mod halfspace_oracle;
#[path = "support/mesh_gauges.rs"]
#[allow(dead_code)] // the welded-mesh gauge is not used here
mod mesh_gauges;

use apex_cones::*;
use chamfer_fixtures::*;
use popkern_blend::{BlendResult, ChamferKind, ChamferOptions, ChamferSpec, chamfer_edges};
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn chamfer_report(body: &BRep, edges: &[(Point3, Point3)], d: f64) -> BlendResult {
    let specs: Vec<ChamferSpec> = edges
        .iter()
        .map(|&(a, b)| ChamferSpec {
            edge: edge(body, a, b),
            kind: ChamferKind::Equal(d),
        })
        .collect();
    chamfer_edges(
        body,
        &specs,
        &ChamferOptions::default(),
        Tolerance::default(),
        None,
    )
    .unwrap_or_else(|e| panic!("declined: {e}"))
}

/// Chamfer `selection` of the apex edges of `body`: exact, and settled
/// by the truncated apex stage whenever the planes do not concur at
/// the old apex.
fn assert_apex_exact(label: &str, body: &BRep, selection: &[usize], d: f64) {
    let (apex, ends) = apex_edges(body);
    let edges: Vec<(Point3, Point3)> = selection.iter().map(|&i| (apex, ends[i])).collect();
    let result = chamfer_report(body, &edges, d);
    let exact = halfspace_oracle::chamfered_volume(body, &edges, d);
    assert_exact(label, &result.brep, exact);
}

#[test]
fn regular_pyramid_apexes() {
    for n in 4..=6usize {
        let body = BRep::n_gonal_pyramid(n, 1.0, 1.2, Tolerance::default()).expect("pyramid");
        for (name, selection) in partial_selections(n) {
            assert_apex_exact(&format!("{n}-gonal pyramid, {name}"), &body, &selection, D);
        }
    }
}

#[test]
fn irregular_cone_apexes() {
    for (k, body) in [(4, cone4()), (5, cone5()), (6, cone6())] {
        for d in [0.05, 0.1] {
            for (name, selection) in partial_selections(k) {
                assert_apex_exact(&format!("{k}-cone d={d}, {name}"), &body, &selection, d);
            }
        }
    }
}

/// The stage is what serves these: one bevel at a 5-valent apex needs
/// a two-node ridge tree, and the counter says the corner was rebuilt.
#[test]
fn truncated_apex_is_reported() {
    let body = cone5();
    let (apex, ends) = apex_edges(&body);
    let result = chamfer_report(&body, &[(apex, ends[0])], D);
    assert_eq!(result.repair_report.truncated_apex_corners, 1);
    let around = chamfer_report(&body, &[(apex, ends[0]), (apex, ends[2])], D);
    assert_eq!(around.repair_report.truncated_apex_corners, 1);
}

#[test]
fn apex_scale_and_rotation_sweep() {
    let rotation =
        Transform3::rotation_around_axis(p(0.1, -0.2, 0.3), Vec3::new(1.0, 2.0, 3.0), 0.7)
            .expect("rotation");
    let base = cone5();
    let (apex, ends) = apex_edges(&base);
    for s in [1e-2, 1.0, 1e2] {
        let place = rotation.compose(&Transform3::scale_uniform(s));
        let body = base.transform(&place).expect("placed cone");
        let at = |q: Point3| place.apply_point(&q);
        for selection in [vec![0], vec![0, 1], vec![0, 2], vec![0, 1, 2, 3]] {
            let edges: Vec<(Point3, Point3)> =
                selection.iter().map(|&i| (at(apex), at(ends[i]))).collect();
            let out = chamfer(&body, &edges, D * s)
                .unwrap_or_else(|e| panic!("scale {s}, {selection:?}: {e}"));
            let exact = halfspace_oracle::chamfered_volume(
                &base,
                &edges_unplaced(apex, &ends, &selection),
                D,
            );
            assert_exact(
                &format!("5-cone {selection:?}, scale {s}"),
                &out,
                exact * s * s * s,
            );
        }
    }
}

fn edges_unplaced(apex: Point3, ends: &[Point3], selection: &[usize]) -> Vec<(Point3, Point3)> {
    selection.iter().map(|&i| (apex, ends[i])).collect()
}
