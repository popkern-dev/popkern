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

#[path = "support/chamfer_fixtures.rs"]
#[allow(dead_code)] // not every fixture is used here
mod chamfer_fixtures;
#[path = "support/halfspace_oracle.rs"]
mod halfspace_oracle;
#[path = "support/mesh_gauges.rs"]
#[allow(dead_code)] // the welded-mesh gauge is not used here
mod mesh_gauges;

use chamfer_fixtures::*;
use popkern_blend::{
    BlendError, BlendResult, ChamferKind, ChamferOptions, ChamferSpec, chamfer_edges,
};
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// The half-space `(x − a) · n ≤ 0` as a box with one face on the plane.
fn half_space(a: Point3, n: Vec3, size: f64) -> BRep {
    let n = n.normalized().expect("normal");
    let helper = if n.z.abs() < 0.9 {
        Vec3::new(0.0, 0.0, 1.0)
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let u = n.cross(&helper).normalized().expect("u");
    let v = n.cross(&u);
    let q = |su: f64, sv: f64| a + u * (su * size) + v * (sv * size);
    let outline = [q(-1.0, -1.0), q(1.0, -1.0), q(1.0, 1.0), q(-1.0, 1.0)];
    BRep::prism(&outline, -n, 2.0 * size, Tolerance::default()).expect("half-space box")
}

const APEX: (f64, f64, f64) = (0.13, -0.07, 1.1);

/// A convex cone with apex [`APEX`] cut from the slab `0 ≤ z ≤ 3` by
/// half-spaces whose outward normals have the given azimuths and
/// elevations (degrees): unequal slopes, apex off centre.
fn cone(faces: &[(f64, f64)]) -> BRep {
    let a = p(APEX.0, APEX.1, APEX.2);
    let mut body = prism(&rect(-3.0, -3.0, 3.0, 3.0), 0.0, 3.0);
    for &(az, el) in faces {
        let (az, el) = (az.to_radians(), el.to_radians());
        let n = Vec3::new(el.cos() * az.cos(), el.cos() * az.sin(), el.sin());
        body = popkern_bool::intersection(body, half_space(a, n, 6.0), Tolerance::default(), None)
            .expect("cone cut");
    }
    body
}

fn cone4() -> BRep {
    cone(&[(5.0, 40.0), (95.0, 52.0), (170.0, 35.0), (265.0, 47.0)])
}

fn cone5() -> BRep {
    cone(&[
        (0.0, 38.0),
        (80.0, 50.0),
        (150.0, 33.0),
        (220.0, 45.0),
        (290.0, 41.0),
    ])
}

fn cone6() -> BRep {
    cone(&[
        (10.0, 42.0),
        (65.0, 36.0),
        (130.0, 50.0),
        (185.0, 39.0),
        (250.0, 44.0),
        (310.0, 34.0),
    ])
}

/// The apex (the highest vertex) and the far ends of its edges, in
/// angular order about it.
fn apex_edges(body: &BRep) -> (Point3, Vec<Point3>) {
    let apex = body
        .vertices
        .values()
        .map(|v| v.position)
        .max_by(|a, b| a.z.total_cmp(&b.z))
        .expect("apex");
    let mut ends: Vec<Point3> = body
        .edges
        .values()
        .filter_map(|e| {
            let he = &body.half_edges[e.half_edge];
            let o = body.vertices[he.origin].position;
            let d = body.vertices[body.half_edges[he.twin].origin].position;
            if o.distance_to(&apex) < 1e-9 {
                Some(d)
            } else if d.distance_to(&apex) < 1e-9 {
                Some(o)
            } else {
                None
            }
        })
        .collect();
    let angle = |q: &Point3| (q.y - apex.y).atan2(q.x - apex.x);
    ends.sort_by(|a, b| angle(a).total_cmp(&angle(b)));
    (apex, ends)
}

/// Every apex selection that leaves an edge sharp, by shape: each
/// single edge, each adjacent pair, a pair around one sharp edge, and
/// all but one.
fn partial_selections(valence: usize) -> Vec<(String, Vec<usize>)> {
    let mut out = Vec::new();
    for i in 0..valence {
        out.push((format!("edge {i}"), vec![i]));
        out.push((format!("pair {i}"), vec![i, (i + 1) % valence]));
    }
    out.push(("around a sharp edge".into(), vec![0, 2]));
    out.push(("all but one".into(), (0..valence - 1).collect()));
    out
}

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

/// Large bevels on the 6-valent cone: two adjacent bevels can consume
/// the whole face between them, which the per-edge pipeline cannot
/// delete (PK-1001). Those decline cleanly; every other selection at
/// this distance serves exactly.
#[test]
fn face_consuming_bevels_decline_cleanly() {
    let body = cone6();
    let (apex, ends) = apex_edges(&body);
    let d = 0.3;
    let mut declined = 0;
    for (name, selection) in partial_selections(6) {
        let edges: Vec<(Point3, Point3)> = selection.iter().map(|&i| (apex, ends[i])).collect();
        match chamfer(&body, &edges, d) {
            Ok(out) => {
                let exact = halfspace_oracle::chamfered_volume(&body, &edges, d);
                assert_exact(&format!("6-cone d={d}, {name}"), &out, exact);
            }
            Err(BlendError::PartialFailure { partial_result, .. }) => {
                declined += 1;
                assert_eq!(
                    partial_result.faces.len(),
                    body.faces.len(),
                    "{name}: the input comes back"
                );
            }
            Err(e) => panic!("{name}: unexpected error {e}"),
        }
    }
    assert!(
        (1..=3).contains(&declined),
        "{declined} face-consuming selections declined"
    );
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
