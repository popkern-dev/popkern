//! PK-1001: planar chamfers whose bevels consume a whole face.
//!
//! The one-pass pipeline trims every face by its rails and cannot
//! delete one, so it declines; the request is served at a fraction of
//! its size and the bevels are grown to their planes, deleting each
//! consumed face (a triangle shrinks to a point, a thin wall's side to
//! a ridge) and flipping each edge a bevel sweeps past, at the exact
//! instants. Every test asserts the growth did the work (its counters)
//! and gates on validation plus an exact volume: the independent
//! half-space clip for convex bodies, an exact cross-section area for
//! prismatic ones (a non-convex L and a two-distance wall).

#[path = "support/apex_cones.rs"]
#[allow(dead_code)] // the 4-valent cone is not used here
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
use popkern_blend::{
    BlendError, BlendResult, ChamferKind, ChamferOptions, ChamferSpec, chamfer_edges,
};
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

fn served(body: &BRep, specs: Vec<ChamferSpec>, label: &str) -> BlendResult {
    chamfer_edges(
        body,
        &specs,
        &ChamferOptions::default(),
        Tolerance::default(),
        None,
    )
    .unwrap_or_else(|e| panic!("{label}: declined: {e}"))
}

fn equal(body: &BRep, edges: &[(Point3, Point3)], d: f64) -> Vec<ChamferSpec> {
    edges
        .iter()
        .map(|&(a, b)| ChamferSpec {
            edge: edge(body, a, b),
            kind: ChamferKind::Equal(d),
        })
        .collect()
}

/// Area of a simple polygon.
fn shoelace(poly: &[(f64, f64)]) -> f64 {
    let n = poly.len();
    let twice: f64 = (0..n)
        .map(|i| poly[i].0 * poly[(i + 1) % n].1 - poly[(i + 1) % n].0 * poly[i].1)
        .sum();
    0.5 * twice.abs()
}

/// The 4 × 2 × 0.3 plate and its front (y = 0) edges.
fn plate() -> BRep {
    prism(&rect(0.0, 0.0, 4.0, 2.0), 0.0, 0.3)
}
const TOP_FRONT: (Point3, Point3) = (
    Point3 {
        x: 0.0,
        y: 0.0,
        z: 0.3,
    },
    Point3 {
        x: 4.0,
        y: 0.0,
        z: 0.3,
    },
);
const BOTTOM_FRONT: (Point3, Point3) = (
    Point3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    },
    Point3 {
        x: 4.0,
        y: 0.0,
        z: 0.0,
    },
);
const TOP_LEFT: (Point3, Point3) = (
    Point3 {
        x: 0.0,
        y: 0.0,
        z: 0.3,
    },
    Point3 {
        x: 0.0,
        y: 2.0,
        z: 0.3,
    },
);

/// Every partial apex selection of the 5- and 6-valent cones at bevel
/// sizes that consume faces serves exactly; the 6-valent cone's
/// adjacent pair at d = 0.3 is served by consuming the narrow face
/// between them.
#[test]
fn cone_bevels_consume_faces_exactly() {
    let mut grown = 0;
    for (k, body, sizes) in [(5, cone5(), [0.45, 0.6]), (6, cone6(), [0.3, 0.45])] {
        let (apex, ends) = apex_edges(&body);
        for d in sizes {
            for (name, selection) in partial_selections(k) {
                let edges: Vec<(Point3, Point3)> =
                    selection.iter().map(|&i| (apex, ends[i])).collect();
                let label = format!("{k}-cone d={d}, {name}");
                let result = served(&body, equal(&body, &edges, d), &label);
                assert_exact(
                    &label,
                    &result.brep,
                    halfspace_oracle::chamfered_volume(&body, &edges, d),
                );
                grown += usize::from(result.repair_report.bevel_consumed_faces > 0);
            }
        }
    }
    assert!(grown >= 10, "only {grown} selections consumed a face");
    let body = cone6();
    let (apex, ends) = apex_edges(&body);
    let result = served(
        &body,
        equal(&body, &[(apex, ends[2]), (apex, ends[3])], 0.3),
        "pair",
    );
    assert_eq!(result.repair_report.bevel_consumed_faces, 1);
    assert_eq!(result.brep.faces.len(), 8, "6 slants + base + 2 bevels − 1");
}

/// A thin wall chamfered top and bottom past half its thickness: the
/// front face shrinks to a segment, which becomes the ridge between the
/// bevels. With the top-left edge too, the corner where the two top
/// bevels meet the sharp faces (four planes) splits as the face goes.
#[test]
fn thin_wall_side_becomes_a_ridge() {
    let body = plate();
    for d in [0.2, 0.25] {
        for (name, edges) in [
            ("top and bottom", vec![TOP_FRONT, BOTTOM_FRONT]),
            ("and top left", vec![TOP_FRONT, BOTTOM_FRONT, TOP_LEFT]),
        ] {
            let label = format!("plate d={d}, {name}");
            let result = served(&body, equal(&body, &edges, d), &label);
            assert_eq!(result.repair_report.bevel_consumed_faces, 1, "{label}");
            assert_exact(
                &label,
                &result.brep,
                halfspace_oracle::chamfered_volume(&body, &edges, d),
            );
        }
    }
}

/// Two-distance bevels (0.2 down the front, 0.1 across top and bottom)
/// consume the 0.3 front face; their ridge is at y = 0.025, z = 0.15.
#[test]
fn two_distance_bevels_consume_a_face() {
    let body = plate();
    let front: FaceHandle = body
        .faces
        .keys()
        .find(|&f| {
            body.face_outer_loop_positions(f)
                .iter()
                .all(|q| q.y.abs() < 1e-12)
        })
        .expect("front face");
    let kind = ChamferKind::TwoDistance {
        face_1: front,
        dist_1: 0.2,
        dist_2: 0.1,
    };
    let specs = [TOP_FRONT, BOTTOM_FRONT]
        .iter()
        .map(|&(a, b)| ChamferSpec {
            edge: edge(&body, a, b),
            kind: kind.clone(),
        })
        .collect();
    let result = served(&body, specs, "two-distance");
    assert_eq!(result.repair_report.bevel_consumed_faces, 1);
    let section = [
        (0.1, 0.0),
        (2.0, 0.0),
        (2.0, 0.3),
        (0.1, 0.3),
        (0.025, 0.15),
    ];
    assert_exact("two-distance", &result.brep, 4.0 * shoelace(&section));
}

/// A non-convex L section: both legs' 0.3 end faces consumed by d = 0.2
/// bevels, far from the concave inner corner.
#[test]
fn non_convex_body_ends_consumed() {
    let body = prism(
        &[
            (0.0, 0.0),
            (3.0, 0.0),
            (3.0, 0.3),
            (0.3, 0.3),
            (0.3, 2.0),
            (0.0, 2.0),
        ],
        0.0,
        4.0,
    );
    let upright = |x: f64, y: f64| (p(x, y, 0.0), p(x, y, 4.0));
    let edges = [
        upright(3.0, 0.0),
        upright(3.0, 0.3),
        upright(0.3, 2.0),
        upright(0.0, 2.0),
    ];
    let result = served(&body, equal(&body, &edges, 0.2), "L");
    assert_eq!(result.repair_report.bevel_consumed_faces, 2);
    let section = [
        (0.0, 0.0),
        (2.8, 0.0),
        (2.95, 0.15),
        (2.8, 0.3),
        (0.3, 0.3),
        (0.3, 1.8),
        (0.15, 1.95),
        (0.0, 1.8),
    ];
    assert_exact("L", &result.brep, 4.0 * shoelace(&section));
}

/// A block with a sloping top: the top-front bevel reaches the low
/// front corner first, where the sharp front-left edge vanishes and
/// flips into an edge between the bevel and the bottom face.
#[test]
fn bevel_flips_an_edge_it_sweeps_past() {
    let outline = [
        p(0.0, 0.0, 0.0),
        p(4.0, 0.0, 0.0),
        p(4.0, 0.0, 1.8),
        p(0.0, 0.0, 1.0),
    ];
    let body = BRep::prism(
        &outline,
        Vec3::new(0.0, 1.0, 0.0),
        2.0,
        Tolerance::default(),
    )
    .expect("sloping block");
    let top_front = [(p(0.0, 0.0, 1.0), p(4.0, 0.0, 1.8))];
    for d in [1.2, 1.6] {
        let label = format!("sloping block d={d}");
        let result = served(&body, equal(&body, &top_front, d), &label);
        assert_eq!(result.repair_report.bevel_growth_flips, 1, "{label}");
        assert_exact(
            &label,
            &result.brep,
            halfspace_oracle::chamfered_volume(&body, &top_front, d),
        );
    }
}

#[test]
fn consumed_face_scale_and_rotation_sweep() {
    let rotation =
        Transform3::rotation_around_axis(p(0.1, -0.2, 0.3), Vec3::new(1.0, 2.0, 3.0), 0.7)
            .expect("rotation");
    let cone = cone6();
    let (apex, ends) = apex_edges(&cone);
    let cases = [
        (
            "6-cone pair",
            cone,
            vec![(apex, ends[2]), (apex, ends[3])],
            0.3,
        ),
        ("plate", plate(), vec![TOP_FRONT, BOTTOM_FRONT], 0.2),
    ];
    for (name, base, edges, d) in cases {
        let exact = halfspace_oracle::chamfered_volume(&base, &edges, d);
        for s in [1e-2, 1.0, 1e2] {
            let place = rotation.compose(&Transform3::scale_uniform(s));
            let body = base.transform(&place).expect("placed body");
            let placed: Vec<(Point3, Point3)> = edges
                .iter()
                .map(|(a, b)| (place.apply_point(a), place.apply_point(b)))
                .collect();
            let label = format!("{name}, scale {s}");
            let result = served(&body, equal(&body, &placed, d * s), &label);
            assert_eq!(result.repair_report.bevel_consumed_faces, 1, "{label}");
            assert_exact(&label, &result.brep, exact * s * s * s);
        }
    }
}

/// A boss standing 0.15 behind the plate's front edge: growing the top
/// bevel from the served half size to d = 0.2 would sweep it under the
/// boss (a far feature, no edge event), so the request declines
/// cleanly with the input back (unguarded, a valid body of volume 2.011
/// served: less than the chamfered plate alone, 2.25).
#[test]
fn growth_into_a_far_feature_declines() {
    let boss = prism(&rect(1.5, 0.15, 2.5, 0.5), 0.3, 0.5);
    let body = popkern_bool::union(plate(), boss, Tolerance::default(), None).expect("union");
    let specs = equal(&body, &[TOP_FRONT, BOTTOM_FRONT], 0.2);
    match chamfer_edges(
        &body,
        &specs,
        &ChamferOptions::default(),
        Tolerance::default(),
        None,
    ) {
        Err(BlendError::PartialFailure { partial_result, .. }) => {
            assert_eq!(
                partial_result.faces.len(),
                body.faces.len(),
                "the input comes back"
            );
        }
        Ok(_) => panic!("served a bevel grown through the boss"),
        Err(e) => panic!("unexpected error {e}"),
    }
    // Short of the boss it serves: 0.575 + 2.4 less two d²/2 wedges.
    let near = served(&body, equal(&body, &[TOP_FRONT, BOTTOM_FRONT], 0.1), "near");
    assert_exact("short of the boss", &near.brep, 2.575 - 4.0 * 0.1 * 0.1);
}
