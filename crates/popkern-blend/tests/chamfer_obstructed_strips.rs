//! PK-1003: a chamfer whose strip runs into another feature declines
//! instead of serving a wrong body.
//!
//! A boss standing on the plate's top face 0.05 behind the chamfered
//! edge: once the setback reaches the boss, the per-edge trims ran the
//! bevel under it and served validator-clean bodies LARGER than the
//! input (d = 0.06: +0.0158), whether the boss's footprint is a hole in
//! the top face or a notch in its boundary. The strip gate declines
//! those edges (a rail through the boss's corner included), bevel
//! growth's sweep guard declines growing into the boss, and a request
//! with one declined edge declines whole (the pipeline used to drop the
//! edge and serve the rest).

#[path = "support/chamfer_fixtures.rs"]
#[allow(dead_code)] // not every fixture is used here
mod chamfer_fixtures;
#[path = "support/mesh_gauges.rs"]
#[allow(dead_code)] // the welded-mesh gauge is not used here
mod mesh_gauges;

use chamfer_fixtures::*;
use popkern_blend::BlendError;
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

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

/// The 4 × 2 × 0.3 plate with a 0.5-tall boss on top from y = 0.05 to
/// `back` (0.5: a hole in the top face; 2.0: a notch in its boundary).
fn plate_with_boss(back: f64) -> BRep {
    let boss = prism(&rect(1.5, 0.05, 2.5, back), 0.3, 0.5);
    popkern_bool::union(
        prism(&rect(0.0, 0.0, 4.0, 2.0), 0.0, 0.3),
        boss,
        Tolerance::default(),
        None,
    )
    .expect("plate with boss")
}

fn assert_declines_whole(label: &str, body: &BRep, edges: &[(Point3, Point3)], d: f64) {
    match chamfer(body, edges, d) {
        Err(BlendError::PartialFailure { partial_result, .. }) => assert_eq!(
            partial_result.faces.len(),
            body.faces.len(),
            "{label}: the input comes back"
        ),
        Ok(out) => panic!(
            "{label}: served (volume {} from {})",
            mesh_gauges::mesh_volume(&out),
            mesh_gauges::mesh_volume(body)
        ),
        Err(e) => panic!("{label}: unexpected error {e}"),
    }
}

#[test]
fn bevel_into_a_boss_declines() {
    for (name, back) in [("hole", 0.5), ("notch", 2.0)] {
        let body = plate_with_boss(back);
        let input = mesh_gauges::mesh_volume(&body);
        // Short of the boss: the plain wedge, d²/2 along 4.
        let d = 0.04;
        let out = chamfer(&body, &[TOP_FRONT], d).expect("short of the boss");
        assert_exact(name, &out, input - 2.0 * d * d);
        for d in [0.05, 0.06, 0.1] {
            assert_declines_whole(&format!("{name}, d={d}"), &body, &[TOP_FRONT], d);
        }
    }
}

/// Top and bottom front edges past the boss (the measured case), and a
/// request mixing the obstructed edge with a free one: all decline.
#[test]
fn one_obstructed_edge_declines_the_request() {
    let body = plate_with_boss(0.5);
    let bottom_front = (p(0.0, 0.0, 0.0), p(4.0, 0.0, 0.0));
    let bottom_back = (p(0.0, 2.0, 0.0), p(4.0, 2.0, 0.0));
    for d in [0.06, 0.08, 0.1, 0.14] {
        assert_declines_whole(
            &format!("top and bottom, d={d}"),
            &body,
            &[TOP_FRONT, bottom_front],
            d,
        );
        assert_declines_whole(
            &format!("with a free edge, d={d}"),
            &body,
            &[TOP_FRONT, bottom_back],
            d,
        );
    }
}
