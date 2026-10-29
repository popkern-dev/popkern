//! PK-1002: a chamfer on every edge of one planar face (a closed ring)
//! serves on every face of a body, not only the face whose loop winding
//! happened to keep the corner edge first in the vertex ring walk.
//!
//! The ring surgery adds a spur at each corner before it re-authors the
//! corner edges; "the edge that is not a rim edge" then also matches
//! the spur, and on bottom and side faces the walk met the spur first,
//! leaving the corner edge's curve at the corner's old position (a
//! shell gap of exactly d). Rings outside the surgery's class (walls
//! meeting the cap at different angles, so the setbacks along a corner
//! edge differ) go to the per-edge path; rings whose bevels consume the
//! cap decline and grow from a smaller ring (PK-1001) instead of
//! serving an inside-out inset. Gates: validation and the exact
//! half-space volume, on every face of a box, a hexagonal prism and a
//! thin plate, and a scale and rotation sweep.

#[path = "support/chamfer_fixtures.rs"]
#[allow(dead_code)] // not every fixture is used here
mod chamfer_fixtures;
#[path = "support/halfspace_oracle.rs"]
mod halfspace_oracle;
#[path = "support/mesh_gauges.rs"]
#[allow(dead_code)] // the welded-mesh gauge is not used here
mod mesh_gauges;

use chamfer_fixtures::*;
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;

/// The outer-loop edges of every face, by endpoints.
fn face_rings(body: &BRep) -> Vec<Vec<(Point3, Point3)>> {
    body.faces
        .keys()
        .map(|f| {
            let ring = body.face_outer_loop_positions(f);
            (0..ring.len())
                .map(|i| (ring[i], ring[(i + 1) % ring.len()]))
                .collect()
        })
        .collect()
}

fn assert_every_ring_exact(label: &str, body: &BRep, d: f64) {
    for (k, ring) in face_rings(body).iter().enumerate() {
        let label = format!("{label}, face {k}");
        let out = chamfer(body, ring, d).unwrap_or_else(|e| panic!("{label}: declined: {e}"));
        assert_exact(
            &label,
            &out,
            halfspace_oracle::chamfered_volume(body, ring, d),
        );
    }
}

fn hexagonal_prism() -> BRep {
    let hexagon: Vec<(f64, f64)> = (0..6)
        .map(|k| {
            let a = (60.0 * f64::from(k) + 7.0).to_radians();
            (1.3 * a.cos() + 0.2, 1.3 * a.sin() - 0.1)
        })
        .collect();
    prism(&hexagon, 0.0, 1.5)
}

#[test]
fn every_face_ring_of_a_box_serves() {
    assert_every_ring_exact("block", &prism(&rect(0.0, 0.0, 4.0, 2.0), 0.0, 1.0), 0.1);
}

#[test]
fn every_face_ring_of_a_hexagonal_prism_serves() {
    assert_every_ring_exact("hexagonal prism", &hexagonal_prism(), 0.15);
}

/// The thin plate's 0.3-tall side faces: rings at d = 0.1 serve by the
/// ring surgery; at d = 0.2 the top and bottom bevels consume the face
/// (PK-1001 growth from the ring served at half size).
#[test]
fn thin_plate_side_rings_serve() {
    let plate = prism(&rect(0.0, 0.0, 4.0, 2.0), 0.0, 0.3);
    for d in [0.1, 0.2] {
        assert_every_ring_exact(&format!("plate d={d}"), &plate, d);
    }
}

#[test]
fn face_ring_scale_and_rotation_sweep() {
    let rotation =
        Transform3::rotation_around_axis(p(0.1, -0.2, 0.3), Vec3::new(1.0, 2.0, 3.0), 0.7)
            .expect("rotation");
    let base = hexagonal_prism();
    let rings = face_rings(&base);
    for s in [1e-2, 1.0, 1e2] {
        let place = rotation.compose(&Transform3::scale_uniform(s));
        let body = base.transform(&place).expect("placed prism");
        for (k, ring) in rings.iter().enumerate() {
            let placed: Vec<(Point3, Point3)> = ring
                .iter()
                .map(|(a, b)| (place.apply_point(a), place.apply_point(b)))
                .collect();
            let label = format!("scale {s}, face {k}");
            let out = chamfer(&body, &placed, 0.15 * s).unwrap_or_else(|e| panic!("{label}: {e}"));
            let exact = halfspace_oracle::chamfered_volume(&base, ring, 0.15);
            assert_exact(&label, &out, exact * s * s * s);
        }
    }
}
