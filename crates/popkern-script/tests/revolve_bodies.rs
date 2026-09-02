//! PK-408: revolves recorded through the recorder serve VALID bodies —
//! a welded ring, a torus from a circle profile, a disc whose profile
//! touches the axis, a capped partial revolve — and a boolean on a
//! revolved ring serves. Before, every revolve replayed "served" and
//! failed Part validation (unwelded faces, inconsistent orientation).

use popkern_history::feature::AxisDesc;
use popkern_script::prelude::*;
use popkern_script::{Plane, RevolveEnd};

fn z_axis_at(x: f64) -> AxisDesc {
    AxisDesc::Literal {
        origin: [x, 0.0, 0.0],
        direction: [0.0, 0.0, 1.0],
    }
}

fn counts(m: &mut Model) -> (usize, usize, usize) {
    let brep = m.to_brep().expect("brep");
    (brep.vertices.len(), brep.edges.len(), brep.faces.len())
}

#[test]
fn full_revolution_of_a_rectangle_is_a_valid_ring() {
    let mut m = Model::new();
    m.profile(Plane::XZ)
        .rect(0.2, 0.1)
        .at(0.6, 3.0)
        .sketch("ring_section")
        .expect("sketch")
        .revolve(&mut m, z_axis_at(0.0), RevolveEnd::FullCircle, "ring")
        .expect("revolve");
    assert_eq!(counts(&mut m), (4, 8, 4));
    let part = m.to_part().expect("ring validates");
    assert_eq!(part.face_count(), 4);
    // Tessellation welds between the NURBS revolution faces are PK-409
    // (measured 789 open mesh edges on this ring); the body is valid.
}

#[test]
fn circle_profile_revolves_to_a_torus_with_its_curvature() {
    let mut m = Model::new();
    m.profile(Plane::XY)
        .circle(0.3)
        .at(1.0, 0.0)
        .sketch("tube_section")
        .expect("sketch")
        .revolve(
            &mut m,
            AxisDesc::Literal {
                origin: [0.0, 0.0, 0.0],
                direction: [0.0, 1.0, 0.0],
            },
            RevolveEnd::FullCircle,
            "torus",
        )
        .expect("revolve");
    let brep = m.to_brep().expect("brep");
    // One periodic face, one seam circle + one profile circle, one vertex.
    assert_eq!(
        (brep.vertices.len(), brep.edges.len(), brep.faces.len()),
        (1, 2, 1)
    );
    let radii: Vec<f64> = brep
        .vertices
        .values()
        .map(|v| (v.position.x.powi(2) + v.position.z.powi(2)).sqrt())
        .collect();
    assert!(
        radii
            .iter()
            .all(|r| (r - 1.3).abs() < 1e-9 || (r - 0.7).abs() < 1e-9),
        "{radii:?}"
    );
    m.to_part().expect("torus validates");
    m.assert_watertight().expect("watertight");
}

#[test]
fn profile_touching_the_axis_revolves_to_a_disc() {
    let mut m = Model::new();
    m.profile(Plane::XZ)
        .rect(0.2, 0.1)
        .at(0.1, 3.0)
        .sketch("disc_section")
        .expect("sketch")
        .revolve(&mut m, z_axis_at(0.0), RevolveEnd::FullCircle, "disc")
        .expect("revolve");
    assert_eq!(counts(&mut m), (4, 5, 3));
    m.to_part().expect("disc validates");
    // Tessellation of the pole faces is PK-409 (measured 99 open mesh
    // edges); the body itself is valid.
}

#[test]
fn partial_revolution_is_capped_and_valid() {
    let mut m = Model::new();
    m.profile(Plane::XZ)
        .rect(0.2, 0.1)
        .at(0.6, 3.0)
        .sketch("arc_section")
        .expect("sketch")
        .revolve(
            &mut m,
            z_axis_at(0.0),
            RevolveEnd::Blind(90.0.into()),
            "quarter",
        )
        .expect("revolve");
    assert_eq!(counts(&mut m), (8, 12, 6));
    m.to_part().expect("quarter validates");
}

#[test]
fn a_boolean_on_a_revolved_ring_serves() {
    let mut m = Model::new();
    let ring = m
        .profile(Plane::XZ)
        .rect(0.2, 0.1)
        .at(0.6, 3.0)
        .sketch("ring_section")
        .expect("sketch")
        .revolve(&mut m, z_axis_at(0.0), RevolveEnd::FullCircle, "ring")
        .expect("revolve");
    let cutter = m
        .tool(|m| {
            m.profile(Plane::xy_at(2.9))
                .circle(0.05)
                .at(0.6, 0.0)
                .sketch("cutter_sketch")?
                .extrude(m, End::Blind(0.3.into()))
        })
        .expect("cutter");
    ring.subtract(&mut m, cutter).expect("subtract serves");
    let part = m.to_part().expect("bored ring validates");
    assert!(
        part.face_count() > 4,
        "the bore adds faces: {}",
        part.face_count()
    );
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("revolve_bodies");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join("bored_ring.pkn");
    m.save(&path).expect("save");
    let mut again = Model::open(&path).expect("open");
    assert_eq!(again.replay().expect("replay").failures().count(), 0);
    assert_eq!(
        again.to_part().expect("part").face_count(),
        part.face_count()
    );
}
