//! Spur-collapse pass tests: an authored MEV spur is retired; real
//! topology is untouched (PK-293).
use super::*;
use popkern_core::geometry::Point3;

#[test]
fn mev_spur_is_collapsed_and_cube_survives() {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).expect("cube");
    let edges_before = brep.edges.len();
    let vertices_before = brep.vertices.len();
    // Author a short spur off an existing vertex with mev.
    let (vh, _) = brep.vertices.iter().next().expect("vertex");
    let he_out = brep.vertices[vh].outgoing;
    let base = brep.vertices[vh].position;
    let tip = Point3::new(base.x + 1.0e-7, base.y, base.z);
    brep.mev(he_out, tip).expect("mev spur");
    assert_eq!(brep.edges.len(), edges_before + 1);

    let report = collapse_degenerate_spur_edges(&mut brep, &tol);
    assert_eq!(report.spurs_removed, 1, "the authored spur must retire");
    assert_eq!(brep.edges.len(), edges_before);
    assert_eq!(brep.vertices.len(), vertices_before);
    brep.validate().expect("cube must stay valid");
}

#[test]
fn long_spur_is_left_visible() {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).expect("cube");
    let (vh, _) = brep.vertices.iter().next().expect("vertex");
    let he_out = brep.vertices[vh].outgoing;
    let base = brep.vertices[vh].position;
    let tip = Point3::new(base.x + 0.25, base.y, base.z);
    brep.mev(he_out, tip).expect("mev spur");
    let report = collapse_degenerate_spur_edges(&mut brep, &tol);
    assert_eq!(
        report.spurs_removed, 0,
        "a macroscopic spur is an upstream bug and stays visible"
    );
}
