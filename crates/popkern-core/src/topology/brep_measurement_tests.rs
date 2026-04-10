//! Tests for BRep geometric deviation measurement.

use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn edge_deviation_on_unit_cube_is_small() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    for (eh, _) in &brep.edges {
        let dev = brep.measure_edge_deviation(eh);
        // Cube edges are linear — pcurves should match exactly.
        assert!(
            dev < 1e-6,
            "Edge deviation {dev} exceeds threshold for cube edge"
        );
    }
}

#[test]
fn vertex_deviation_on_unit_cube_is_small() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    for (vh, _) in &brep.vertices {
        let dev = brep.measure_vertex_deviation(vh);
        // Cube vertices should sit exactly on their edge curves.
        assert!(
            dev < 1e-10,
            "Vertex deviation {dev} exceeds threshold for cube vertex"
        );
    }
}

#[test]
fn edge_deviation_returns_zero_for_invalid_handle() {
    let brep = BRep::new(Tolerance::default());
    // Use a default (invalid) handle.
    let eh = slotmap::KeyData::from_ffi(u64::MAX).into();
    assert_eq!(brep.measure_edge_deviation(eh), 0.0);
}

#[test]
fn vertex_deviation_returns_zero_for_invalid_handle() {
    let brep = BRep::new(Tolerance::default());
    let vh = slotmap::KeyData::from_ffi(u64::MAX).into();
    assert_eq!(brep.measure_vertex_deviation(vh), 0.0);
}
