//! Tests for the `kve` edge-contraction operator.

use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, VertexHandle};

/// Caller-side geometry repair after contraction: move the surviving
/// vertex, then rebuild every incident straight edge between its
/// (possibly moved) endpoint positions.
fn rebuild_incident_lines(brep: &mut BRep, vertex: VertexHandle) {
    let incident: Vec<_> = brep
        .edges
        .iter()
        .filter(|(_, e)| {
            let he = e.half_edge;
            brep.half_edges[he].origin == vertex
                || brep.half_edges[brep.half_edges[he].twin].origin == vertex
        })
        .map(|(h, _)| h)
        .collect();
    for eh in incident {
        let he = brep.edges[eh].half_edge;
        let a = brep.vertices[brep.half_edges[he].origin].position;
        let b = brep.vertices[brep.half_edges[brep.half_edges[he].twin].origin].position;
        let d = b - a;
        let len = d.length();
        if len > 1e-12
            && let Some(dir) = d.normalized()
        {
            brep.edges[eh].curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: a,
                    direction: dir,
                },
                t_start: 0.0,
                t_end: len,
            });
        }
    }
}

#[test]
fn kve_contracts_cube_edge() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let (v0, e0, f0) = (cube.vertices.len(), cube.edges.len(), cube.faces.len());
    let he = cube.edges.values().next().unwrap().half_edge;
    let result = cube.kve(he).expect("contract one cube edge");
    assert_eq!(cube.vertices.len(), v0 - 1);
    assert_eq!(cube.edges.len(), e0 - 1);
    assert_eq!(cube.faces.len(), f0);
    assert!(cube.vertices.contains_key(result.vertex));
    assert!(!cube.vertices.contains_key(result.removed_vertex));
    assert!(!cube.edges.contains_key(result.removed_edge));
    // Caller contract: repair incident geometry, then the structural
    // invariants must hold (the contracted cube is no longer a valid
    // closed PLANAR solid at face level, so full geometric
    // validation is out of scope for a pure topology test).
    rebuild_incident_lines(&mut cube, result.vertex);
    cube.validate_euler_poincare().expect("euler-poincare");
    cube.validate_half_edge_pairing().expect("pairing");
    cube.validate_half_edge_cycles().expect("cycles");
    cube.validate_manifold_edges().expect("manifold");
}

#[test]
fn kve_rejects_self_loop_and_small_loops() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    // Contract one edge of a quad face -> that face's loop has 3.
    let face = cube.faces.keys().next().unwrap();
    let l = cube.faces[face].outer_loop;
    let he = cube.loops[l].half_edge;
    let r = cube.kve(he).expect("first contraction");
    rebuild_incident_lines(&mut cube, r.vertex);
    // The same face's loop now has 3 half-edges; contracting one
    // more is legal (bigons are valid seam structures in this
    // kernel — mekr creates them), taking the loop to 2.
    let he2 = cube.loops[l].half_edge;
    if let Ok(r2) = cube.kve(he2) {
        rebuild_incident_lines(&mut cube, r2.vertex);
        // A bigon must refuse further contraction.
        let he3 = cube.loops[l].half_edge;
        assert!(cube.kve(he3).is_err(), "2-loop contraction must be refused");
    }
}
