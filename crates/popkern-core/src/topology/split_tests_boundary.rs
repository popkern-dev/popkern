//! Splitting a free (boundary) edge.

use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// PK-530: a sheet-boundary edge whose curve stores a DESCENDING
/// parameter range (authored against its half-edge, as face removal
/// leaves them) must split into two edges that both carry curves
/// ending on their vertices. The NURBS conversion once took the
/// descending range as given, the curve split failed, and the new
/// edge silently kept no curve.
#[test]
fn split_boundary_edge_with_descending_curve_range_keeps_both_curves() {
    use crate::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
    let mut brep = BRep::open_box(Tolerance::default()).unwrap();
    let edge = brep.free_edges().into_iter().next().unwrap();
    let he = brep.edges[edge].half_edge;
    let a = brep.vertices[brep.half_edges[he].origin].position;
    let b = brep.vertices[brep.half_edges[brep.half_edges[he].next].origin].position;
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: a,
            direction: b - a,
        },
        t_start: 1.0,
        t_end: 0.0,
    });
    let split = brep.split_edge(edge, 0.3).unwrap();
    for e in [split.edge_a, split.edge_b] {
        let data = &brep.edges[e];
        let curve = data.curve.as_ref().expect("both halves carry a curve");
        let h = &brep.half_edges[data.half_edge];
        let ends = [
            brep.vertices[h.origin].position,
            brep.vertices[brep.half_edges[h.next].origin].position,
        ];
        for t in [curve.t_start, curve.t_end] {
            let p = curve.geometry.point_at(t);
            assert!(
                ends.iter().any(|q| (*q - p).length() < 1e-9),
                "curve end {p:?} is not on the edge's vertices {ends:?}"
            );
        }
    }
}
