//! PK-1189: a closed edge splits at the point asked for, never its mirror
//! image. Its curve's ends are one point, so the curve's direction along
//! the half-edge is read as `split_edge` reads it.

use crate::geometry::Transform3;
use crate::geometry::curve::CurveEval;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn a_closed_edge_splits_where_asked() {
    let tol = Tolerance::default();
    let upright = BRep::cylinder(1.0, 2.0, tol).unwrap();
    let flipped = upright
        .transform(&Transform3::rotation_x(std::f64::consts::PI))
        .unwrap();
    for body in [upright, flipped] {
        let closed: Vec<_> = body
            .edges
            .iter()
            .filter(|(_, e)| {
                let h = &body.half_edges[e.half_edge];
                h.origin == body.half_edges[h.next].origin
            })
            .map(|(e, _)| e)
            .collect();
        assert!(!closed.is_empty());
        for edge in closed {
            for (fraction, against) in
                [(0.125, false), (0.875, false), (0.125, true), (0.875, true)]
            {
                let mut b = body.clone();
                if against {
                    // The same circle run the other way: its curve against
                    // the primary half-edge, its pcurves along it.
                    let c = b.edges[edge].curve.clone().unwrap();
                    let n = c.geometry.to_nurbs_in_range((c.t_start, c.t_end));
                    let (a, z) = n.parameter_range();
                    let mut r = n.clone();
                    r.control_points.reverse();
                    r.weights.reverse();
                    r.knots = n.knots.iter().rev().map(|k| a + z - k).collect();
                    let (ws, we) = (a, z);
                    b.edges[edge].curve = Some(crate::geometry::curve::CurveOnEdge {
                        geometry: crate::geometry::curve::CurveGeometry::Nurbs(r),
                        t_start: ws,
                        t_end: we,
                    });
                }
                let c = b.edges[edge].curve.clone().unwrap();
                let point = c
                    .geometry
                    .point_at(c.t_start + (c.t_end - c.t_start) * fraction);
                let split = b.split_edge_at_point(edge, point, &tol).unwrap();
                let landed = b.vertices[split.vertex].position.distance_to(&point);
                assert!(
                    landed < 1e-9,
                    "split at fraction {fraction} (against {against}) landed {landed:.3e} away"
                );
            }
        }
    }
}
