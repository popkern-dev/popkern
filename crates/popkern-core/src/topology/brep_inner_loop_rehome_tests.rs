//! PK-1189: a hole the split's new edge crosses goes to the side holding
//! most of it, not the side its first edge's middle happens to lie on.

use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{NurbsCurve, Point3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle};

/// A closed degree-1 curve through `n` points of a circle.
fn polygon(points: impl Fn(f64) -> Point3, n: usize) -> NurbsCurve {
    let mut control_points: Vec<Point3> = (0..n)
        .map(|k| points(std::f64::consts::TAU * k as f64 / n as f64))
        .collect();
    control_points.push(control_points[0]);
    let m = control_points.len();
    let mut knots = vec![0.0];
    knots.extend((0..m).map(|k| k as f64 / (m - 1) as f64));
    knots.push(1.0);
    NurbsCurve {
        degree: 1,
        knots,
        weights: vec![1.0; m],
        control_points,
    }
}

#[test]
fn a_hole_crossing_the_split_goes_to_the_side_holding_most_of_it() {
    let tol = Tolerance::default();
    let mut b = BRep::unit_cube(tol).unwrap();
    // The top face, its plane chart.
    let top: FaceHandle = b
        .faces
        .keys()
        .find(|&f| {
            b.face_outer_loop_positions(f)
                .iter()
                .all(|p| (p.z - 1.0).abs() < 1e-12)
        })
        .unwrap();
    let SurfaceGeometry::Plane {
        origin,
        u_axis,
        v_axis,
        ..
    } = b.faces[top].surface.as_ref().unwrap().geometry.clone()
    else {
        panic!("top is not planar")
    };
    // A hole of radius 0.2 about (0.4, 0.5) in eight edges, split along
    // x = 0.56 (eight tenths of its radius out): the first edge's middle
    // (0.571, the old single probe) lies beyond the split, four of the
    // hole's 24 quarter stations do, the rest and its centre inside.
    let (cx, cy, r) = (0.4, 0.5, 0.2);
    let at = |th: f64| Point3::new(cx + r * th.cos(), cy + r * th.sin(), 1.0);
    let curve_3d = polygon(at, 64);
    let curve_uv = polygon(
        |th| {
            let d = at(th) - origin;
            Point3::new(d.dot(&u_axis), d.dot(&v_axis), 0.0)
        },
        64,
    );
    let hole = b
        .add_inner_loop(top, &curve_3d, &curve_uv, 8, &tol)
        .unwrap()
        .loop_;
    let _ = b.recompute_face_pcurves(top);
    // Split the top along x = 0.56 between its edges at y = 0 and y = 1.
    let x_cut = cx + 0.8 * r;
    let mut cut = Vec::new();
    for y in [0.0, 1.0] {
        let edge = b
            .edges
            .iter()
            .find(|(_, e)| {
                let h = &b.half_edges[e.half_edge];
                let (p, q) = (
                    b.vertices[h.origin].position,
                    b.vertices[b.half_edges[h.next].origin].position,
                );
                (p.z - 1.0).abs() < 1e-12
                    && (q.z - 1.0).abs() < 1e-12
                    && (p.y - y).abs() < 1e-12
                    && (q.y - y).abs() < 1e-12
            })
            .map(|(e, _)| e)
            .unwrap();
        cut.push(
            b.split_edge_at_point(edge, Point3::new(x_cut, y, 1.0), &tol)
                .unwrap()
                .vertex,
        );
    }
    let from = |b: &BRep, v| {
        b.half_edges
            .iter()
            .find(|(_, h)| h.origin == v && b.loops.get(h.loop_).is_some_and(|l| l.face == top))
            .map(|(h, _)| h)
            .unwrap()
    };
    let (h0, h1) = (from(&b, cut[0]), from(&b, cut[1]));
    let split = b.mef(h0, h1).unwrap();
    b.rehome_inner_loops_after_split(top, split.face);
    // The side holding the hole's centre (x < 0.56) owns it.
    let owner = b.loops[hole].face;
    let owner_xs: Vec<f64> = b
        .face_outer_loop_positions(owner)
        .iter()
        .map(|p| p.x)
        .collect();
    assert!(
        owner_xs.iter().all(|&x| x <= x_cut + 1e-12),
        "the hole went to the side beyond the split (outer loop xs {owner_xs:?})"
    );
}
