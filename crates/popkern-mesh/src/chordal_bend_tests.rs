use super::*;

const T: [Uv; 3] = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];

fn mids() -> [Option<Uv>; 3] {
    [None; 3]
}

#[test]
fn a_straight_net_lies_in_its_triangle() {
    assert!(net_inside(T, mids(), [None; 3]));
}

#[test]
fn a_bowed_edge_needs_the_triangle_across_it() {
    // Edge 0 (the u-axis) bowed down by 0.1: outside T.
    let mut c = mids();
    c[0] = Some((0.5, -0.2));
    assert!(!net_inside(T, c, [None; 3]));
    let below = [(0.0, 0.0), (1.0, -1.0), (1.0, 0.0)];
    assert!(net_inside(T, c, [Some(below), None, None]));
    // Bowed into T, nothing across is needed.
    c[0] = Some((0.5, 0.2));
    assert!(net_inside(T, c, [None; 3]));
}

#[test]
fn a_control_past_the_corner_wedge_is_refused() {
    // Slid along its edge beyond the vertex, the control makes the
    // hull reach round the corner, out of the triangles it may use.
    let mut c = mids();
    c[0] = Some((1.3, -0.05));
    let below = [(0.0, 0.0), (2.0, -1.0), (2.0, 0.0)];
    assert!(!net_inside(T, c, [Some(below), None, None]));
}

#[test]
fn a_bend_pass_keeps_edges_straight_beside_unbent_triangles() {
    // Two triangles of a flat plane: nothing to gain, nothing changes.
    let s = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: popkern_core::geometry::Vec3::Z,
        u_axis: popkern_core::geometry::Vec3::new(1.0, 0.0, 0.0),
        v_axis: popkern_core::geometry::Vec3::new(0.0, 1.0, 0.0),
    };
    let at = |q: Uv| Point3::new(q.0, q.1, 0.0);
    let a = BendTriangle {
        weld: [0, 1, 2],
        w: T,
        p: T.map(at),
    };
    let b = BendTriangle {
        weld: [1, 3, 2],
        w: [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        p: [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)].map(at),
    };
    let mut bounds = vec![1.0, 1.0];
    let mut calls = Vec::new();
    bend_pass(&s, &[Some(a), Some(b)], &mut bounds, 0.5, |w, _, c| {
        calls.push((w, c));
        Some(2.0)
    });
    // A bend that does not improve is rejected and the bounds stand.
    assert_eq!(bounds, vec![1.0, 1.0]);
    assert!(!calls.is_empty());
}
