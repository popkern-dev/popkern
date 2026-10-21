use super::*;

fn plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: popkern_core::geometry::Vec3::Z,
        u_axis: popkern_core::geometry::Vec3::new(1.0, 0.0, 0.0),
        v_axis: popkern_core::geometry::Vec3::new(0.0, 1.0, 0.0),
    }
}

fn tri(weld: [usize; 3], w: [Uv; 3]) -> BendTriangle {
    BendTriangle {
        weld,
        w,
        p: w.map(|q| Point3::new(q.0, q.1, 0.0)),
    }
}

/// The unit square as two triangles sharing the diagonal 1–2.
fn square() -> [Option<BendTriangle>; 2] {
    [
        Some(tri([0, 1, 2], [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)])),
        Some(tri([1, 3, 2], [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])),
    ]
}

#[test]
fn a_shared_edge_subdivides_to_the_same_floats_from_both_sides() {
    let [Some(a), Some(b)] = square() else {
        return;
    };
    // Edge 1–2 is a's edge 1 (third coordinate 0) and b's edge 2.
    for s in 0..=8 {
        let from_a = sub_vertex(&a, [0, 8 - s, s], 8);
        let from_b = sub_vertex(&b, [8 - s, 0, s], 8);
        let (ya, yb) = (from_a.0, from_b.0);
        assert!(ya.0.to_bits() == yb.0.to_bits() && ya.1.to_bits() == yb.1.to_bits());
    }
}

#[test]
fn a_region_is_refined_into_oriented_pieces_and_improves() {
    let tris = square();
    let mut bounds = vec![1.0, 1.0];
    let mut pieces = 0;
    refine_pass(
        &plane(),
        &tris,
        &mut bounds,
        &HashSet::new(),
        0.5,
        |w, _| {
            assert!(orient(w[0], w[1], w[2]) > 0.0);
            pieces += 1;
            Some(0.25)
        },
    );
    // The first level (2: four pieces a triangle) meets the tolerance.
    assert_eq!(bounds, vec![0.25, 0.25]);
    assert_eq!(pieces, 8);
}

#[test]
fn a_region_that_does_not_improve_keeps_its_bounds() {
    let tris = square();
    let mut bounds = vec![1.0, 1.0];
    refine_pass(
        &plane(),
        &tris,
        &mut bounds,
        &HashSet::new(),
        0.5,
        |_, _| Some(2.0),
    );
    assert_eq!(bounds, vec![1.0, 1.0]);
}

#[test]
fn a_triangle_bent_toward_a_neighbour_outside_stays_out() {
    // Only the first is over; the diagonal it shares with the second
    // (under, so outside the region) was left bent: refining it would
    // straighten the second's bent edge.
    let tris = square();
    let mut bounds = vec![1.0, 0.1];
    let bent: HashSet<EdgeKey> = [(1, 2)].into_iter().collect();
    let mut calls = 0;
    refine_pass(&plane(), &tris, &mut bounds, &bent, 0.5, |_, _| {
        calls += 1;
        Some(0.0)
    });
    assert_eq!((bounds, calls), (vec![1.0, 0.1], 0));
}

#[test]
fn a_foot_on_the_region_boundary_is_pulled_strictly_inside() {
    // A foot on the boundary would flatten its piece (a fold); halfway
    // back it lies inside, and a foot outside every triangle falls back
    // to the sub-vertex itself only when no step toward it is inside.
    let shapes = [[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]];
    let y = (0.25, 0.25);
    let q = inside(y, (0.0, 0.25), &shapes);
    assert!(q.0 > 0.0 && q.0 < 0.25 && q.1 == 0.25, "{q:?}");
    assert_eq!(inside(y, (0.2, 0.2), &shapes), (0.2, 0.2));
    let edge = (0.5, 0.0);
    assert_eq!(inside(edge, (0.5, -0.1), &shapes), edge);
}
