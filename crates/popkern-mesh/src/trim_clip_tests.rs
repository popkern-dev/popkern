use super::*;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};

/// Create a simple MeshVertex at given position with UV.
fn make_vertex(x: f64, y: f64, z: f64, u: f64, v: f64) -> MeshVertex {
    MeshVertex {
        position: Point3::new(x, y, z),
        normal: Vec3::new(0.0, 0.0, 1.0),
        uv: (u, v),
    }
}

/// Create a Z-up plane surface for testing.
fn test_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

/// Helper: create a unit-square outer trim loop.
fn unit_square_loop() -> TrimLoop {
    TrimLoop {
        points: vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        points_3d: Vec::new(),
        is_outer: true,
        from_pcurves: true,
    }
}

#[test]
fn clip_triangle_against_vertical_line() {
    let tri_uv = [(0.0, 0.0), (2.0, 0.0), (1.0, 2.0)];
    let hits = clip_triangle_against_segment(&tri_uv, (1.0, -1.0), (1.0, 3.0));
    assert!(!hits.is_empty(), "should find at least one intersection");
    let edge0_hits: Vec<_> = hits.iter().filter(|h| h.0 == 0).collect();
    assert_eq!(edge0_hits.len(), 1, "edge 0 should have 1 intersection");
    assert!(
        (edge0_hits[0].1 - 0.5).abs() < 1e-10,
        "edge 0 t should be 0.5"
    );
    assert!(
        (edge0_hits[0].2.0 - 1.0).abs() < 1e-10,
        "edge 0 x should be 1.0"
    );
}

#[test]
fn apply_trim_preserves_inside_triangles() {
    let loops = vec![unit_square_loop()];
    let surface = test_plane();

    // Triangle fully inside the unit square.
    let vertices = vec![
        make_vertex(0.3, 0.3, 0.0, 0.3, 0.3),
        make_vertex(0.5, 0.3, 0.0, 0.5, 0.3),
        make_vertex(0.4, 0.5, 0.0, 0.4, 0.5),
    ];
    let triangles = vec![[0, 1, 2]];

    let (out_verts, out_tris) = apply_trim(&vertices, &triangles, &loops, &surface, true);

    assert_eq!(out_tris.len(), 1, "inside triangle should be preserved");
    assert_eq!(out_tris[0], [0, 1, 2]);
    assert_eq!(out_verts.len(), 3);
}

#[test]
fn apply_trim_removes_outside_triangles() {
    let loops = vec![unit_square_loop()];
    let surface = test_plane();

    // Triangle fully outside the unit square.
    let vertices = vec![
        make_vertex(1.5, 1.5, 0.0, 1.5, 1.5),
        make_vertex(2.0, 1.5, 0.0, 2.0, 1.5),
        make_vertex(1.7, 2.0, 0.0, 1.7, 2.0),
    ];
    let triangles = vec![[0, 1, 2]];

    let (_out_verts, out_tris) = apply_trim(&vertices, &triangles, &loops, &surface, true);

    assert!(out_tris.is_empty(), "outside triangle should be removed");
}

#[test]
fn segment_segment_intersection_basic() {
    // Crossing segments.
    let result = segment_segment_intersection((0.0, 0.0), (2.0, 2.0), (0.0, 2.0), (2.0, 0.0));
    assert!(result.is_some(), "crossing segments should intersect");
    let (t, pt) = result.unwrap();
    assert!((t - 0.5).abs() < 1e-10);
    assert!((pt.0 - 1.0).abs() < 1e-10);
    assert!((pt.1 - 1.0).abs() < 1e-10);
}

#[test]
fn segment_segment_intersection_parallel() {
    // Parallel segments.
    let result = segment_segment_intersection((0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0));
    assert!(result.is_none(), "parallel segments should not intersect");
}

#[test]
fn segment_segment_intersection_non_overlapping() {
    // Non-overlapping segments.
    let result = segment_segment_intersection((0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 1.0));
    assert!(
        result.is_none(),
        "non-overlapping segments should not intersect"
    );
}

#[test]
fn split_triangle_at_edge_vertices_no_splits() {
    let tri = [0, 1, 2];
    let splits = vec![(0, vec![]), (1, vec![]), (2, vec![])];
    let mut output = Vec::new();
    split_triangle_at_edge_vertices(&tri, &splits, &mut output);
    assert_eq!(output.len(), 1);
    assert_eq!(output[0], [0, 1, 2]);
}

#[test]
fn split_triangle_at_edge_vertices_with_split() {
    // Triangle [0,1,2] with vertex 3 inserted on edge 0 (between 0 and 1).
    let tri = [0, 1, 2];
    let splits = vec![(0, vec![3]), (1, vec![]), (2, vec![])];
    let mut output = Vec::new();
    split_triangle_at_edge_vertices(&tri, &splits, &mut output);
    // Ring: [0, 3, 1, 2] -> fan from 0: [0,3,1], [0,1,2]
    assert_eq!(output.len(), 2);
    assert_eq!(output[0], [0, 3, 1]);
    assert_eq!(output[1], [0, 1, 2]);
}

#[test]
fn apply_trim_works_without_pcurves() {
    // Simulate a projection-based (non-pcurve) trim loop.
    // The clipping pipeline should still produce reasonable output.
    let trim_loop = TrimLoop {
        points: vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        points_3d: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        is_outer: true,
        from_pcurves: false,
    };
    let loops = vec![trim_loop];
    let surface = test_plane();

    // A triangle that straddles the trim boundary (extends beyond x=1).
    let vertices = vec![
        make_vertex(0.5, 0.3, 0.0, 0.5, 0.3),
        make_vertex(1.5, 0.3, 0.0, 1.5, 0.3),
        make_vertex(0.8, 0.7, 0.0, 0.8, 0.7),
    ];
    let triangles = vec![[0, 1, 2]];

    let (out_verts, out_tris) = apply_trim(&vertices, &triangles, &loops, &surface, true);

    // The clipping pipeline should produce at least one triangle from the
    // part of the original triangle that lies inside the boundary.
    assert!(
        !out_tris.is_empty(),
        "clipping should produce at least one triangle for the inside region"
    );

    // Verify no vertices have extreme coordinates (no wild geometry).
    for v in &out_verts {
        assert!(
            v.position.x.abs() < 10.0 && v.position.y.abs() < 10.0 && v.position.z.abs() < 10.0,
            "vertex position should not be extreme: {:?}",
            v.position
        );
    }
}
