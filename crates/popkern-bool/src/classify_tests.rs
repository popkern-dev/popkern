// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
#[cfg(test)]
use super::*;
use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Translate all vertices and geometry of a BRep by a vector offset.
fn translate(mut brep: BRep, offset: Vec3) -> BRep {
    use popkern_core::geometry::{curve::CurveGeometry, surface::SurfaceGeometry};
    for (_, vertex) in &mut brep.vertices {
        vertex.position = vertex.position + offset;
    }
    for (_, face) in &mut brep.faces {
        if let Some(ref mut sof) = face.surface {
            if let SurfaceGeometry::Plane { ref mut origin, .. } = sof.geometry {
                *origin = *origin + offset;
            }
        }
    }
    for (_, edge) in &mut brep.edges {
        if let Some(ref mut coe) = edge.curve {
            if let CurveGeometry::Line { ref mut origin, .. } = coe.geometry {
                *origin = *origin + offset;
            }
        }
    }
    brep
}

// -- ray_intersects_aabb

#[test]
fn ray_through_aabb_center_hits() {
    let aabb =
        Aabb::from_points(&[Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)]).expect("aabb");
    let origin = Point3::new(-1.0, 0.5, 0.5);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    assert!(ray_intersects_aabb(&origin, &dir, &aabb));
}

#[test]
fn ray_missing_aabb_returns_false() {
    let aabb =
        Aabb::from_points(&[Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)]).expect("aabb");
    let origin = Point3::new(-1.0, 5.0, 0.5);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    assert!(!ray_intersects_aabb(&origin, &dir, &aabb));
}

#[test]
fn ray_behind_aabb_returns_false() {
    let aabb =
        Aabb::from_points(&[Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)]).expect("aabb");
    // Origin beyond the box, ray pointing away.
    let origin = Point3::new(2.0, 0.5, 0.5);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    assert!(!ray_intersects_aabb(&origin, &dir, &aabb));
}

#[test]
fn ray_from_inside_aabb_hits() {
    let aabb =
        Aabb::from_points(&[Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)]).expect("aabb");
    let origin = Point3::new(0.5, 0.5, 0.5);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    assert!(ray_intersects_aabb(&origin, &dir, &aabb));
}

#[test]
fn ray_parallel_to_slab_outside_misses() {
    let aabb =
        Aabb::from_points(&[Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)]).expect("aabb");
    // Origin outside y-slab, direction has zero y component.
    let origin = Point3::new(0.5, 2.0, 0.5);
    let dir = Vec3::new(1.0, 0.0, 0.0);
    assert!(!ray_intersects_aabb(&origin, &dir, &aabb));
}

// -- point_in_solid

fn unit_cube_with_aabbs() -> (BRep, Vec<(FaceHandle, Aabb)>) {
    let cube = BRep::unit_cube(tol()).expect("unit cube");
    let aabbs = crate::intersect::compute_face_aabbs(&cube);
    (cube, aabbs)
}

#[test]
fn point_inside_unit_cube() {
    let (cube, aabbs) = unit_cube_with_aabbs();
    let pt = Point3::new(0.5, 0.5, 0.5);
    let result = point_in_solid(&pt, &cube, &aabbs, &tol()).expect("classification should succeed");
    assert_eq!(result, PointClassification::Inside);
}

#[test]
fn point_outside_unit_cube() {
    let (cube, aabbs) = unit_cube_with_aabbs();
    let pt = Point3::new(2.0, 0.0, 0.0);
    let result = point_in_solid(&pt, &cube, &aabbs, &tol()).expect("classification should succeed");
    assert_eq!(result, PointClassification::Outside);
}

#[test]
fn point_far_outside_unit_cube() {
    let (cube, aabbs) = unit_cube_with_aabbs();
    let pt = Point3::new(100.0, 100.0, 100.0);
    let result = point_in_solid(&pt, &cube, &aabbs, &tol()).expect("classification should succeed");
    assert_eq!(result, PointClassification::Outside);
}

#[test]
fn point_negative_outside_unit_cube() {
    let (cube, aabbs) = unit_cube_with_aabbs();
    let pt = Point3::new(-1.0, -1.0, -1.0);
    let result = point_in_solid(&pt, &cube, &aabbs, &tol()).expect("classification should succeed");
    assert_eq!(result, PointClassification::Outside);
}

// -- classify_fragments

#[test]
fn classify_cube_faces_against_disjoint_cube() {
    let cube_a = BRep::unit_cube(tol()).expect("cube a");
    let cube_b = translate(
        BRep::unit_cube(tol()).expect("cube b"),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let aabbs_b = crate::intersect::compute_face_aabbs(&cube_b);
    let empty_edges = HashSet::new();

    let result = classify_fragments(&cube_a, &cube_b, &aabbs_b, &empty_edges, &tol())
        .expect("classification should succeed");

    // All faces of cube_a should be outside cube_b.
    assert_eq!(result.len(), 6);
    for (_, class) in &result {
        assert_eq!(*class, Classification::Outside);
    }
}

#[test]
fn classify_small_cube_inside_large_cube() {
    // Small cube at (0.25..0.75) inside unit cube (0..1).
    let outer = BRep::unit_cube(tol()).expect("outer cube");
    let inner = translate(
        {
            let mut c = BRep::unit_cube(tol()).expect("inner cube");
            // Scale to 0.5 size: vertices go from 0..1 to 0..0.5
            for (_, v) in &mut c.vertices {
                v.position =
                    Point3::new(v.position.x * 0.5, v.position.y * 0.5, v.position.z * 0.5);
            }
            c
        },
        Vec3::new(0.25, 0.25, 0.25),
    );

    let aabbs_outer = crate::intersect::compute_face_aabbs(&outer);
    let empty_edges = HashSet::new();

    let result = classify_fragments(&inner, &outer, &aabbs_outer, &empty_edges, &tol())
        .expect("classification should succeed");

    // All faces of the inner cube should be inside the outer cube.
    assert_eq!(result.len(), 6);
    for (_, class) in &result {
        assert_eq!(
            *class,
            Classification::Inside,
            "inner face should be inside outer cube"
        );
    }
}

#[test]
fn classify_fragments_empty_brep() {
    let empty = BRep::new(tol());
    let cube = BRep::unit_cube(tol()).expect("cube");
    let aabbs = crate::intersect::compute_face_aabbs(&cube);
    let empty_edges = HashSet::new();

    let result = classify_fragments(&empty, &cube, &aabbs, &empty_edges, &tol())
        .expect("classification should succeed");

    assert!(result.is_empty());
}

// -- pick_sample_point

#[test]
fn pick_sample_point_on_cube_face_is_inside_bbox() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    for fh in cube.faces.keys() {
        let pt = pick_sample_point(&cube, fh);
        assert!(pt.x >= -0.1 && pt.x <= 1.1, "x out of range: {}", pt.x);
        assert!(pt.y >= -0.1 && pt.y <= 1.1, "y out of range: {}", pt.y);
        assert!(pt.z >= -0.1 && pt.z <= 1.1, "z out of range: {}", pt.z);
    }
}

// -- face_parameter_range

#[test]
fn face_parameter_range_is_finite_for_cube() {
    let cube = BRep::unit_cube(tol()).expect("cube");
    for fh in cube.faces.keys() {
        let ((u_min, u_max), (v_min, v_max)) = face_parameter_range(&cube, fh);
        assert!(u_min.is_finite(), "u_min not finite");
        assert!(u_max.is_finite(), "u_max not finite");
        assert!(v_min.is_finite(), "v_min not finite");
        assert!(v_max.is_finite(), "v_max not finite");
        assert!(u_min < u_max, "u_min >= u_max");
        assert!(v_min < v_max, "v_min >= v_max");
    }
}

#[test]
fn face_parameter_range_face_with_no_surface_returns_default() {
    use popkern_core::topology::{Face, HalfEdgeHandle, Loop};
    let mut brep = BRep::new(tol());
    // Create a face with no surface.
    let loop_h = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: popkern_core::topology::FaceHandle::default(),
    });
    let face_h = brep.faces.insert(Face {
        outer_loop: loop_h,
        inner_loops: Vec::new(),
        shell: popkern_core::topology::ShellHandle::default(),
        surface: None,
        name: None,
        tolerance: None,
    });
    brep.loops[loop_h].face = face_h;
    let ((u_min, u_max), (v_min, v_max)) = face_parameter_range(&brep, face_h);
    // Should return the fallback (-10, 10) x (-10, 10).
    assert_eq!(u_min, -10.0);
    assert_eq!(u_max, 10.0);
    assert_eq!(v_min, -10.0);
    assert_eq!(v_max, 10.0);
}

// Edge-case / fallback tests live in classify_edge_tests.rs.
