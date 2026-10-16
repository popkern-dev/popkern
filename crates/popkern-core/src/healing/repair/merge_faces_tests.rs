use super::*;
use crate::constructors::builder::BRepBuilder;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

// -----------------------------------------------------------------------
// Original tests (migrated from inline module)
// -----------------------------------------------------------------------

#[test]
fn unit_cube_has_no_coplanar_pairs() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let pairs = detect_coplanar_pairs(&brep, &tol);
    assert!(
        pairs.is_empty(),
        "all cube faces are perpendicular, expected 0 coplanar pairs, got {}",
        pairs.len()
    );
}

#[test]
fn merge_coplanar_faces_on_cube_detects_zero() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let result = merge_coplanar_faces(&mut brep, &tol);
    assert_eq!(result.pairs_detected, 0);
    assert_eq!(result.pairs_merged, 0);
}

#[test]
fn empty_brep_has_no_coplanar_pairs() {
    let brep = BRep::new(Tolerance::default());
    let tol = Tolerance::default();
    let pairs = detect_coplanar_pairs(&brep, &tol);
    assert!(pairs.is_empty());
}

#[test]
fn canonical_pair_is_stable() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let faces: Vec<FaceHandle> = brep.faces.keys().collect();
    if faces.len() >= 2 {
        let a = faces[0];
        let b = faces[1];
        assert_eq!(canonical_pair(a, b), canonical_pair(b, a));
    }
}

#[test]
fn surface_types_compatible_same() {
    let plane_a = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let plane_b = SurfaceGeometry::Plane {
        origin: Point3::new(1.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert!(surface_types_compatible(&plane_a, &plane_b));
}

#[test]
fn surface_types_incompatible_cross() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    assert!(!surface_types_compatible(&plane, &cyl));
}

// -----------------------------------------------------------------------
// Helpers: build BRep with coplanar adjacent faces
// -----------------------------------------------------------------------

/// Build a BRep with two coplanar adjacent faces sharing one edge.
/// Both faces lie in the z=0 plane with the same normal.
fn make_coplanar_pair() -> BRep {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let v4 = b.add_vertex(Point3::new(2.0, 0.0, 0.0));
    let v5 = b.add_vertex(Point3::new(2.0, 1.0, 0.0));

    let solid = b.add_solid();
    let shell = b.add_shell(solid);

    // Face A: v0 -> v1 -> v2 -> v3
    let (_, he01, _) = b.add_edge(v0, v1);
    let (_, he12, he21) = b.add_edge(v1, v2); // shared edge
    let (_, he23, _) = b.add_edge(v2, v3);
    let (_, he30, _) = b.add_edge(v3, v0);

    // Face B: v1 -> v4 -> v5 -> v2 (uses twin of shared edge)
    let (_, he14, _) = b.add_edge(v1, v4);
    let (_, he45, _) = b.add_edge(v4, v5);
    let (_, he52, _) = b.add_edge(v5, v2);

    let fa = b.add_face(shell);
    b.add_loop(&[he01, he12, he23, he30], fa).unwrap();

    let fb = b.add_face(shell);
    b.add_loop(&[he21, he14, he45, he52], fb).unwrap();

    // Both faces get the same z=0 plane surface
    let plane = SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        orientation: true,
    };
    b.set_face_surface(fa, plane.clone());
    b.set_face_surface(fb, plane);

    b.build_unchecked()
}

/// Build a BRep with two faces on parallel but distinct planes (not coplanar).
fn make_parallel_but_offset_pair() -> BRep {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.5));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let v4 = b.add_vertex(Point3::new(2.0, 0.0, 0.0));
    let v5 = b.add_vertex(Point3::new(2.0, 1.0, 0.5));

    let solid = b.add_solid();
    let shell = b.add_shell(solid);

    let (_, he01, _) = b.add_edge(v0, v1);
    let (_, he12, he21) = b.add_edge(v1, v2);
    let (_, he23, _) = b.add_edge(v2, v3);
    let (_, he30, _) = b.add_edge(v3, v0);
    let (_, he14, _) = b.add_edge(v1, v4);
    let (_, he45, _) = b.add_edge(v4, v5);
    let (_, he52, _) = b.add_edge(v5, v2);

    let fa = b.add_face(shell);
    b.add_loop(&[he01, he12, he23, he30], fa).unwrap();
    let fb = b.add_face(shell);
    b.add_loop(&[he21, he14, he45, he52], fb).unwrap();

    let plane_a = SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        orientation: true,
    };
    let plane_b = SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 1.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        orientation: true,
    };
    b.set_face_surface(fa, plane_a);
    b.set_face_surface(fb, plane_b);

    b.build_unchecked()
}

// -----------------------------------------------------------------------
// Tests with coplanar face geometry
// -----------------------------------------------------------------------

#[test]
fn detect_coplanar_pair_on_same_plane() {
    let brep = make_coplanar_pair();
    let tol = Tolerance::default();
    let pairs = detect_coplanar_pairs(&brep, &tol);
    assert_eq!(
        pairs.len(),
        1,
        "two adjacent faces on the same plane should be detected as coplanar"
    );
}

#[test]
fn merge_coplanar_pair() {
    let mut brep = make_coplanar_pair();
    let initial_faces = brep.faces.len();
    let tol = Tolerance::default();
    let result = merge_coplanar_faces(&mut brep, &tol);
    assert_eq!(result.pairs_detected, 1);
    assert_eq!(result.pairs_merged, 1);
    assert_eq!(
        brep.faces.len(),
        initial_faces - 1,
        "merge should remove one face"
    );
}

#[test]
fn parallel_but_offset_planes_not_coplanar() {
    let brep = make_parallel_but_offset_pair();
    let tol = Tolerance::default();
    let pairs = detect_coplanar_pairs(&brep, &tol);
    assert!(
        pairs.is_empty(),
        "parallel planes at different heights should not be coplanar"
    );
}

#[test]
fn surfaces_coplanar_same_plane_returns_true() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    assert!(surfaces_coplanar(&plane, &plane, 1e-6));
}

#[test]
fn surfaces_coplanar_anti_parallel_normals() {
    let plane_a = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let plane_b = SurfaceGeometry::Plane {
        origin: Point3::new(1.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, -1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, -1.0, 0.0),
    };
    assert!(surfaces_coplanar(&plane_a, &plane_b, 1e-6));
}

#[test]
fn surfaces_coplanar_non_planar_always_true() {
    let cyl_a = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let cyl_b = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 5.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 2.0,
    };
    assert!(surfaces_coplanar(&cyl_a, &cyl_b, 1e-6));
}

#[test]
fn surface_types_compatible_all_same_types() {
    let sphere_a = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: crate::geometry::Vec3::Z,
    };
    let sphere_b = SurfaceGeometry::Sphere {
        center: Point3::new(1.0, 0.0, 0.0),
        radius: 2.0,
        axis: crate::geometry::Vec3::Z,
    };
    assert!(surface_types_compatible(&sphere_a, &sphere_b));

    let cone_a = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: 0.5,
    };
    let cone_b = SurfaceGeometry::Cone {
        apex: Point3::new(1.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 1.0, 0.0),
        half_angle: 0.3,
    };
    assert!(surface_types_compatible(&cone_a, &cone_b));

    let torus_a = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 1.0,
    };
    let torus_b = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 3.0,
        minor_radius: 0.5,
    };
    assert!(surface_types_compatible(&torus_a, &torus_b));
}

#[test]
fn surface_types_incompatible_nurbs() {
    let nurbs = SurfaceGeometry::Nurbs(crate::geometry::nurbs_surface::NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0, 1.0, 1.0, 1.0],
    });
    assert!(!surface_types_compatible(&nurbs, &nurbs));
}

#[test]
fn merge_preserves_remaining_topology() {
    let mut brep = make_coplanar_pair();
    let tol = Tolerance::default();
    let initial_edges = brep.edges.len();

    let result = merge_coplanar_faces(&mut brep, &tol);
    assert_eq!(result.pairs_merged, 1);

    assert_eq!(brep.faces.len(), 1);
    assert_eq!(brep.edges.len(), initial_edges - 1);
}
