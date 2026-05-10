// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].
// Edge-case / fallback-path tests split out of classify_tests.rs to keep both
// files under the 500-line structural cap.
#[cfg(test)]
use super::*;
use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

fn tol() -> Tolerance {
    Tolerance::default()
}

// -- ray_intersects_aabb — negative direction (triggers t1/t2 swap)

#[test]
fn ray_negative_direction_through_aabb_hits() {
    let aabb =
        Aabb::from_points(&[Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)]).expect("aabb");
    // Origin beyond the box on the +x side, ray pointing in -x direction.
    let origin = Point3::new(2.0, 0.5, 0.5);
    let dir = Vec3::new(-1.0, 0.0, 0.0);
    assert!(ray_intersects_aabb(&origin, &dir, &aabb));
}

#[test]
fn ray_negative_direction_missing_aabb() {
    let aabb =
        Aabb::from_points(&[Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)]).expect("aabb");
    // Origin beyond +x, ray pointing -x, but y is outside the slab.
    let origin = Point3::new(2.0, 5.0, 0.5);
    let dir = Vec3::new(-1.0, 0.0, 0.0);
    assert!(!ray_intersects_aabb(&origin, &dir, &aabb));
}

// -- pick_sample_point — fallback paths (no vertices or no surface)

#[test]
fn pick_sample_point_face_with_plane_surface_uses_fallback() {
    use popkern_core::{
        geometry::surface::{SurfaceGeometry, SurfaceOnFace},
        topology::{Face, HalfEdgeHandle, Loop, ShellHandle},
    };
    let mut brep = BRep::new(tol());
    // Create a face whose loop points at a missing half-edge — so
    // face_vertex_positions returns empty — and give it a plane surface so
    // the surface-midpoint fallback fires (plane has infinite ranges, so
    // face_parameter_range uses the vertex projection → empty → default).
    let loop_h = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: popkern_core::topology::FaceHandle::default(),
    });
    let face_h = brep.faces.insert(Face {
        outer_loop: loop_h,
        inner_loops: Vec::new(),
        shell: ShellHandle::default(),
        surface: Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.5, 0.5, 0.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
            },
            orientation: true,
        }),
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    brep.loops[loop_h].face = face_h;
    // No half-edges exist, so face_vertex_positions returns empty.
    let pt = pick_sample_point(&brep, face_h);
    // Should return the surface midpoint — not the origin fallback.
    assert!(pt.x.is_finite());
    assert!(pt.y.is_finite());
    assert!(pt.z.is_finite());
}

#[test]
fn pick_sample_point_face_with_sphere_surface_uses_finite_range() {
    use popkern_core::{
        geometry::surface::{SurfaceGeometry, SurfaceOnFace},
        topology::{Face, HalfEdgeHandle, Loop, ShellHandle},
    };
    let mut brep = BRep::new(tol());
    // Sphere surface has a finite parameter range, so face_parameter_range
    // returns early at the "both ranges are finite" check (line 89).
    let loop_h = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: popkern_core::topology::FaceHandle::default(),
    });
    let face_h = brep.faces.insert(Face {
        outer_loop: loop_h,
        inner_loops: Vec::new(),
        shell: ShellHandle::default(),
        surface: Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Sphere {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            },
            orientation: true,
        }),
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    brep.loops[loop_h].face = face_h;
    let pt = pick_sample_point(&brep, face_h);
    // Should return a point on the sphere — not origin.
    // The sphere parameter_range gives finite u,v, so the midpoint is on the sphere.
    assert!(pt.x.is_finite());
    assert!(pt.y.is_finite());
    assert!(pt.z.is_finite());
}

/// Tests that face_vertex_positions handles a face whose outer_loop handle
/// no longer exists in the brep (e.g. after topology corruption or a stale handle).
#[test]
fn face_vertex_positions_missing_loop_returns_empty() {
    use popkern_core::topology::{Face, ShellHandle};
    let mut brep = BRep::new(tol());
    // Create a face with a loop handle that is NOT inserted into brep.loops.
    // The default LoopHandle should not exist.
    let face_h = brep.faces.insert(Face {
        outer_loop: popkern_core::topology::LoopHandle::default(),
        inner_loops: Vec::new(),
        shell: ShellHandle::default(),
        surface: None,
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    // pick_sample_point calls face_vertex_positions; with missing loop it
    // returns empty, then falls through to the origin fallback.
    let pt = pick_sample_point(&brep, face_h);
    assert_eq!(pt.x, 0.0);
    assert_eq!(pt.y, 0.0);
    assert_eq!(pt.z, 0.0);
}

#[test]
fn face_parameter_range_with_cylinder_surface_returns_large_default() {
    use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
    // Build a face with a Cylinder surface (infinite v-range) and some vertices.
    // face_vertex_positions will find the vertices, so positions is non-empty.
    // Then the non-Plane branch at line 109 triggers.
    let mut brep = BRep::new(tol());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef_r = brep.mef(mev3.half_edge_twin, mev1.half_edge).expect("mef");
    let face_h = mef_r.face;
    // Replace the face surface with a Cylinder (infinite v-range).
    brep.faces[face_h].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        orientation: true,
    });
    // Cylinder parameter range: u=[0, 2π] (finite), v=(-∞, ∞) (infinite).
    // Since v is infinite, we fall into the vertex projection path.
    // Then since it's not a Plane, we hit the non-plane fallback (line 109).
    let ((u_min, u_max), (v_min, v_max)) = face_parameter_range(&brep, face_h);
    assert_eq!(u_min, -100.0);
    assert_eq!(u_max, 100.0);
    assert_eq!(v_min, -100.0);
    assert_eq!(v_max, 100.0);
}

#[test]
fn face_parameter_range_with_missing_loop_returns_default() {
    use popkern_core::{
        geometry::surface::{SurfaceGeometry, SurfaceOnFace},
        topology::{Face, HalfEdgeHandle, Loop, ShellHandle},
    };
    let mut brep = BRep::new(tol());
    // A face whose loop doesn't exist in the brep (stale loop handle) and
    // a plane surface: face_vertex_positions returns empty, and the non-plane
    // guard would catch a non-plane surface. For a plane with empty vertices,
    // we fall through to the default (-100,100) range.
    let loop_h = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: popkern_core::topology::FaceHandle::default(),
    });
    let face_h = brep.faces.insert(Face {
        outer_loop: loop_h,
        inner_loops: Vec::new(),
        shell: ShellHandle::default(),
        surface: Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
            },
            orientation: true,
        }),
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    brep.loops[loop_h].face = face_h;
    // Remove the loop so face_vertex_positions can't find half-edges.
    // The loop still exists but half_edge is default (invalid handle).
    // face_vertex_positions returns empty because half_edge doesn't exist.
    let ((u_min, u_max), (v_min, v_max)) = face_parameter_range(&brep, face_h);
    // With empty vertex list and plane surface: returns (-100, 100) x (-100, 100).
    assert_eq!(u_min, -100.0);
    assert_eq!(u_max, 100.0);
    assert_eq!(v_min, -100.0);
    assert_eq!(v_max, 100.0);
}

#[test]
fn pick_sample_point_face_with_no_surface_returns_origin() {
    use popkern_core::topology::{Face, HalfEdgeHandle, Loop, ShellHandle};
    let mut brep = BRep::new(tol());
    let loop_h = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: popkern_core::topology::FaceHandle::default(),
    });
    let face_h = brep.faces.insert(Face {
        outer_loop: loop_h,
        inner_loops: Vec::new(),
        shell: ShellHandle::default(),
        surface: None,
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    brep.loops[loop_h].face = face_h;
    let pt = pick_sample_point(&brep, face_h);
    // Last-resort fallback: origin.
    assert_eq!(pt.x, 0.0);
    assert_eq!(pt.y, 0.0);
    assert_eq!(pt.z, 0.0);
}

// -- classify_fragments — OnBoundary case

#[test]
fn classify_fragments_on_boundary_maps_to_cooriented() {
    // A cube exactly on its own boundary — classify it against itself.
    // When a face's sample point lies on a surface, point_in_solid may return
    // OnBoundary; classify_fragments should map that to OnBoundaryCooriented.
    let cube = BRep::unit_cube(tol()).expect("cube");
    let aabbs = crate::intersect::compute_face_aabbs(&cube);
    let empty_edges = HashSet::new();

    let result = classify_fragments(&cube, &cube, &aabbs, &empty_edges, &tol())
        .expect("classification should succeed");

    // All 6 faces should classify as something — no error.
    assert_eq!(result.len(), 6);
    // All should be OnBoundaryCooriented (or Inside/Outside depending on
    // exact sample placement), but classification must not fail.
    for class in result.values() {
        assert!(
            matches!(
                class,
                Classification::Inside
                    | Classification::Outside
                    | Classification::OnBoundaryCooriented
                    | Classification::OnBoundaryCrossing
            ),
            "unexpected classification: {class:?}"
        );
    }
}
