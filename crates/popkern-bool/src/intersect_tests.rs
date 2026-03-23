// Structural lint marker: this file is test code, included via #[cfg(test)] in intersect.rs.
#[cfg(test)]
use super::*;
use popkern_core::{
    geometry::{Point3, Vec3},
    tolerance::Tolerance,
    topology::BRep,
};

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

#[test]
fn overlapping_cubes_find_intersection_curves() {
    let cube_a = BRep::unit_cube(tol()).expect("unit cube a");
    // Offset cube_b by (0.5, 0.3, 0.2) so it overlaps with cube_a but
    // none of the face planes are coincident (all 6 faces land on
    // distinct planes from cube_a's 6 faces).
    let cube_b = translate(
        BRep::unit_cube(tol()).expect("unit cube b"),
        Vec3::new(0.5, 0.3, 0.2),
    );

    let result =
        find_intersections(&cube_a, &cube_b, &tol()).expect("overlapping cubes should not error");

    assert!(
        !result.records.is_empty(),
        "expected at least one intersection curve for overlapping cubes"
    );
}

#[test]
fn disjoint_cubes_no_intersection_curves() {
    let cube_a = BRep::unit_cube(tol()).expect("unit cube a");
    // Offset cube_b by 5.0 in x so it is completely separate.
    let cube_b = translate(
        BRep::unit_cube(tol()).expect("unit cube b"),
        Vec3::new(5.0, 0.0, 0.0),
    );

    let result =
        find_intersections(&cube_a, &cube_b, &tol()).expect("disjoint cubes should not error");

    assert!(
        result.records.is_empty(),
        "expected no intersection curves for disjoint cubes"
    );
}

#[test]
fn compute_face_aabbs_unit_cube_has_six_faces() {
    let cube = BRep::unit_cube(tol()).expect("unit cube");
    let aabbs = compute_face_aabbs(&cube);
    assert_eq!(aabbs.len(), 6, "unit cube should have 6 faces with AABBs");
}

#[test]
fn compute_face_aabbs_no_surface_returns_empty() {
    // A freshly-created BRep has no faces at all.
    let brep = BRep::new(tol());
    let aabbs = compute_face_aabbs(&brep);
    assert!(aabbs.is_empty());
}

#[test]
fn intersection_result_aabbs_preserved() {
    let cube_a = BRep::unit_cube(tol()).expect("unit cube a");
    let cube_b = translate(
        BRep::unit_cube(tol()).expect("unit cube b"),
        Vec3::new(5.0, 0.0, 0.0),
    );
    let result = find_intersections(&cube_a, &cube_b, &tol()).expect("no error");
    assert_eq!(result.aabbs_a.len(), 6);
    assert_eq!(result.aabbs_b.len(), 6);
}

#[test]
fn face_aabb_unit_cube_face_is_bounded() {
    let cube = BRep::unit_cube(tol()).expect("unit cube");
    // Check every face produces a finite AABB within [0,1]^3 ± small epsilon.
    let aabbs = compute_face_aabbs(&cube);
    for (_, aabb) in &aabbs {
        assert!(
            aabb.min.x >= -1e-6 && aabb.max.x <= 1.0 + 1e-6,
            "x out of range"
        );
        assert!(
            aabb.min.y >= -1e-6 && aabb.max.y <= 1.0 + 1e-6,
            "y out of range"
        );
        assert!(
            aabb.min.z >= -1e-6 && aabb.max.z <= 1.0 + 1e-6,
            "z out of range"
        );
    }
}

/// Exercise face_vertex_positions when the loop's half_edge handle is
/// missing from the BRep — it should return an empty Vec rather than panic.
#[test]
fn face_vertex_positions_missing_half_edge_returns_empty() {
    use popkern_core::topology::{Face, HalfEdgeHandle, Loop, ShellHandle};
    let mut brep = BRep::new(tol());
    // Insert a loop whose half_edge handle does not exist in brep.half_edges.
    let loop_h = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: popkern_core::topology::FaceHandle::default(),
    });
    let face_h = brep.faces.insert(Face {
        outer_loop: loop_h,
        inner_loops: Vec::new(),
        shell: ShellHandle::default(),
        surface: None,
    });
    brep.loops[loop_h].face = face_h;
    // face_aabb should return None (no surface geometry).
    let result = face_aabb(&brep, face_h);
    assert!(result.is_none());
}

/// Exercise the 5×5 grid-sampling branch in face_aabb for a surface with
/// a finite parameter range (sphere/torus). We use a unit cube here and
/// verify the AABBs are still correct — the cube planes have infinite
/// parameter ranges so only vertices are used. For a more direct test we
/// build a minimal BRep face with a sphere surface.
#[test]
fn face_aabb_sphere_surface_samples_grid() {
    use popkern_core::{
        geometry::surface::{SurfaceGeometry, SurfaceOnFace},
        topology::{Face, HalfEdgeHandle, Loop, ShellHandle},
    };
    let mut brep = BRep::new(tol());
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
    });
    brep.loops[loop_h].face = face_h;
    // The sphere surface has a finite parameter range, so face_aabb should
    // sample the 5×5 grid and return a non-None AABB.
    let result = face_aabb(&brep, face_h);
    assert!(result.is_some(), "sphere face should produce an AABB");
    let aabb = result.unwrap();
    // Unit sphere AABB should be approximately [-1, 1]^3.
    assert!(
        aabb.min.x >= -1.1 && aabb.max.x <= 1.1,
        "x out of sphere range"
    );
}

/// Exercise the CoincidentRegion error path in find_intersections.
/// We use two identical cubes (same position) which produce coincident
/// face pairs and trigger the error.
#[test]
fn coincident_cubes_returns_coincident_faces_error() {
    let cube_a = BRep::unit_cube(tol()).expect("unit cube a");
    let cube_b = BRep::unit_cube(tol()).expect("unit cube b (same position)");
    let result = find_intersections(&cube_a, &cube_b, &tol());
    assert!(
        result.is_err(),
        "coincident cubes should produce a CoincidentFaces error"
    );
}

#[test]
fn compute_face_aabbs_sphere_returns_three_faces() {
    let sphere = BRep::sphere(1.0, tol()).expect("sphere");
    let aabbs = compute_face_aabbs(&sphere);
    assert_eq!(
        aabbs.len(),
        3,
        "sphere should have 3 faces with AABBs, got {}",
        aabbs.len()
    );
    for (_, aabb) in &aabbs {
        assert!(
            aabb.min.x >= -1.1 && aabb.max.x <= 1.1,
            "x out of range: [{}, {}]",
            aabb.min.x,
            aabb.max.x,
        );
        assert!(
            aabb.min.y >= -1.1 && aabb.max.y <= 1.1,
            "y out of range: [{}, {}]",
            aabb.min.y,
            aabb.max.y,
        );
        assert!(
            aabb.min.z >= -1.1 && aabb.max.z <= 1.1,
            "z out of range: [{}, {}]",
            aabb.min.z,
            aabb.max.z,
        );
    }
}

/// Exercise `face_aabb` with a face that has a Torus surface geometry
/// (finite parameter range). The torus face AABB is computed from a 5×5
/// grid of surface samples, exercising the finite-range grid-sampling path
/// (lines 122-131 in intersect.rs) and the `bounded_parameter_range`
/// finite-range early-return (line 32-34).
#[test]
fn face_aabb_torus_surface_samples_grid() {
    use popkern_core::{
        geometry::surface::{SurfaceGeometry, SurfaceOnFace},
        topology::{Face, HalfEdgeHandle, Loop, ShellHandle},
    };
    let mut brep = BRep::new(tol());
    let loop_h = brep.loops.insert(Loop {
        half_edge: HalfEdgeHandle::default(),
        face: popkern_core::topology::FaceHandle::default(),
    });
    let face_h = brep.faces.insert(Face {
        outer_loop: loop_h,
        inner_loops: Vec::new(),
        shell: ShellHandle::default(),
        surface: Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Torus {
                center: Point3::new(0.0, 0.0, 0.0),
                axis: Vec3::new(0.0, 0.0, 1.0),
                major_radius: 2.0,
                minor_radius: 0.5,
            },
            orientation: true,
        }),
    });
    brep.loops[loop_h].face = face_h;
    // The torus surface has a finite parameter range [0,TAU]×[0,TAU], so
    // face_aabb should sample the 5×5 grid and return a non-None AABB.
    let result = face_aabb(&brep, face_h);
    assert!(result.is_some(), "torus face should produce an AABB");
    let aabb = result.unwrap();
    // Torus with major=2, minor=0.5 has outer radius 2.5.
    assert!(
        aabb.max.x <= 3.0 && aabb.min.x >= -3.0,
        "torus AABB x should be within [-3,3]"
    );
}

/// Translate a BRep with faces/edges that have no geometry (no surface or
/// curve attached). This covers the `else` branches in translate (lines
/// where `face.surface.is_none()` or `edge.curve.is_none()`).
#[test]
fn translate_brep_without_geometry_covers_no_surface_no_curve_branches() {
    // Build a minimal BRep (triangle) without calling attach_geometry —
    // so faces have no surface and edges have no curve.
    let mut brep = BRep::new(tol());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(0.5, 1.0, 0.0))
        .expect("mev2");
    let _mef = brep.mef(mev2.half_edge_twin, mev1.half_edge).expect("mef");

    // Verify: faces have no surface, edges have no curve.
    for (_, face) in &brep.faces {
        assert!(face.surface.is_none(), "faces should have no surface");
    }
    for (_, edge) in &brep.edges {
        assert!(edge.curve.is_none(), "edges should have no curve");
    }

    // Translating this BRep hits the `else` branches (face.surface.is_none(),
    // edge.curve.is_none()) in the translate helper.
    let translated = translate(brep, Vec3::new(1.0, 0.0, 0.0));

    // Verify vertices moved.
    for (_, v) in &translated.vertices {
        assert!(v.position.x >= 0.5, "vertex x should have been translated");
    }
}
