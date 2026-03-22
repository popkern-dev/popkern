#[cfg(test)]
use crate::error::ValidationError;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// Build a degree-1 bilinear NURBS patch (a flat quad) spanning a finite
/// parameter region.  The four corner control points are given in the order
/// `(u0,v0), (u0,v1), (u1,v0), (u1,v1)` (row-major, u varies in each row).
fn bilinear_patch(p00: Point3, p01: Point3, p10: Point3, p11: Point3) -> SurfaceGeometry {
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p00, p01, p10, p11],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0, 1.0, 1.0, 1.0],
    })
}

#[test]
fn cube_passes_full_self_intersection() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    assert!(cube.validate_self_intersection().is_ok());
}

#[test]
fn prism_passes_full_self_intersection() {
    let profile = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.5, 1.0, 0.0),
    ];
    let prism = BRep::prism(
        &profile,
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        Tolerance::default(),
    )
    .unwrap();
    assert!(prism.validate_self_intersection().is_ok());
}

#[test]
fn empty_brep_passes() {
    let brep = BRep::new(Tolerance::default());
    assert!(brep.validate_self_intersection().is_ok());
}

#[test]
fn single_face_brep_passes() {
    // A brep with only one face cannot self-intersect.
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let _mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
    assert!(brep.validate_self_intersection().is_ok());
}

/// Faces without surface geometry are skipped gracefully.
#[test]
fn face_without_surface_skips_gracefully() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face_keys: Vec<_> = brep.faces.keys().collect();

    // Remove surface from one face.
    brep.faces[face_keys[0]].surface = None;

    // Should still pass — faces without geometry are skipped.
    assert!(brep.validate_self_intersection().is_ok());
}

/// Test that adjacent faces in a cube are not reported as self-intersecting.
#[test]
fn adjacent_faces_not_reported() {
    // The cube has geometry on all faces; no non-adjacent faces share a plane.
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    assert!(cube.validate_self_intersection().is_ok());
}

/// Test the adjacency walking with a valid cube.
#[test]
fn unit_cube_with_geometry_no_self_intersection() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    assert!(cube.validate_self_intersection().is_ok());
}

/// Faces whose surfaces have no geometry (all None) cannot yield any BVH pairs.
#[test]
fn all_faces_without_surfaces_pass() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Remove all surface geometry.
    for fk in brep.faces.keys().collect::<Vec<_>>() {
        brep.faces[fk].surface = None;
    }
    assert!(brep.validate_self_intersection().is_ok());
}

/// Build a minimal two-face BRep (no shared edges) and attach crossing NURBS
/// bilinear patches so the self-intersection check fires the error path.
///
/// Topology: two completely disconnected faces created via mvfs twice —
/// they share no edges, so the adjacency set is empty.  `SurfaceGeometry::Nurbs`
/// is used because it has finite parameter ranges (from the knot vectors),
/// which means the AABB sampling in `validate_self_intersection` produces real
/// numbers, enabling the BVH to pair the two faces.
///
/// Patch A: flat quad in the XY plane (z=0), u∈[0,1], v∈[0,1].
///   corners: (0,0,0), (0,1,0), (1,0,0), (1,1,0).
/// Patch B: flat quad in the XZ plane (y=0), u∈[0,1], v∈[0,1].
///   corners: (0,0,0), (0,0,1), (1,0,0), (1,0,1).
///
/// Both patches share the segment x∈[0,1], y=0, z=0, so their AABBs overlap
/// and `intersect_surfaces_with_handler` finds an intersection line.
#[test]
fn self_intersecting_faces_detected() {
    let mut brep = BRep::new(Tolerance::default());

    // Face A — isolated topology node; carries patch A geometry.
    let mvfs_a = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev_a1 = brep
        .mev_at_vertex(mvfs_a.vertex, mvfs_a.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev_a2 = brep
        .mev(mev_a1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev_a3 = brep
        .mev(mev_a2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let mef_a = brep.mef(mev_a3.half_edge_twin, mev_a1.half_edge).unwrap();

    // Face B — isolated topology node; carries patch B geometry.
    let mvfs_b = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev_b1 = brep
        .mev_at_vertex(mvfs_b.vertex, mvfs_b.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev_b2 = brep
        .mev(mev_b1.half_edge_twin, Point3::new(1.0, 0.0, 1.0))
        .unwrap();
    let mev_b3 = brep
        .mev(mev_b2.half_edge_twin, Point3::new(0.0, 0.0, 1.0))
        .unwrap();
    let mef_b = brep.mef(mev_b3.half_edge_twin, mev_b1.half_edge).unwrap();

    // Patch A: flat quad in the XY plane (z = 0).
    //   row-major order: (u0,v0), (u0,v1), (u1,v0), (u1,v1)
    brep.faces[mef_a.face].surface = Some(SurfaceOnFace {
        geometry: bilinear_patch(
            Point3::new(0.0, 0.0, 0.0), // u=0, v=0
            Point3::new(0.0, 1.0, 0.0), // u=0, v=1
            Point3::new(1.0, 0.0, 0.0), // u=1, v=0
            Point3::new(1.0, 1.0, 0.0), // u=1, v=1
        ),
        orientation: true,
    });
    // Patch B: flat quad in the XZ plane (y = 0).
    //   intersects patch A along x∈[0,1], y=0, z=0.
    brep.faces[mef_b.face].surface = Some(SurfaceOnFace {
        geometry: bilinear_patch(
            Point3::new(0.0, 0.0, 0.0), // u=0, v=0
            Point3::new(0.0, 0.0, 1.0), // u=0, v=1
            Point3::new(1.0, 0.0, 0.0), // u=1, v=0
            Point3::new(1.0, 0.0, 1.0), // u=1, v=1
        ),
        orientation: true,
    });

    let result = brep.validate_self_intersection();
    assert!(
        matches!(result, Err(ValidationError::SelfIntersection { .. })),
        "expected SelfIntersection error, got: {result:?}"
    );
}

/// Verify that the tighter 2× tolerance catches a near-intersection that the
/// old 10× tolerance would have missed.
///
/// Two non-adjacent flat patches are placed 3× `tolerance.linear` apart in Z.
/// With 2× tolerance the intersection test uses a coarse tolerance of 2×, which
/// is smaller than the 3× gap — so the surfaces are not within that band and
/// the test **passes** (no intersection flagged).  The test serves as a
/// regression guard: if the multiplier were raised back to 10×, the 3× gap
/// would fall inside the 10× band and the intersection test might fire
/// spuriously, reversing the result.  For completeness we also check that
/// two patches separated by only 1× tolerance ARE flagged under the 2× regime.
#[test]
fn tighter_tolerance_catches_near_intersection() {
    let tol = Tolerance::default();
    let gap_3x = tol.linear * 3.0; // outside 2× band → should NOT be flagged
    let gap_half = tol.linear * 0.5; // inside 2× band → should be flagged

    // Helper: build a two-disconnected-face BRep with patch A in the XY plane
    // (z=0) and patch B offset by `z_offset` in Z.
    let make_brep = |z_offset: f64| -> BRep {
        let mut brep = BRep::new(tol);

        let mvfs_a = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev_a1 = brep
            .mev_at_vertex(mvfs_a.vertex, mvfs_a.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev_a2 = brep
            .mev(mev_a1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev_a3 = brep
            .mev(mev_a2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let mef_a = brep.mef(mev_a3.half_edge_twin, mev_a1.half_edge).unwrap();

        let mvfs_b = brep.mvfs(Point3::new(0.0, 0.0, z_offset)).unwrap();
        let mev_b1 = brep
            .mev_at_vertex(mvfs_b.vertex, mvfs_b.loop_, Point3::new(1.0, 0.0, z_offset))
            .unwrap();
        let mev_b2 = brep
            .mev(mev_b1.half_edge_twin, Point3::new(1.0, 1.0, z_offset))
            .unwrap();
        let mev_b3 = brep
            .mev(mev_b2.half_edge_twin, Point3::new(0.0, 1.0, z_offset))
            .unwrap();
        let mef_b = brep.mef(mev_b3.half_edge_twin, mev_b1.half_edge).unwrap();

        // Patch A: XY plane at z = 0.
        brep.faces[mef_a.face].surface = Some(SurfaceOnFace {
            geometry: bilinear_patch(
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ),
            orientation: true,
        });
        // Patch B: XY plane at z = z_offset (parallel, offset).
        brep.faces[mef_b.face].surface = Some(SurfaceOnFace {
            geometry: bilinear_patch(
                Point3::new(0.0, 0.0, z_offset),
                Point3::new(0.0, 1.0, z_offset),
                Point3::new(1.0, 0.0, z_offset),
                Point3::new(1.0, 1.0, z_offset),
            ),
            orientation: true,
        });

        brep
    };

    // Two parallel, well-separated patches (3× gap) should NOT be reported
    // as self-intersecting under the new 2× tolerance regime.
    let brep_far = make_brep(gap_3x);
    assert!(
        brep_far.validate_self_intersection().is_ok(),
        "parallel patches 3× apart should not be flagged as self-intersecting"
    );

    // Two parallel patches separated by only half a tolerance unit are so
    // close that the intersection check (using 2× coarse tolerance) treats
    // them as intersecting — they should be flagged.
    let brep_near = make_brep(gap_half);
    assert!(
        matches!(
            brep_near.validate_self_intersection(),
            Err(ValidationError::SelfIntersection { .. })
        ),
        "parallel patches 0.5× apart should be flagged as self-intersecting"
    );
}

/// When one face in a non-adjacent pair has `surface = None` the pair is
/// silently skipped (the `None => continue` branch in the intersection loop).
/// Build the same two-disconnected-face topology but leave face B's surface
/// as None — the check must pass.
#[test]
fn non_adjacent_pair_with_one_missing_surface_passes() {
    let mut brep = BRep::new(Tolerance::default());

    // Face A — bilinear NURBS patch in XY plane.
    let mvfs_a = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev_a1 = brep
        .mev_at_vertex(mvfs_a.vertex, mvfs_a.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev_a2 = brep
        .mev(mev_a1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .unwrap();
    let mev_a3 = brep
        .mev(mev_a2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .unwrap();
    let mef_a = brep.mef(mev_a3.half_edge_twin, mev_a1.half_edge).unwrap();

    // Face B — disconnected, surface stays None.
    let mvfs_b = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev_b1 = brep
        .mev_at_vertex(mvfs_b.vertex, mvfs_b.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let mev_b2 = brep
        .mev(mev_b1.half_edge_twin, Point3::new(1.0, 0.0, 1.0))
        .unwrap();
    let mev_b3 = brep
        .mev(mev_b2.half_edge_twin, Point3::new(0.0, 0.0, 1.0))
        .unwrap();
    let mef_b = brep.mef(mev_b3.half_edge_twin, mev_b1.half_edge).unwrap();

    // Only face A gets a surface; face B stays None.
    brep.faces[mef_a.face].surface = Some(SurfaceOnFace {
        geometry: bilinear_patch(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ),
        orientation: true,
    });
    let _ = mef_b; // face B surface stays None

    assert!(
        brep.validate_self_intersection().is_ok(),
        "pair with one missing surface should not be reported as self-intersecting"
    );
}
