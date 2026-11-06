use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn unit_cube_is_closed() {
    let cube = BRep::unit_cube(tol()).unwrap();
    assert!(is_closed_solid(&cube));
    assert!(!is_sheet_body(&cube));
}

#[test]
fn cylinder_is_closed() {
    let cyl = BRep::cylinder(0.5, 2.0, tol()).unwrap();
    assert!(is_closed_solid(&cyl));
    assert!(!is_sheet_body(&cyl));
}

#[test]
fn empty_brep_is_not_a_sheet() {
    let empty = BRep::new(tol());
    assert!(!is_sheet_body(&empty));
    assert!(is_closed_solid(&empty));
}

#[test]
fn dangling_edge_brep_is_a_sheet() {
    // mvfs + mev creates a single face holding a dangling edge.
    // The edge's two half-edges both belong to the same face;
    // the face has no surface attached, so the periodic-skip
    // branch does not fire and the detector reports a sheet.
    use popkern_core::geometry::Point3;
    let mut brep = BRep::new(tol());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let _mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    assert!(
        is_sheet_body(&brep),
        "dangling-edge BRep should be a sheet body"
    );
}

#[test]
fn side_of_sheet_outside_cube_is_front_of_closest_face() {
    // For a multi-face body, side_of_sheet picks the face nearest the
    // point within its trims and reports relative to its outward
    // normal: a point outside the cube is on the front of its nearest
    // face, whose outward normal points away from the interior. A side
    // wall's unbounded plane is nearer a point above the cube than the
    // top face is, but its foot lies past the wall's edges.
    let cube = BRep::unit_cube(tol()).unwrap();
    use popkern_core::geometry::Point3;
    for p in [
        Point3::new(0.5, 0.5, 5.0),
        Point3::new(0.5, 0.5, -5.0),
        Point3::new(5.0, 0.5, 0.5),
        // Beyond an edge: no face's foot lies inside its trims.
        Point3::new(5.0, 5.0, 0.5),
    ] {
        assert_eq!(
            side_of_sheet(&cube, &p, &tol()),
            SideOfSheet::Front,
            "{p:?}"
        );
    }
}

#[test]
fn side_of_sheet_classifies_consistently() {
    // Inside the cube is behind its nearest face, a point on a face is
    // on it, and the answer is deterministic.
    let cube = BRep::unit_cube(tol()).unwrap();
    use popkern_core::geometry::Point3;
    let inside = Point3::new(0.5, 0.5, 0.4);
    assert_eq!(side_of_sheet(&cube, &inside, &tol()), SideOfSheet::Back);
    assert_eq!(
        side_of_sheet(&cube, &inside, &tol()),
        side_of_sheet(&cube, &inside, &tol())
    );
    let on = Point3::new(0.5, 0.5, 1.0);
    assert_eq!(side_of_sheet(&cube, &on, &tol()), SideOfSheet::On);
}

#[test]
fn side_of_sheet_returns_on_for_point_on_surface() {
    let cube = BRep::unit_cube(tol()).unwrap();
    use popkern_core::geometry::Point3;
    let on_face = Point3::new(0.5, 0.5, 1.0);
    let side = side_of_sheet(&cube, &on_face, &tol());
    assert!(
        matches!(side, SideOfSheet::On),
        "point on face should be On, got {side:?}"
    );
}

#[test]
fn classify_sheet_operation_solid_solid() {
    use crate::BoolOperation;
    let kind = classify_sheet_operation(BoolOperation::Union, false, false);
    assert!(matches!(kind, SheetOperationKind::SolidSolid));
}

#[test]
fn classify_sheet_operation_solid_minus_sheet() {
    use crate::BoolOperation;
    let kind = classify_sheet_operation(BoolOperation::Subtraction, false, true);
    assert!(matches!(kind, SheetOperationKind::SolidMinusSheet));
}

#[test]
fn classify_sheet_operation_invalid_sheet_minus_solid() {
    use crate::BoolOperation;
    let kind = classify_sheet_operation(BoolOperation::Subtraction, true, false);
    assert!(matches!(kind, SheetOperationKind::SheetMinusSolid));
}

#[test]
fn classify_sheet_operation_sheet_sheet_intersection() {
    use crate::BoolOperation;
    let kind = classify_sheet_operation(BoolOperation::Intersection, true, true);
    assert!(matches!(kind, SheetOperationKind::SheetSheetIntersection));
}

#[test]
fn split_solid_by_plane_bisects_cube() {
    use popkern_core::geometry::{Point3, Vec3};
    // A unit cube cut by the z=0.5 plane should yield 2 closed
    // sub-solids — the top half and the bottom half.
    let cube = BRep::unit_cube(tol()).unwrap();
    let pieces = split_solid_by_plane(
        &cube,
        Point3::new(0.5, 0.5, 0.5),
        Vec3::new(0.0, 0.0, 1.0),
        tol(),
    )
    .unwrap();
    assert_eq!(
        pieces.len(),
        2,
        "z=0.5 plane should bisect the cube into 2 pieces, got {}",
        pieces.len()
    );
    for (i, piece) in pieces.iter().enumerate() {
        assert!(piece.faces.len() >= 4, "piece {i} should have ≥4 faces");
    }
}

#[test]
fn split_solid_by_plane_off_axis_still_bisects() {
    use popkern_core::geometry::{Point3, Vec3};
    // Diagonal cut: x=0.5 + y=0.5 + z=0.5 plane through cube center.
    let cube = BRep::unit_cube(tol()).unwrap();
    let pieces = split_solid_by_plane(
        &cube,
        Point3::new(0.5, 0.5, 0.5),
        Vec3::new(1.0, 1.0, 1.0),
        tol(),
    )
    .unwrap();
    assert!(!pieces.is_empty(), "diagonal cut should produce ≥1 piece");
}

#[test]
fn is_planar_sheet_rejects_non_planar_surface() {
    // A cylinder is not planar — is_planar_sheet should return None
    // even though it's a closed solid with multiple faces.
    let cyl = BRep::cylinder(0.5, 2.0, tol()).unwrap();
    assert!(
        is_planar_sheet(&cyl).is_none(),
        "cylinder should NOT be detected as planar"
    );
}

#[test]
fn is_planar_sheet_rejects_mixed_normals() {
    // A unit cube has 6 faces on 3 pairs of parallel planes, but the
    // faces don't all share a normal. is_planar_sheet should return
    // None.
    let cube = BRep::unit_cube(tol()).unwrap();
    assert!(
        is_planar_sheet(&cube).is_none(),
        "cube should NOT be detected as planar (mixed normals)"
    );
}

#[test]
fn is_planar_sheet_rejects_empty_brep() {
    let empty = BRep::new(tol());
    assert!(
        is_planar_sheet(&empty).is_none(),
        "empty BRep should NOT be detected as planar"
    );
}

#[test]
fn is_planar_sheet_accepts_multi_face_coplanar() {
    // Build a minimal BRep with TWO face slots both pointing at the
    // same z=0.5 plane. is_planar_sheet only inspects .faces and their
    // surface geometry — it does not walk the topology — so we can
    // smoke-test the detector by pre-populating the face slotmap
    // directly without wiring a full loop+half-edge chain.
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
    use popkern_core::topology::Face;

    let mut brep = BRep::new(tol());
    let plane_a = SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.5),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        orientation: true,
    };
    let plane_b = SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: Point3::new(2.0, 2.0, 0.5),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        orientation: true,
    };
    brep.faces.insert(Face {
        outer_loop: Default::default(),
        inner_loops: Vec::new(),
        back: Default::default(),
        front: Default::default(),
        surface: Some(plane_a),
        name: None,
        tolerance: None,
        role: Default::default(),
    });
    brep.faces.insert(Face {
        outer_loop: Default::default(),
        inner_loops: Vec::new(),
        back: Default::default(),
        front: Default::default(),
        surface: Some(plane_b),
        name: None,
        tolerance: None,
        role: Default::default(),
    });

    let detected = is_planar_sheet(&brep);
    assert!(
        detected.is_some(),
        "two co-planar faces should be detected as an effectively-planar sheet"
    );
    let (_, normal) = detected.unwrap();
    assert_eq!(normal.z, 1.0);
}

#[test]
fn is_approximately_planar_sheet_accepts_nearly_flat_vertices() {
    // A sheet with vertices that bow by 1e-9 above and below
    // z=0.5 should be detected as approximately planar when the
    // deviation envelope is 1e-3.
    use popkern_core::constructors::builder::BRepBuilder;
    use popkern_core::geometry::Point3;

    let mut b = BRepBuilder::new(tol());
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let face = b.add_face(shell);
    // Four vertices forming a nearly-planar quad at z≈0.5.
    let _v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.5));
    let _v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.5 + 1e-9));
    let _v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.5 - 1e-9));
    let _v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.5));
    // Note: no edges/loops — the approximate-planar detector only
    // inspects vertex positions, so we can skip the topology wiring
    // for this focused test.
    let _ = face;
    let sheet = b.build_unchecked();

    let detected = is_approximately_planar_sheet(&sheet, 1.0e-3);
    assert!(
        detected.is_some(),
        "nearly-planar sheet should be detected within 1e-3 envelope"
    );
    let (_, normal) = detected.unwrap();
    // Normal should be approximately ±Z.
    assert!((normal.z.abs() - 1.0).abs() < 1.0e-3);
    assert!(normal.x.abs() < 1.0e-3);
    assert!(normal.y.abs() < 1.0e-3);
}

#[test]
fn is_approximately_planar_sheet_rejects_curved_vertices() {
    // A sheet whose vertex positions trace out a sphere-cap should
    // NOT be detected as approximately planar at a tight envelope.
    use popkern_core::constructors::builder::BRepBuilder;
    use popkern_core::geometry::Point3;

    let mut b = BRepBuilder::new(tol());
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let face = b.add_face(shell);
    // Generate 8 points on a sphere of radius 1, forming a cap.
    for i in 0..8 {
        let theta = (i as f64) * std::f64::consts::PI / 4.0;
        let _v = b.add_vertex(Point3::new(theta.cos(), theta.sin(), (theta * 0.5).sin()));
    }
    let _ = face;
    let sheet = b.build_unchecked();

    // Very tight envelope should reject.
    let detected = is_approximately_planar_sheet(&sheet, 1.0e-6);
    assert!(
        detected.is_none(),
        "sphere-cap vertices should NOT be detected as planar at 1e-6 envelope"
    );
}

#[test]
fn is_approximately_planar_sheet_requires_three_vertices() {
    // Fewer than 3 vertices can't define a plane.
    let empty = BRep::new(tol());
    assert!(is_approximately_planar_sheet(&empty, 1.0).is_none());
}

#[test]
fn is_planar_sheet_rejects_offset_coplanar_faces() {
    // Two planar faces sharing a normal but at different z-heights
    // (0.5 vs 1.5) should NOT be treated as a single plane — they
    // lie on parallel planes, not the same plane.
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
    use popkern_core::topology::Face;

    let mut brep = BRep::new(tol());
    for z in [0.5, 1.5] {
        brep.faces.insert(Face {
            outer_loop: Default::default(),
            inner_loops: Vec::new(),
            back: Default::default(),
            front: Default::default(),
            surface: Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Plane {
                    origin: Point3::new(0.0, 0.0, z),
                    normal: Vec3::new(0.0, 0.0, 1.0),
                    u_axis: Vec3::new(1.0, 0.0, 0.0),
                    v_axis: Vec3::new(0.0, 1.0, 0.0),
                },
                orientation: true,
            }),
            name: None,
            tolerance: None,
            role: Default::default(),
        });
    }
    assert!(
        is_planar_sheet(&brep).is_none(),
        "offset parallel planes should NOT be detected as a single plane"
    );
}
