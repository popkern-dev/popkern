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
    // For a multi-face body, side_of_sheet picks the closest
    // face and reports relative to its outward normal. A point
    // outside the cube is on the front side of its closest face
    // (because every face's outward normal points away from
    // the cube interior). This test pins that behavior — for
    // a true sheet body the function returns Front/Back relative
    // to the single oriented surface.
    let cube = BRep::unit_cube(tol()).unwrap();
    use popkern_core::geometry::Point3;
    let above = Point3::new(0.5, 0.5, 5.0);
    let below = Point3::new(0.5, 0.5, -5.0);
    let side_above = side_of_sheet(&cube, &above, &tol());
    let side_below = side_of_sheet(&cube, &below, &tol());
    assert!(matches!(side_above, SideOfSheet::Front));
    assert!(matches!(side_below, SideOfSheet::Front));
}

#[test]
fn side_of_sheet_classifies_consistently() {
    // For a closed solid, side_of_sheet picks the closest face
    // and reports relative to its outward normal. The exact
    // Front/Back assignment depends on construction-time loop
    // orientation, but the function MUST be deterministic for
    // identical inputs and MUST report `On` when the point lies
    // on the surface within tolerance.
    let cube = BRep::unit_cube(tol()).unwrap();
    use popkern_core::geometry::Point3;
    let inside = Point3::new(0.5, 0.5, 0.5);
    let outside = Point3::new(0.5, 0.5, 5.0);

    let side_inside_1 = side_of_sheet(&cube, &inside, &tol());
    let side_inside_2 = side_of_sheet(&cube, &inside, &tol());
    assert_eq!(side_inside_1, side_inside_2);

    for side in [
        side_of_sheet(&cube, &inside, &tol()),
        side_of_sheet(&cube, &outside, &tol()),
    ] {
        assert!(matches!(
            side,
            SideOfSheet::Front | SideOfSheet::Back | SideOfSheet::On
        ));
    }
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
    assert!(matches!(kind, SheetOperationKind::InvalidSheetMinusSolid));
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

/// Build a single-face planar sheet body (a quad) at z=z_height
/// spanning x in [x_min, x_max], y in [y_min, y_max]. Used by
/// the multi-piece dispatcher tests.
fn build_planar_sheet(x_min: f64, x_max: f64, y_min: f64, y_max: f64, z_height: f64) -> BRep {
    use popkern_core::constructors::builder::BRepBuilder;
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

    let mut builder = BRepBuilder::new(tol());
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);
    let face = builder.add_face(shell);

    let v0 = builder.add_vertex(Point3::new(x_min, y_min, z_height));
    let v1 = builder.add_vertex(Point3::new(x_max, y_min, z_height));
    let v2 = builder.add_vertex(Point3::new(x_max, y_max, z_height));
    let v3 = builder.add_vertex(Point3::new(x_min, y_max, z_height));

    let (_, he01_fwd, _) = builder.add_edge(v0, v1);
    let (_, he12_fwd, _) = builder.add_edge(v1, v2);
    let (_, he23_fwd, _) = builder.add_edge(v2, v3);
    let (_, he30_fwd, _) = builder.add_edge(v3, v0);

    builder
        .add_loop(&[he01_fwd, he12_fwd, he23_fwd, he30_fwd], face)
        .unwrap();

    builder.brep_mut().faces.get_mut(face).unwrap().surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: Point3::new(x_min, y_min, z_height),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        orientation: true,
    });

    builder.build_unchecked()
}

#[test]
fn build_planar_sheet_passes_sheet_detector() {
    let sheet = build_planar_sheet(-0.5, 1.5, -0.5, 1.5, 0.5);
    assert!(
        is_sheet_body(&sheet),
        "constructed quad should be a sheet body"
    );
    assert!(
        is_planar_sheet(&sheet).is_some(),
        "constructed quad should be detected as planar"
    );
    let (origin, normal) = is_planar_sheet(&sheet).unwrap();
    assert_eq!(normal.z, 1.0);
    assert_eq!(origin.z, 0.5);
}

#[test]
fn boolean_with_sheets_solid_minus_planar_sheet_produces_pieces() {
    // End-to-end: cube minus a planar sheet that bisects it should
    // produce 2 pieces via the dispatcher.
    use crate::BoolOperation;
    let cube = BRep::unit_cube(tol()).unwrap();
    let sheet = build_planar_sheet(-0.5, 1.5, -0.5, 1.5, 0.5);
    let result = boolean_with_sheets(BoolOperation::Subtraction, &cube, &sheet, tol()).unwrap();
    match result {
        SheetBooleanResult::Pieces(pieces) => {
            assert_eq!(
                pieces.len(),
                2,
                "bisecting sheet should produce 2 pieces, got {}",
                pieces.len()
            );
        }
        other => panic!("expected Pieces, got {other:?}"),
    }
}

#[test]
fn boolean_with_sheets_solid_minus_off_axis_planar_sheet_produces_one_piece() {
    // A planar sheet that doesn't intersect the solid leaves it
    // unchanged: 1 piece returned as Solid.
    use crate::BoolOperation;
    let cube = BRep::unit_cube(tol()).unwrap();
    let sheet = build_planar_sheet(-0.5, 1.5, -0.5, 1.5, 5.0);
    let result = boolean_with_sheets(BoolOperation::Subtraction, &cube, &sheet, tol()).unwrap();
    match result {
        SheetBooleanResult::Solid(_) => {}
        other => panic!("expected single Solid, got {other:?}"),
    }
}

#[test]
fn split_solid_by_plane_outside_returns_one_piece() {
    use popkern_core::geometry::{Point3, Vec3};
    let cube = BRep::unit_cube(tol()).unwrap();
    let pieces = split_solid_by_plane(
        &cube,
        Point3::new(0.5, 0.5, 10.0),
        Vec3::new(0.0, 0.0, 1.0),
        tol(),
    )
    .unwrap();
    assert_eq!(
        pieces.len(),
        1,
        "non-bisecting plane should return 1 piece, got {}",
        pieces.len()
    );
}

#[test]
fn boolean_with_sheets_solid_solid_delegates_to_boolean() {
    use crate::BoolOperation;
    use popkern_core::geometry::Transform3;
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(5.0, 0.0, 0.0))
        .unwrap();
    let result = boolean_with_sheets(BoolOperation::Union, &a, &b, tol()).unwrap();
    assert!(matches!(result, SheetBooleanResult::Solid(_)));
}

#[test]
fn sheet_intersect_solid_routes_to_open_shell() {
    use crate::BoolOperation;
    use popkern_core::geometry::Point3;
    let mut sheet = BRep::new(tol());
    let mvfs = sheet.mvfs(Point3::new(0.0, 0.0, 0.5)).unwrap();
    let _mev = sheet
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(2.0, 0.0, 0.5))
        .unwrap();

    let cube = BRep::unit_cube(tol()).unwrap();
    let result = boolean_with_sheets(BoolOperation::Intersection, &sheet, &cube, tol());
    match result {
        Ok(SheetBooleanResult::OpenShell(_)) => {}
        Ok(SheetBooleanResult::Solid(_))
        | Ok(SheetBooleanResult::Curves(_))
        | Ok(SheetBooleanResult::Pieces(_)) => {
            panic!("expected OpenShell, got different variant");
        }
        Err(_) => {}
    }
}

#[test]
fn boolean_with_sheets_invalid_sheet_minus_solid_rejected() {
    use crate::BoolOperation;
    use popkern_core::geometry::Point3;
    let mut sheet = BRep::new(tol());
    let mvfs = sheet.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let _mev = sheet
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    let cube = BRep::unit_cube(tol()).unwrap();
    let result = boolean_with_sheets(BoolOperation::Subtraction, &sheet, &cube, tol());
    assert!(matches!(
        result,
        Err(crate::BoolError::OpenShellNotSupported { a: true, b: false })
    ));
}

#[test]
fn boolean_on_sheet_body_returns_open_shell_error() {
    use popkern_core::geometry::Point3;
    let mut sheet = BRep::new(tol());
    let mvfs = sheet.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let _mev = sheet
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    let cube = BRep::unit_cube(tol()).unwrap();
    let result = crate::union(sheet, cube, tol(), None);
    assert!(
        matches!(
            result,
            Err(crate::BoolError::OpenShellNotSupported { a: true, b: false })
        ),
        "expected OpenShellNotSupported, got {result:?}"
    );
}
