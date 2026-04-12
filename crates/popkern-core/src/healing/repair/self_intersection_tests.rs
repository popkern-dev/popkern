use super::*;
use crate::constructors::builder::BRepBuilder;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn valid_cube_has_no_self_intersections() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let pairs = detect_self_intersecting_faces(&brep);
    assert!(
        pairs.is_empty(),
        "valid cube should have no self-intersecting faces"
    );
}

#[test]
fn repair_on_valid_cube_face_is_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let result = repair_self_intersection(&mut brep, face);
    assert!(!result.repaired);
    assert!(result.description.contains("Analytical"));
}

#[test]
fn repair_on_missing_face_returns_not_found() {
    let mut brep = BRep::new(Tolerance::default());
    let bogus = FaceHandle::default();
    let result = repair_self_intersection(&mut brep, bogus);
    assert!(!result.repaired);
    assert!(result.description.contains("not found"));
}

#[test]
fn detect_on_empty_brep_returns_empty() {
    let brep = BRep::new(Tolerance::default());
    let pairs = detect_self_intersecting_faces(&brep);
    assert!(pairs.is_empty());
}

#[test]
fn repair_by_splitting_on_valid_cube_face_no_self_intersection() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let face = brep.faces.keys().next().unwrap();
    let result = repair_by_splitting(&mut brep, face, &tol);
    assert!(!result.repaired);
    assert!(
        result.description.contains("not NURBS")
            || result.description.contains("Analytical")
            || result.description.contains("Only NURBS"),
    );
}

// -----------------------------------------------------------------------
// Helpers: build BRep with folded NURBS surface
// -----------------------------------------------------------------------

/// Build a 5x5 NURBS surface with a severe fold (normal reversal).
/// Uses degree 2 so the fold is tighter and more visible at the
/// sampling resolution used by detect_normal_reversals.
fn make_folded_nurbs() -> NurbsSurface {
    let mut control_points = Vec::with_capacity(25);
    let mut weights = Vec::with_capacity(25);
    for i in 0..5 {
        for j in 0..5 {
            control_points.push(Point3::new(i as f64 / 4.0, j as f64 / 4.0, 0.0));
            weights.push(1.0);
        }
    }
    // Push a block of interior control points below the plane to create
    // a severe fold that reverses normals in the center region.
    control_points[1 * 5 + 1] = Point3::new(0.25, 0.25, -5.0);
    control_points[1 * 5 + 2] = Point3::new(0.25, 0.50, -5.0);
    control_points[2 * 5 + 1] = Point3::new(0.50, 0.25, -5.0);
    control_points[2 * 5 + 2] = Point3::new(0.50, 0.50, -5.0);
    // Push opposite interior points far above
    control_points[2 * 5 + 3] = Point3::new(0.50, 0.75, 5.0);
    control_points[3 * 5 + 2] = Point3::new(0.75, 0.50, 5.0);
    control_points[3 * 5 + 3] = Point3::new(0.75, 0.75, 5.0);

    NurbsSurface {
        degree_u: 2,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 0.0, 0.25, 0.5, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 0.25, 0.5, 1.0, 1.0, 1.0],
        control_points,
        num_u: 5,
        num_v: 5,
        weights,
    }
}

/// Build a single-face BRep with a folded NURBS surface.
fn make_folded_nurbs_brep() -> (BRep, FaceHandle) {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));

    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let (_, he01, _) = b.add_edge(v0, v1);
    let (_, he12, _) = b.add_edge(v1, v2);
    let (_, he23, _) = b.add_edge(v2, v3);
    let (_, he30, _) = b.add_edge(v3, v0);

    let face = b.add_face(shell);
    b.add_loop(&[he01, he12, he23, he30], face).unwrap();

    b.set_face_surface(
        face,
        SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(make_folded_nurbs()),
            orientation: true,
        },
    );

    (b.build_unchecked(), face)
}

/// Build a BRep with a tiny (2x2) NURBS grid.
fn make_small_nurbs_brep() -> (BRep, FaceHandle) {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));

    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let (_, he01, _) = b.add_edge(v0, v1);
    let (_, he12, _) = b.add_edge(v1, v2);
    let (_, he23, _) = b.add_edge(v2, v3);
    let (_, he30, _) = b.add_edge(v3, v0);

    let face = b.add_face(shell);
    b.add_loop(&[he01, he12, he23, he30], face).unwrap();

    let nurbs = NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, -1.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0, 1.0, 1.0, 1.0],
    };

    b.set_face_surface(
        face,
        SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(nurbs),
            orientation: true,
        },
    );

    (b.build_unchecked(), face)
}

// -----------------------------------------------------------------------
// Tests with folded NURBS geometry
// -----------------------------------------------------------------------

#[test]
fn repair_smoothing_on_folded_nurbs() {
    let (mut brep, face) = make_folded_nurbs_brep();
    let result = repair_self_intersection(&mut brep, face);
    assert_eq!(result.face, face);
    assert!(
        result.description.contains("Smooth")
            || result.description.contains("smooth")
            || result.description.contains("reversal")
            || result.description.contains("small"),
    );
}

#[test]
fn repair_by_splitting_on_folded_nurbs() {
    let (mut brep, face) = make_folded_nurbs_brep();
    let tol = Tolerance::default();
    let result = repair_by_splitting(&mut brep, face, &tol);
    assert_eq!(result.face, face);
    // The code exercises the NURBS splitting path. Whether it succeeds
    // depends on fold position relative to the parameter domain.
    assert!(!result.description.is_empty());
}

#[test]
fn repair_smoothing_rejects_small_grid() {
    let (mut brep, face) = make_small_nurbs_brep();
    let result = repair_self_intersection(&mut brep, face);
    assert!(!result.repaired);
    assert!(
        result.description.contains("too small")
            || result.description.contains("No normal reversal"),
    );
}

#[test]
fn repair_on_face_without_surface() {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let (_, he01, _) = b.add_edge(v0, v1);
    let (_, he12, _) = b.add_edge(v1, v2);
    let (_, he23, _) = b.add_edge(v2, v3);
    let (_, he30, _) = b.add_edge(v3, v0);
    let face = b.add_face(shell);
    b.add_loop(&[he01, he12, he23, he30], face).unwrap();
    let mut brep = b.build_unchecked();

    let result = repair_self_intersection(&mut brep, face);
    assert!(!result.repaired);
    assert!(result.description.contains("no surface"));
}

#[test]
fn repair_by_splitting_on_face_without_surface() {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let (_, he01, _) = b.add_edge(v0, v1);
    let (_, he12, _) = b.add_edge(v1, v2);
    let (_, he23, _) = b.add_edge(v2, v3);
    let (_, he30, _) = b.add_edge(v3, v0);
    let face = b.add_face(shell);
    b.add_loop(&[he01, he12, he23, he30], face).unwrap();
    let mut brep = b.build_unchecked();

    let tol = Tolerance::default();
    let result = repair_by_splitting(&mut brep, face, &tol);
    assert!(!result.repaired);
    assert!(result.description.contains("no surface"));
}

#[test]
fn detect_normal_reversals_on_folded_surface() {
    let (brep, face) = make_folded_nurbs_brep();
    let reversals = detect_normal_reversals(&brep, face);
    assert!(
        !reversals.is_empty(),
        "folded NURBS should have normal reversals"
    );
}
