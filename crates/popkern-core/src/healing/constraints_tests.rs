//! Tests for constraint types, ConstraintSet, and check_constraints.

use super::*;
use crate::geometry::{Point3, Vec3};
use crate::healing::diagnostics::EntityRef;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn empty_constraint_set_no_violations() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let cs = ConstraintSet::new();
    let violations = check_constraints(&brep, &cs, &Tolerance::default());
    assert!(violations.is_empty());
}

#[test]
fn empty_constraint_set_no_violations_on_empty_brep() {
    let brep = BRep::new(Tolerance::default());
    let cs = ConstraintSet::new();
    let violations = check_constraints(&brep, &cs, &Tolerance::default());
    assert!(violations.is_empty());
}

#[test]
fn constraint_set_add_len_is_empty() {
    let mut cs = ConstraintSet::new();
    assert!(cs.is_empty());
    assert_eq!(cs.len(), 0);

    cs.add(GeometricConstraint::FixedPoint {
        vertex: Default::default(),
        position: Point3::new(0.0, 0.0, 0.0),
        max_displacement: 1e-6,
    });
    assert!(!cs.is_empty());
    assert_eq!(cs.len(), 1);

    cs.add(GeometricConstraint::Coplanarity {
        faces: vec![],
        max_deviation: 1e-6,
    });
    assert_eq!(cs.len(), 2);
}

#[test]
fn constraint_set_iter() {
    let mut cs = ConstraintSet::new();
    cs.add(GeometricConstraint::FixedPoint {
        vertex: Default::default(),
        position: Point3::new(1.0, 2.0, 3.0),
        max_displacement: 0.001,
    });
    cs.add(GeometricConstraint::Coplanarity {
        faces: vec![],
        max_deviation: 0.01,
    });

    let count = cs.iter().count();
    assert_eq!(count, 2);
}

#[test]
fn fixed_point_satisfied_on_cube_vertex() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    // Find the first vertex and use its actual position.
    let (vh, v) = brep.vertices.iter().next().unwrap();
    let mut cs = ConstraintSet::new();
    cs.add(GeometricConstraint::FixedPoint {
        vertex: vh,
        position: v.position,
        max_displacement: 1e-6,
    });

    let violations = check_constraints(&brep, &cs, &Tolerance::default());
    assert!(
        violations.is_empty(),
        "constraint should be satisfied when position matches"
    );
}

#[test]
fn fixed_point_violated_after_vertex_move() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    let (vh, v) = brep.vertices.iter().next().unwrap();
    let original = v.position;

    let mut cs = ConstraintSet::new();
    cs.add(GeometricConstraint::FixedPoint {
        vertex: vh,
        position: original,
        max_displacement: 1e-6,
    });

    // Move the vertex far away.
    brep.vertices[vh].position = Point3::new(100.0, 100.0, 100.0);

    let violations = check_constraints(&brep, &cs, &Tolerance::default());
    assert_eq!(violations.len(), 1, "should detect the moved vertex");
    assert_eq!(violations[0].constraint_index, 0);
    assert!(violations[0].violation_magnitude > 0.0);
    assert!(violations[0].description.contains("FixedPoint"));
}

#[test]
fn fixed_point_within_tolerance() {
    let mut brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    let (vh, v) = brep.vertices.iter().next().unwrap();
    let original = v.position;

    let mut cs = ConstraintSet::new();
    cs.add(GeometricConstraint::FixedPoint {
        vertex: vh,
        position: original,
        max_displacement: 0.01,
    });

    // Move vertex slightly within tolerance.
    brep.vertices[vh].position = Point3::new(original.x + 0.005, original.y, original.z);

    let violations = check_constraints(&brep, &cs, &Tolerance::default());
    assert!(violations.is_empty(), "small displacement within tolerance");
}

#[test]
fn coplanarity_on_cube_face() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    // A single cube face is trivially coplanar.
    let fh = brep.faces.iter().next().unwrap().0;
    let mut cs = ConstraintSet::new();
    cs.add(GeometricConstraint::Coplanarity {
        faces: vec![fh],
        max_deviation: 1e-6,
    });

    let violations = check_constraints(&brep, &cs, &Tolerance::default());
    assert!(violations.is_empty());
}

#[test]
fn reflect_point_across_xy_plane() {
    let point = Point3::new(1.0, 2.0, 3.0);
    let plane_pt = Point3::new(0.0, 0.0, 0.0);
    let plane_n = Vec3::new(0.0, 0.0, 1.0);

    let reflected = reflect_point(&point, &plane_pt, &plane_n);
    assert!((reflected.x - 1.0).abs() < 1e-12);
    assert!((reflected.y - 2.0).abs() < 1e-12);
    assert!((reflected.z - (-3.0)).abs() < 1e-12);
}

#[test]
fn reflect_point_across_offset_plane() {
    // Reflect (1, 2, 5) across z=2 plane.
    let point = Point3::new(1.0, 2.0, 5.0);
    let plane_pt = Point3::new(0.0, 0.0, 2.0);
    let plane_n = Vec3::new(0.0, 0.0, 1.0);

    let reflected = reflect_point(&point, &plane_pt, &plane_n);
    assert!((reflected.x - 1.0).abs() < 1e-12);
    assert!((reflected.y - 2.0).abs() < 1e-12);
    assert!((reflected.z - (-1.0)).abs() < 1e-12);
}

#[test]
fn best_fit_plane_on_coplanar_points() {
    let points = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(1.0, 1.0, 0.0),
    ];

    let (centroid, normal) = best_fit_plane(&points).unwrap();
    // Centroid should be at (0.5, 0.5, 0.0).
    assert!((centroid.x - 0.5).abs() < 1e-10);
    assert!((centroid.y - 0.5).abs() < 1e-10);
    assert!((centroid.z).abs() < 1e-10);
    // Normal should be (0, 0, +/-1).
    assert!(normal.z.abs() > 0.99);
}

#[test]
fn best_fit_plane_too_few_points() {
    let points = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)];
    assert!(best_fit_plane(&points).is_none());
}

#[test]
fn distance_constraint_satisfied() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    let mut verts = brep.vertices.iter();
    let (v1, _) = verts.next().unwrap();
    let (v2, _) = verts.next().unwrap();

    let actual_dist = brep.vertices[v1]
        .position
        .distance_to(&brep.vertices[v2].position);

    let mut cs = ConstraintSet::new();
    cs.add(GeometricConstraint::Distance {
        entity_a: EntityRef::Vertex(v1),
        entity_b: EntityRef::Vertex(v2),
        target_distance: actual_dist,
        tolerance: 0.01,
    });

    let violations = check_constraints(&brep, &cs, &Tolerance::default());
    assert!(violations.is_empty());
}

#[test]
fn distance_constraint_violated() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    let mut verts = brep.vertices.iter();
    let (v1, _) = verts.next().unwrap();
    let (v2, _) = verts.next().unwrap();

    let mut cs = ConstraintSet::new();
    cs.add(GeometricConstraint::Distance {
        entity_a: EntityRef::Vertex(v1),
        entity_b: EntityRef::Vertex(v2),
        target_distance: 999.0, // Way off.
        tolerance: 0.001,
    });

    let violations = check_constraints(&brep, &cs, &Tolerance::default());
    assert_eq!(violations.len(), 1);
    assert!(violations[0].description.contains("Distance"));
}

#[test]
fn entity_representative_point_vertex() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let (vh, v) = brep.vertices.iter().next().unwrap();
    let pt = entity_representative_point(&brep, &EntityRef::Vertex(vh));
    assert!(pt.is_some());
    let pt = pt.unwrap();
    assert!((pt.x - v.position.x).abs() < 1e-12);
}

#[test]
fn entity_representative_point_face() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let (fh, _) = brep.faces.iter().next().unwrap();
    let pt = entity_representative_point(&brep, &EntityRef::Face(fh));
    assert!(pt.is_some());
}

#[test]
fn entity_representative_point_edge() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);
    let (eh, _) = brep.edges.iter().next().unwrap();
    let pt = entity_representative_point(&brep, &EntityRef::Edge(eh));
    assert!(pt.is_some());
}

#[test]
fn entity_representative_point_invalid_handle() {
    let brep = BRep::new(Tolerance::default());
    let pt = entity_representative_point(&brep, &EntityRef::Vertex(Default::default()));
    assert!(pt.is_none());
}

#[test]
fn tangency_constraint_on_adjacent_cube_faces() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    let faces: Vec<_> = brep.faces.iter().map(|(fh, _)| fh).collect();
    if faces.len() >= 2 {
        // Cube faces are perpendicular (90 degrees), so a tight tangency
        // constraint should fail.
        let mut cs = ConstraintSet::new();
        cs.add(GeometricConstraint::Tangency {
            face_a: faces[0],
            face_b: faces[1],
            max_angle: 0.01, // Very tight — cube faces are 90 degrees apart.
        });

        let violations = check_constraints(&brep, &cs, &Tolerance::default());
        // Cube faces meet at 90 degrees; tangency at 0.01 rad should fail.
        assert!(
            !violations.is_empty(),
            "perpendicular faces should violate tight tangency"
        );
    }
}

#[test]
fn multiple_constraints_mixed_results() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    let (vh, v) = brep.vertices.iter().next().unwrap();

    let mut cs = ConstraintSet::new();

    // This one should pass.
    cs.add(GeometricConstraint::FixedPoint {
        vertex: vh,
        position: v.position,
        max_displacement: 1e-6,
    });

    // This one should fail (target distance is impossible).
    let mut verts = brep.vertices.iter();
    let (v1, _) = verts.next().unwrap();
    let (v2, _) = verts.next().unwrap();
    cs.add(GeometricConstraint::Distance {
        entity_a: EntityRef::Vertex(v1),
        entity_b: EntityRef::Vertex(v2),
        target_distance: 999.0,
        tolerance: 0.001,
    });

    let violations = check_constraints(&brep, &cs, &Tolerance::default());
    assert_eq!(violations.len(), 1);
    assert_eq!(violations[0].constraint_index, 1);
}

#[test]
fn constraint_violation_has_description() {
    let brep = BRep::unit_cube(Tolerance::default())
        .map_or_else(|_| BRep::new(Tolerance::default()), |b| b);

    let (vh, _) = brep.vertices.iter().next().unwrap();
    let mut cs = ConstraintSet::new();
    cs.add(GeometricConstraint::FixedPoint {
        vertex: vh,
        position: Point3::new(999.0, 999.0, 999.0),
        max_displacement: 0.001,
    });

    let violations = check_constraints(&brep, &cs, &Tolerance::default());
    assert_eq!(violations.len(), 1);
    assert!(!violations[0].description.is_empty());
}
