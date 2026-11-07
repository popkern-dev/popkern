//! Wire ∩ / − / ∪ solid against analytic lengths.

use std::f64::consts::TAU;

use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

use super::boolean_with_wire;
use crate::error::{BoolError, BoolOperation};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// The polyline through `points` as a wire body.
fn polyline(points: &[Point3]) -> BRep {
    let mut b = BRepBuilder::new(tol());
    let vs: Vec<_> = points.iter().map(|&p| b.add_vertex(p)).collect();
    for (k, w) in vs.windows(2).enumerate() {
        let (eh, _, _) = b.add_edge(w[0], w[1]);
        let segment = CurveOnEdge::line_segment(points[k], points[k + 1], 1e-9).expect("segment");
        b.set_edge_curve(eh, segment);
    }
    b.build_unchecked()
}

/// The full circle of `radius` about `center` in a z plane, one closed
/// edge.
fn circle(center: Point3, radius: f64) -> BRep {
    let geometry = CurveGeometry::Circle {
        center,
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius,
    };
    let mut b = BRepBuilder::new(tol());
    let v = b.add_vertex(geometry.point_at(0.0));
    let (eh, _, _) = b.add_edge(v, v);
    b.set_edge_curve(
        eh,
        CurveOnEdge {
            geometry,
            t_start: 0.0,
            t_end: TAU,
        },
    );
    b.build_unchecked()
}

/// The box `[0, 2]³`.
fn cube() -> BRep {
    let square = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 2.0, 0.0),
        Point3::new(0.0, 2.0, 0.0),
    ];
    BRep::prism(&square, Vec3::new(0.0, 0.0, 1.0), 2.0, tol()).expect("cube")
}

/// Total length of a wire body's edges.
fn length(wire: &BRep) -> f64 {
    wire.edges
        .values()
        .map(|e| {
            let c = e.curve.as_ref().expect("curve");
            match &c.geometry {
                CurveGeometry::Circle { radius, .. } => (c.t_end - c.t_start).abs() * radius,
                g => g.point_at(c.t_start).distance_to(&g.point_at(c.t_end)),
            }
        })
        .sum()
}

fn run(op: BoolOperation, a: &BRep, b: &BRep) -> BRep {
    boolean_with_wire(op, a, b, tol()).expect("wire boolean")
}

#[test]
fn a_segment_through_a_box_keeps_its_inside_and_outside_parts() {
    let wire = polyline(&[Point3::new(-1.0, 1.0, 1.0), Point3::new(3.0, 1.0, 1.0)]);
    let inside = run(BoolOperation::Intersection, &wire, &cube());
    assert_eq!(inside.body_kind(), BodyKind::Wire);
    assert_eq!(inside.edges.len(), 1);
    assert!((length(&inside) - 2.0).abs() < 1e-9, "{}", length(&inside));
    // Either order gives the intersection.
    let flipped = run(BoolOperation::Intersection, &cube(), &wire);
    assert!((length(&flipped) - 2.0).abs() < 1e-9);
    let outside = run(BoolOperation::Subtraction, &wire, &cube());
    assert_eq!(outside.edges.len(), 2);
    assert!(
        (length(&outside) - 2.0).abs() < 1e-9,
        "{}",
        length(&outside)
    );
    // A wire removes no volume.
    let solid = run(BoolOperation::Subtraction, &cube(), &wire);
    assert_eq!(solid.body_kind(), BodyKind::Solid);
    assert_eq!(solid.faces.len(), cube().faces.len());
    // The union holds the solid and the outside pieces.
    let union = run(BoolOperation::Union, &cube(), &wire);
    assert_eq!(union.body_kind(), BodyKind::General);
    assert_eq!(union.faces.len(), cube().faces.len());
    assert_eq!(union.edges.len(), cube().edges.len() + 2);
}

#[test]
fn a_bent_wire_stays_joined_at_its_kept_vertex() {
    let wire = polyline(&[
        Point3::new(-1.0, 1.0, 1.0),
        Point3::new(1.0, 1.0, 1.0),
        Point3::new(1.0, 1.0, 3.0),
    ]);
    let inside = run(BoolOperation::Intersection, &wire, &cube());
    assert_eq!((inside.edges.len(), inside.vertices.len()), (2, 3));
    assert!((length(&inside) - 2.0).abs() < 1e-9);
    let outside = run(BoolOperation::Subtraction, &wire, &cube());
    assert_eq!((outside.edges.len(), outside.vertices.len()), (2, 4));
    assert!((length(&outside) - 2.0).abs() < 1e-9);
}

#[test]
fn a_wire_on_the_boundary_belongs_to_the_intersection() {
    let wire = polyline(&[Point3::new(-1.0, 1.0, 2.0), Point3::new(3.0, 1.0, 2.0)]);
    let inside = run(BoolOperation::Intersection, &wire, &cube());
    assert!((length(&inside) - 2.0).abs() < 1e-9, "{}", length(&inside));
    let outside = run(BoolOperation::Subtraction, &wire, &cube());
    assert!(
        (length(&outside) - 2.0).abs() < 1e-9,
        "{}",
        length(&outside)
    );
}

/// A radius-1.5 circle about (1, 0) crosses a radius-1 cylinder about
/// the axis where cos θ = −0.75: the arc inside spans the angles with
/// cos θ below it.
#[test]
fn a_circle_crossing_a_cylinder_keeps_the_analytic_arc() {
    let cylinder = BRep::cylinder(1.0, 4.0, tol()).expect("cylinder");
    let z = cylinder
        .vertices
        .values()
        .map(|v| v.position.z)
        .fold(f64::INFINITY, f64::min)
        + 1.0;
    let wire = circle(Point3::new(1.0, 0.0, z), 1.5);
    let inside = run(BoolOperation::Intersection, &wire, &cylinder);
    let arc = 1.5 * (TAU - 2.0 * (-0.75f64).acos());
    assert!(
        (length(&inside) - arc).abs() < 1e-7,
        "{} vs {arc}",
        length(&inside)
    );
    let outside = run(BoolOperation::Subtraction, &wire, &cylinder);
    assert!((length(&outside) - (1.5 * TAU - arc)).abs() < 1e-7);
    // A circle wholly outside is kept whole, on its one vertex.
    let apart = circle(Point3::new(0.0, 0.0, z), 2.0);
    let kept = run(BoolOperation::Subtraction, &apart, &cylinder);
    assert_eq!((kept.edges.len(), kept.vertices.len()), (1, 1));
    assert!((length(&kept) - 2.0 * TAU).abs() < 1e-9);
    assert!(
        run(BoolOperation::Intersection, &apart, &cylinder)
            .edges
            .is_empty()
    );
}

#[test]
fn only_a_wire_and_a_solid_are_accepted() {
    let wire = polyline(&[Point3::new(-1.0, 1.0, 1.0), Point3::new(3.0, 1.0, 1.0)]);
    for (a, b) in [(&cube(), &cube()), (&wire, &wire)] {
        assert!(matches!(
            boolean_with_wire(BoolOperation::Intersection, a, b, tol()),
            Err(BoolError::UnsupportedBodyKind { .. })
        ));
    }
}
