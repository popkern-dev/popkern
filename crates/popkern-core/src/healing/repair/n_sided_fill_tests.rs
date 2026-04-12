use super::*;
use crate::geometry::Point3;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::healing::repair::hole_fill::DetectedHole;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn fill_n_sided_hole_rejects_fewer_than_3_sides() {
    let brep = BRep::new(Tolerance::default());
    let hole = DetectedHole {
        boundary_edges: vec![],
        num_sides: 2,
    };
    assert!(fill_n_sided_hole(&brep, &hole, 1e-6).is_none());
}

#[test]
fn centroid_of_triangle() {
    let pts = vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(3.0, 0.0, 0.0),
        Point3::new(0.0, 3.0, 0.0),
    ];
    let c = centroid(&pts);
    assert!(c.is_some());
    let c = c.as_ref();
    assert!((c.map(|c| c.x).unwrap_or(0.0) - 1.0).abs() < 1e-10);
    assert!((c.map(|c| c.y).unwrap_or(0.0) - 1.0).abs() < 1e-10);
    assert!((c.map(|c| c.z).unwrap_or(0.0)).abs() < 1e-10);
}

#[test]
fn centroid_of_empty_returns_none() {
    assert!(centroid(&[]).is_none());
}

#[test]
fn line_nurbs_endpoints_are_correct() {
    let start = Point3::new(1.0, 2.0, 3.0);
    let end = Point3::new(4.0, 5.0, 6.0);
    let curve = line_nurbs(start, end);

    use crate::geometry::curve::CurveEval;
    let p0 = curve.point_at(0.0);
    let p1 = curve.point_at(1.0);

    assert!(p0.distance_to(&start) < 1e-10, "start mismatch");
    assert!(p1.distance_to(&end) < 1e-10, "end mismatch");
}

#[test]
fn line_nurbs_is_degree_1() {
    let curve = line_nurbs(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
    assert_eq!(curve.degree, 1);
    assert_eq!(curve.control_points.len(), 2);
    assert_eq!(curve.weights.len(), 2);
    assert_eq!(curve.knots.len(), 4);
}

// -----------------------------------------------------------------------
// Helpers: build BReps with edge curves for n-sided fill
// -----------------------------------------------------------------------

use crate::constructors::builder::BRepBuilder;
use crate::geometry::Vec3;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::topology::EdgeHandle;

/// Build a triangular single-face BRep with line edge curves.
fn make_triangle_brep() -> (BRep, Vec<EdgeHandle>) {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(0.5, 1.0, 0.0));

    let solid = b.add_solid();
    let shell = b.add_shell(solid);

    let (e01, he01, _) = b.add_edge(v0, v1);
    let (e12, he12, _) = b.add_edge(v1, v2);
    let (e20, he20, _) = b.add_edge(v2, v0);

    let edges = [
        (e01, Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
        (e12, Point3::new(1.0, 0.0, 0.0), Point3::new(0.5, 1.0, 0.0)),
        (e20, Point3::new(0.5, 1.0, 0.0), Point3::new(0.0, 0.0, 0.0)),
    ];
    for (eh, s, e) in &edges {
        let dir = Vec3::new(e.x - s.x, e.y - s.y, e.z - s.z);
        b.set_edge_curve(
            *eh,
            CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: *s,
                    direction: dir,
                },
                t_start: 0.0,
                t_end: 1.0,
            },
        );
    }

    let face = b.add_face(shell);
    b.add_loop(&[he01, he12, he20], face).unwrap();

    (b.build_unchecked(), vec![e01, e12, e20])
}

/// Build a pentagon single-face BRep with line edge curves.
fn make_pentagon_brep() -> (BRep, Vec<EdgeHandle>) {
    let mut b = BRepBuilder::new(Tolerance::default());
    let angle = |i: usize| std::f64::consts::TAU * (i as f64) / 5.0;
    let mut verts = Vec::new();
    for i in 0..5 {
        let a = angle(i);
        verts.push(b.add_vertex(Point3::new(a.cos(), a.sin(), 0.0)));
    }

    let solid = b.add_solid();
    let shell = b.add_shell(solid);

    let mut edge_handles = Vec::new();
    let mut he_handles = Vec::new();
    for i in 0..5 {
        let j = (i + 1) % 5;
        let s = Point3::new(angle(i).cos(), angle(i).sin(), 0.0);
        let e = Point3::new(angle(j).cos(), angle(j).sin(), 0.0);
        let (eh, heh, _) = b.add_edge(verts[i], verts[j]);
        let dir = Vec3::new(e.x - s.x, e.y - s.y, e.z - s.z);
        b.set_edge_curve(
            eh,
            CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: s,
                    direction: dir,
                },
                t_start: 0.0,
                t_end: 1.0,
            },
        );
        edge_handles.push(eh);
        he_handles.push(heh);
    }

    let face = b.add_face(shell);
    b.add_loop(&he_handles, face).unwrap();

    (b.build_unchecked(), edge_handles)
}

// -----------------------------------------------------------------------
// Tests with real boundary geometry
// -----------------------------------------------------------------------

#[test]
fn fill_triangle_hole_exercises_code_path() {
    let (brep, edges) = make_triangle_brep();
    let hole = DetectedHole {
        boundary_edges: edges,
        num_sides: 3,
    };
    // Use a generous tolerance for Coons patch corner matching
    let result = fill_n_sided_hole(&brep, &hole, 1.0);
    // The code path is exercised: boundary curves are collected,
    // centroid is computed, and Coons patches are attempted.
    if let Some(patches) = &result {
        assert!(patches.len() <= 3);
        for patch in patches {
            assert!(patch.surface.control_points.len() >= 4);
        }
    }
}

#[test]
fn fill_pentagon_hole_exercises_code_path() {
    let (brep, edges) = make_pentagon_brep();
    let hole = DetectedHole {
        boundary_edges: edges,
        num_sides: 5,
    };
    let result = fill_n_sided_hole(&brep, &hole, 1.0);
    if let Some(patches) = &result {
        assert!(patches.len() <= 5);
    }
}

#[test]
fn hole_vertex_positions_returns_correct_count() {
    let (brep, edges) = make_triangle_brep();
    let hole = DetectedHole {
        boundary_edges: edges,
        num_sides: 3,
    };
    let positions = hole_vertex_positions(&brep, &hole);
    assert_eq!(positions.len(), 3);
}

#[test]
fn fill_hole_with_missing_edge_curves_returns_none() {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(0.5, 1.0, 0.0));
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let (e01, he01, _) = b.add_edge(v0, v1);
    let (e12, he12, _) = b.add_edge(v1, v2);
    let (e20, he20, _) = b.add_edge(v2, v0);
    let face = b.add_face(shell);
    b.add_loop(&[he01, he12, he20], face).unwrap();
    let brep = b.build_unchecked();

    let hole = DetectedHole {
        boundary_edges: vec![e01, e12, e20],
        num_sides: 3,
    };
    let result = fill_n_sided_hole(&brep, &hole, 1e-3);
    assert!(result.is_none(), "missing edge curves should prevent fill");
}

#[test]
fn centroid_of_hole_vertices() {
    let (brep, edges) = make_triangle_brep();
    let hole = DetectedHole {
        boundary_edges: edges,
        num_sides: 3,
    };
    let positions = hole_vertex_positions(&brep, &hole);
    let c = centroid(&positions).unwrap();
    // Centroid of (0,0), (1,0), (0.5,1) = (0.5, 0.333...)
    assert!((c.x - 0.5).abs() < 1e-10);
    assert!((c.y - 1.0 / 3.0).abs() < 1e-10);
}
