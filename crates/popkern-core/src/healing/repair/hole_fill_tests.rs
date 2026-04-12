use super::*;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn valid_cube_has_no_holes() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let holes = detect_holes(&brep);
    assert!(holes.is_empty(), "valid cube should have no boundary holes");
}

#[test]
fn valid_cube_has_no_free_edges() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let free = find_free_edges(&brep);
    assert!(
        free.is_empty(),
        "valid cube should have no free boundary edges"
    );
}

#[test]
fn fill_holes_on_valid_cube_is_noop() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let result = fill_holes(&mut brep, HoleFillStrategy::Extension, &tol);
    assert_eq!(result.holes_detected, 0);
    assert_eq!(result.holes_filled, 0);
    assert_eq!(result.strategy, HoleFillStrategy::Extension);
}

#[test]
fn empty_brep_has_no_holes() {
    let brep = BRep::new(Tolerance::default());
    let holes = detect_holes(&brep);
    assert!(holes.is_empty());
}

#[test]
fn fill_by_extension_on_valid_cube_returns_none() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let holes = detect_holes(&brep);
    assert!(holes.is_empty());
    let fake_hole = DetectedHole {
        boundary_edges: Vec::new(),
        num_sides: 0,
    };
    let result = fill_by_extension(&brep, &fake_hole, &tol);
    assert!(
        result.is_none(),
        "extension should return None for empty hole boundary"
    );
}

#[test]
fn auto_select_small_hole_prefers_extension() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let hole = DetectedHole {
        boundary_edges: Vec::new(),
        num_sides: 3,
    };
    assert_eq!(
        auto_select_strategy(&brep, &hole, &tol),
        HoleFillStrategy::Extension
    );
}

#[test]
fn auto_select_large_hole_prefers_energy_minimizing() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let hole = DetectedHole {
        boundary_edges: Vec::new(),
        num_sides: 5,
    };
    assert_eq!(
        auto_select_strategy(&brep, &hole, &tol),
        HoleFillStrategy::EnergyMinimizing
    );
}

#[test]
fn auto_select_4sided_hole_without_extension() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let hole = DetectedHole {
        boundary_edges: Vec::new(),
        num_sides: 4,
    };
    assert_eq!(
        auto_select_strategy(&brep, &hole, &tol),
        HoleFillStrategy::EnergyMinimizing
    );
}

#[test]
fn energy_minimizing_reports_detection_only() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let result = fill_holes(&mut brep, HoleFillStrategy::EnergyMinimizing, &tol);
    assert_eq!(result.holes_detected, 0);
    assert_eq!(result.holes_filled, 0);
    assert_eq!(result.strategy, HoleFillStrategy::EnergyMinimizing);
}

// -----------------------------------------------------------------------
// Tests with defective geometry (single-face BRep = boundary hole)
// -----------------------------------------------------------------------

use crate::constructors::builder::BRepBuilder;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};

/// Build a single-face BRep where all 4 edges are free boundaries.
/// Manually assigns twin half-edges to the face's loop so that
/// `find_free_edges` (which checks face_a == face_b) detects them.
fn make_single_face_brep() -> BRep {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let (_, he01, tw01) = b.add_edge(v0, v1);
    let (_, he12, tw12) = b.add_edge(v1, v2);
    let (_, he23, tw23) = b.add_edge(v2, v3);
    let (_, he30, tw30) = b.add_edge(v3, v0);
    let face = b.add_face(shell);
    let lh = b.add_loop(&[he01, he12, he23, he30], face).unwrap();
    // Assign twins to the same loop so they are recognized as free edges
    let brep = b.brep_mut();
    brep.half_edges[tw01].loop_ = lh;
    brep.half_edges[tw12].loop_ = lh;
    brep.half_edges[tw23].loop_ = lh;
    brep.half_edges[tw30].loop_ = lh;
    b.build_unchecked()
}

/// Build a single-face BRep with a planar surface and line edge curves.
fn make_single_face_with_geometry() -> BRep {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let solid = b.add_solid();
    let shell = b.add_shell(solid);

    let (e01, he01, tw01) = b.add_edge(v0, v1);
    let (e12, he12, tw12) = b.add_edge(v1, v2);
    let (e23, he23, tw23) = b.add_edge(v2, v3);
    let (e30, he30, tw30) = b.add_edge(v3, v0);

    let edges_verts = [
        (e01, Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
        (e12, Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        (e23, Point3::new(1.0, 1.0, 0.0), Point3::new(0.0, 1.0, 0.0)),
        (e30, Point3::new(0.0, 1.0, 0.0), Point3::new(0.0, 0.0, 0.0)),
    ];
    for (eh, s, e) in &edges_verts {
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
    let lh = b.add_loop(&[he01, he12, he23, he30], face).unwrap();
    let brep = b.brep_mut();
    brep.half_edges[tw01].loop_ = lh;
    brep.half_edges[tw12].loop_ = lh;
    brep.half_edges[tw23].loop_ = lh;
    brep.half_edges[tw30].loop_ = lh;

    b.set_face_surface(
        face,
        SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
            orientation: true,
        },
    );

    b.build_unchecked()
}

#[test]
fn single_face_has_free_edges() {
    let brep = make_single_face_brep();
    let free = find_free_edges(&brep);
    assert_eq!(free.len(), 4);
}

#[test]
fn detect_holes_on_single_face() {
    let brep = make_single_face_brep();
    let holes = detect_holes(&brep);
    assert!(!holes.is_empty(), "single face boundary should be a hole");
    assert_eq!(holes[0].num_sides, 4);
}

#[test]
fn fill_holes_on_single_face_detects_hole() {
    let mut brep = make_single_face_brep();
    let tol = Tolerance::default();
    let result = fill_holes(&mut brep, HoleFillStrategy::Extension, &tol);
    assert!(result.holes_detected > 0);
}

#[test]
fn fill_by_extension_with_surface_geometry() {
    let brep = make_single_face_with_geometry();
    let tol = Tolerance::default();
    let holes = detect_holes(&brep);
    assert!(!holes.is_empty());
    let _result = fill_by_extension(&brep, &holes[0], &tol);
}

#[test]
fn auto_select_4sided_hole_with_real_edges() {
    let brep = make_single_face_with_geometry();
    let tol = Tolerance::default();
    let holes = detect_holes(&brep);
    assert!(!holes.is_empty());
    assert_eq!(holes[0].num_sides, 4);
    let strategy = auto_select_strategy(&brep, &holes[0], &tol);
    assert!(
        strategy == HoleFillStrategy::Extension || strategy == HoleFillStrategy::EnergyMinimizing,
    );
}

#[test]
fn build_vertex_edge_adjacency_for_single_face() {
    let brep = make_single_face_brep();
    let free = find_free_edges(&brep);
    let adjacency = build_vertex_edge_adjacency(&brep, &free);
    assert_eq!(adjacency.len(), 4);
    for (_vh, edges) in &adjacency {
        assert_eq!(edges.len(), 2);
    }
}

#[test]
fn sample_edge_points_without_curve_uses_interpolation() {
    let brep = make_single_face_brep();
    let edge = brep.edges.keys().next().unwrap();
    let pts = sample_edge_points(&brep, edge, 4);
    assert_eq!(pts.len(), 4);
}

#[test]
fn sample_edge_points_with_curve() {
    let brep = make_single_face_with_geometry();
    let edge = brep.edges.keys().next().unwrap();
    let pts = sample_edge_points(&brep, edge, 6);
    assert_eq!(pts.len(), 6);
}

#[test]
fn fill_holes_energy_minimizing_on_defective_brep() {
    let mut brep = make_single_face_brep();
    let tol = Tolerance::default();
    let result = fill_holes(&mut brep, HoleFillStrategy::EnergyMinimizing, &tol);
    assert!(result.holes_detected > 0);
    assert_eq!(result.holes_filled, 0);
}
