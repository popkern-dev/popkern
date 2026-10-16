//! Tests for import-specific defect detection.

use super::*;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn valid_cube_no_import_defects() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_all(&brep, &tol, &mut report);
    // The unit cube has planar faces (non-periodic), so no lost
    // periodic topology, and valid orientation.
    assert!(
        report.is_empty(),
        "Valid cube should have no import defects, got {}",
        report.len()
    );
}

#[test]
fn orientation_inconsistency_detected() {
    use crate::geometry::Point3;
    use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face_handle = brep.faces.keys().next().unwrap();
    // Attach a degenerate surface to trigger orientation failure.
    brep.faces[face_handle].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 0.0,
            axis: crate::geometry::Vec3::Z,
        },
        orientation: true,
    });

    let mut report = HealingReport::new();
    check_orientation_inconsistencies(&brep, &mut report);
    assert!(
        !report.is_empty(),
        "Degenerate surface should trigger orientation defect"
    );
    assert_eq!(report.defects()[0].category, DefectCategory::Import);
}

#[test]
fn empty_brep_no_import_defects() {
    let brep = BRep::new(Tolerance::default());
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_all(&brep, &tol, &mut report);
    assert!(report.is_empty());
}

#[test]
fn non_periodic_face_skipped() {
    // A cube has only planar faces; none should trigger periodic checks.
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut report = HealingReport::new();
    check_lost_periodic_topology(&brep, &mut report);
    assert!(report.is_empty());
}

#[test]
fn lost_singularity_topology_on_cube_produces_no_defect() {
    // A cube has planar faces; no sphere/cone surfaces.
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut report = HealingReport::new();
    check_lost_singularity_topology(&brep, &mut report);
    assert!(
        report.is_empty(),
        "Cube with planar faces should have no singularity topology defects"
    );
}

#[test]
fn coordinate_artifacts_on_valid_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_coordinate_artifacts(&brep, &tol, &mut report);
    // Unit cube vertices are at 0.0 and 1.0 — none in the suspicious range.
    assert!(
        report.is_empty(),
        "Valid cube should have no coordinate artifact defects"
    );
}

#[test]
fn empty_brep_no_singularity_or_artifact_defects() {
    let brep = BRep::new(Tolerance::default());
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_lost_singularity_topology(&brep, &mut report);
    check_coordinate_artifacts(&brep, &tol, &mut report);
    assert!(report.is_empty());
}

/// Knot multiplicities on a knot vector with no interior knots.
#[test]
fn knot_multiplicities_no_interior() {
    // Degree 2, knot vector [0,0,0,1,1,1] has no interior knots.
    let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let result = knot_multiplicities(&knots, 2);
    assert!(
        result.is_empty(),
        "Should have no interior knots, got {:?}",
        result
    );
}

/// Knot multiplicities detect a C0 break (multiplicity == degree).
#[test]
fn knot_multiplicities_c0_break() {
    // Degree 3, knot vector [0,0,0,0, 0.5,0.5,0.5, 1,1,1,1]
    // Interior knot 0.5 with multiplicity 3 == degree 3 => C0.
    let knots = vec![0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0];
    let result = knot_multiplicities(&knots, 3);
    assert_eq!(
        result.len(),
        1,
        "Should find one distinct interior knot group"
    );
    assert_eq!(result[0].1, 3, "Multiplicity should be 3");
}

/// Knot multiplicities on a simple cubic with one interior knot.
#[test]
fn knot_multiplicities_single_interior() {
    // Degree 3, knot vector [0,0,0,0, 0.5, 1,1,1,1]
    // Interior knot 0.5 with multiplicity 1 => C2 (no warning).
    let knots = vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0];
    let result = knot_multiplicities(&knots, 3);
    assert_eq!(result.len(), 1);
    assert_eq!(result[0].1, 1, "Multiplicity should be 1");
}

/// Cube with no NURBS geometry produces no representation discontinuities.
#[test]
fn cube_no_representation_discontinuities() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_representation_discontinuities(&brep, &tol, &mut report);
    assert!(
        report.is_empty(),
        "Cube with analytic geometry should have no representation discontinuities"
    );
}

// ── Corrupted BRep tests exercising error-path code ──

#[test]
fn tolerance_mismatch_large_vertex_tolerance() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let vh = brep.vertices.keys().next().unwrap();
    // 1.0 >> default linear 1e-6 => triggers "much larger" branch
    brep.vertices[vh].tolerance = Some(1.0);
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_tolerance_mismatch(&brep, &tol, &mut report);
    assert!(!report.is_empty(), "Should detect large tolerance mismatch");
    assert_eq!(report.defects()[0].category, DefectCategory::Import);
    assert_eq!(report.defects()[0].severity, Severity::Warning);
    assert!(report.defects()[0].measured_value.is_some());
}

#[test]
fn tolerance_mismatch_tiny_edge_tolerance() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let eh = brep.edges.keys().next().unwrap();
    // 1e-15 << default linear 1e-6 * 0.01 => triggers "much tighter" branch
    brep.edges[eh].tolerance = Some(1e-15);
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_tolerance_mismatch(&brep, &tol, &mut report);
    assert!(!report.is_empty(), "Should detect tight tolerance mismatch");
    assert_eq!(report.defects()[0].severity, Severity::Info);
}

#[test]
fn tolerance_mismatch_large_face_tolerance() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let fh = brep.faces.keys().next().unwrap();
    brep.faces[fh].tolerance = Some(1.0);
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_tolerance_mismatch(&brep, &tol, &mut report);
    assert!(!report.is_empty(), "Should detect face tolerance mismatch");
}

#[test]
fn coordinate_artifacts_detected_near_origin() {
    use crate::geometry::Point3;

    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let vh = brep.vertices.keys().next().unwrap();
    // Place vertex in the suspicious range: (0, linear_tol) in all axes
    // Default linear tolerance is 1e-6, so put coords in (1e-15, 1e-6).
    let tiny = 5e-10;
    brep.vertices[vh].position = Point3::new(tiny, tiny, tiny);
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_coordinate_artifacts(&brep, &tol, &mut report);
    assert!(
        !report.is_empty(),
        "Should detect coordinate artifact near origin"
    );
    assert_eq!(report.defects()[0].category, DefectCategory::Import);
    assert_eq!(report.defects()[0].severity, Severity::Info);
    assert!(report.defects()[0].measured_value.is_some());
}

#[test]
fn lost_periodic_topology_on_cylinder_face() {
    use crate::geometry::Point3;
    use crate::geometry::Vec3;
    use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let fh = brep.faces.keys().next().unwrap();
    // Replace with a cylindrical surface (periodic in u).
    // No seam edge exists in the cube topology, so the check will flag it.
    brep.faces[fh].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        orientation: true,
    });
    let mut report = HealingReport::new();
    check_lost_periodic_topology(&brep, &mut report);
    assert!(
        !report.is_empty(),
        "Should detect missing seam edge on periodic surface"
    );
    assert_eq!(report.defects()[0].category, DefectCategory::Import);
    assert!(
        report.defects()[0].description.contains("periodic")
            || report.defects()[0].description.contains("Periodic")
    );
}

#[test]
fn lost_singularity_topology_on_sphere_face() {
    use crate::geometry::Point3;
    use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let fh = brep.faces.keys().next().unwrap();
    // Place a sphere centered at origin with radius 1.
    brep.faces[fh].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
            axis: crate::geometry::Vec3::Z,
        },
        orientation: true,
    });
    // Move a vertex to the south pole (0, 0, -1) so it's "near a pole".
    // This triggers the near_pole && !has_zero_length_edge path.
    let positions = brep.face_outer_loop_positions(fh);
    if !positions.is_empty() {
        // Walk the outer loop to find a vertex we can move.
        let face_data = &brep.faces[fh];
        let loop_data = &brep.loops[face_data.outer_loop];
        let he = &brep.half_edges[loop_data.half_edge];
        brep.vertices[he.origin].position = Point3::new(0.0, 0.0, -1.0);
    }
    let mut report = HealingReport::new();
    check_lost_singularity_topology(&brep, &mut report);
    assert!(
        !report.is_empty(),
        "Should detect missing singularity topology at pole"
    );
    assert!(report.defects()[0].description.contains("pole"));
}

#[test]
fn lost_singularity_topology_on_cone_face() {
    use crate::geometry::Point3;
    use crate::geometry::Vec3;
    use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let fh = brep.faces.keys().next().unwrap();
    // Use a non-vertex-coincident apex so that moving one vertex there
    // doesn't accidentally create a zero-length edge with a neighbor.
    let apex = Point3::new(7.77, 8.88, 9.99);
    brep.faces[fh].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Cone {
            apex,
            axis: Vec3::new(0.0, 0.0, 1.0),
            half_angle: 0.5,
        },
        orientation: true,
    });
    // Move a vertex in this face's outer loop to the apex.
    let face_data = &brep.faces[fh];
    let loop_data = &brep.loops[face_data.outer_loop];
    let he = &brep.half_edges[loop_data.half_edge];
    brep.vertices[he.origin].position = apex;
    let mut report = HealingReport::new();
    check_lost_singularity_topology(&brep, &mut report);
    assert!(
        !report.is_empty(),
        "Should detect missing singularity topology at cone apex"
    );
}

#[test]
fn representation_discontinuity_c0_on_nurbs_edge() {
    use crate::geometry::Point3;
    use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
    use crate::geometry::nurbs_curve::NurbsCurve;

    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let eh = brep.edges.keys().next().unwrap();
    // Degree 3 NURBS with a C0 break (multiplicity 3 at interior knot 0.5).
    // Knots: [0,0,0,0, 0.5,0.5,0.5, 1,1,1,1]  (11 knots)
    // Needs 11 - 4 = 7 control points.
    let nurbs = NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.2, 0.0, 0.0),
            Point3::new(0.4, 0.0, 0.0),
            Point3::new(0.5, 0.0, 0.0),
            Point3::new(0.6, 0.0, 0.0),
            Point3::new(0.8, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0; 7],
    };
    brep.edges[eh].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs),
        t_start: 0.0,
        t_end: 1.0,
    });
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_representation_discontinuities(&brep, &tol, &mut report);
    assert!(
        !report.is_empty(),
        "Should detect C0 discontinuity in NURBS edge"
    );
    assert!(report.defects()[0].description.contains("C0"));
}

#[test]
fn representation_discontinuity_c1_on_nurbs_edge() {
    use crate::geometry::Point3;
    use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
    use crate::geometry::nurbs_curve::NurbsCurve;

    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let eh = brep.edges.keys().next().unwrap();
    // Degree 3, multiplicity 2 at interior knot => C1 discontinuity.
    // Knots: [0,0,0,0, 0.5,0.5, 1,1,1,1]  (10 knots)
    // Needs 10 - 4 = 6 control points.
    let nurbs = NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.2, 0.0, 0.0),
            Point3::new(0.4, 0.0, 0.0),
            Point3::new(0.6, 0.0, 0.0),
            Point3::new(0.8, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0; 6],
    };
    brep.edges[eh].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(nurbs),
        t_start: 0.0,
        t_end: 1.0,
    });
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_representation_discontinuities(&brep, &tol, &mut report);
    assert!(
        !report.is_empty(),
        "Should detect C1 discontinuity in NURBS edge"
    );
    assert!(report.defects()[0].description.contains("C1"));
}
