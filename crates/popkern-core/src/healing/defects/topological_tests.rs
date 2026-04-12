use super::*;
use crate::geometry::Point3;
use crate::healing::diagnostics::Severity;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn valid_cube_produces_no_topological_defects() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_all(&brep, &tol, &mut report);
    assert!(
        report.is_empty(),
        "Valid cube should have no topological defects, got {}",
        report.len()
    );
}

#[test]
fn broken_twin_produces_defect() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let he_handle = brep.half_edges.keys().next().unwrap();
    brep.half_edges[he_handle].twin = Default::default();

    let mut report = HealingReport::new();
    check_half_edge_pairing(&brep, &mut report);
    assert!(!report.is_empty(), "Broken twin should produce a defect");
    assert_eq!(report.defects()[0].category, DefectCategory::Topological);
    assert_eq!(report.defects()[0].severity, Severity::Error);
}

#[test]
fn corrupted_vertex_produces_defect() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let v_handle = brep.vertices.keys().next().unwrap();
    brep.vertices[v_handle].position = Point3::new(999.0, 999.0, 999.0);

    let mut report = HealingReport::new();
    check_shell_closure(&brep, &mut report);
    assert!(
        !report.is_empty(),
        "Corrupted vertex should trigger a shell closure defect"
    );
}

#[test]
fn orientation_defect_is_warning() {
    use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face_handle = brep.faces.keys().next().unwrap();
    brep.faces[face_handle].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 0.0,
        },
        orientation: true,
    });

    let mut report = HealingReport::new();
    check_face_orientation(&brep, &mut report);
    assert!(!report.is_empty());
    assert_eq!(report.defects()[0].severity, Severity::Warning);
}

#[test]
fn extract_entities_covers_all_variants() {
    // BrokenTwin
    let e = ValidationError::BrokenTwin {
        half_edge: Default::default(),
    };
    let ents = extract_entities(&e);
    assert_eq!(ents.len(), 1);
    assert!(matches!(ents[0], EntityRef::HalfEdge(_)));

    // NonManifoldEdge
    let e = ValidationError::NonManifoldEdge {
        edge: Default::default(),
    };
    let ents = extract_entities(&e);
    assert!(matches!(ents[0], EntityRef::Edge(_)));

    // InconsistentOrientation
    let e = ValidationError::InconsistentOrientation {
        face: Default::default(),
    };
    let ents = extract_entities(&e);
    assert!(matches!(ents[0], EntityRef::Face(_)));

    // EulerPoincare (no entities)
    let e = ValidationError::EulerPoincare {
        solid: Default::default(),
        expected: 0,
        actual: 2,
    };
    let ents = extract_entities(&e);
    assert!(ents.is_empty());
}

#[test]
fn measurements_extracted_from_gap_error() {
    let e = ValidationError::ShellGap {
        edge: Default::default(),
        gap: 0.05,
    };
    let (measured, _threshold) = extract_measurements(&e);
    assert_eq!(measured, Some(0.05));
}

#[test]
fn empty_brep_no_defects() {
    let brep = BRep::new(Tolerance::default());
    let mut report = HealingReport::new();
    check_all(&brep, &Tolerance::default(), &mut report);
    assert!(report.is_empty());
}

#[test]
fn boundary_simplicity_on_valid_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_boundary_simplicity(&brep, &tol, &mut report);
    assert!(
        report.is_empty(),
        "Valid cube loops should not self-intersect, got {}",
        report.len()
    );
}

#[test]
fn boundary_simplicity_on_empty_brep() {
    let brep = BRep::new(Tolerance::default());
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_boundary_simplicity(&brep, &tol, &mut report);
    assert!(report.is_empty());
}

#[test]
fn inner_loop_orientation_on_valid_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut report = HealingReport::new();
    check_inner_loop_orientation(&brep, &mut report);
    // A cube has no inner loops, so no defects.
    assert!(
        report.is_empty(),
        "Cube with no inner loops should have no inner loop orientation defects"
    );
}

#[test]
fn inner_loop_orientation_on_empty_brep() {
    let brep = BRep::new(Tolerance::default());
    let mut report = HealingReport::new();
    check_inner_loop_orientation(&brep, &mut report);
    assert!(report.is_empty());
}

// ── Corrupted BRep tests exercising defect-detection error paths ──

#[test]
fn corrupted_vertex_triggers_closure_defect() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Find a vertex at the start of a loop -- the loop's starting half-edge
    // determines first_pos. Move that vertex far away so the closure gap
    // (distance from last curve endpoint to first vertex) is large.
    let loop_h = brep.loops.keys().next().unwrap();
    let start_he = brep.loops[loop_h].half_edge;
    let start_v = brep.half_edges[start_he].origin;
    brep.vertices[start_v].position = Point3::new(100.0, 100.0, 100.0);
    let mut report = HealingReport::new();
    check_boundary_closure(&brep, &Tolerance::default(), &mut report);
    assert!(!report.is_empty(), "Should detect loop closure failure");
    assert_eq!(report.defects()[0].category, DefectCategory::Topological);
    assert_eq!(report.defects()[0].severity, Severity::Error);
    // OpenLoop variant should carry a measured gap value.
    assert!(report.defects()[0].measured_value.is_some());
}

#[test]
fn corrupted_cycle_triggers_cycle_defect() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let he = brep.half_edges.keys().next().unwrap();
    brep.half_edges[he].next = he; // self-loop breaks cycle
    let mut report = HealingReport::new();
    check_half_edge_cycles(&brep, &mut report);
    assert!(!report.is_empty(), "Should detect broken half-edge cycle");
    assert_eq!(report.defects()[0].category, DefectCategory::Topological);
}

#[test]
fn corrupted_outgoing_triggers_vertex_consistency_defect() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Find two distinct vertices and cross-wire the outgoing pointers.
    let mut verts = brep.vertices.keys();
    let v1 = verts.next().unwrap();
    let v2 = verts.next().unwrap();
    // Point v1's outgoing to v2's outgoing half-edge (wrong origin).
    let v2_outgoing = brep.vertices[v2].outgoing;
    brep.vertices[v1].outgoing = v2_outgoing;
    let mut report = HealingReport::new();
    check_vertex_consistency(&brep, &mut report);
    assert!(
        !report.is_empty(),
        "Should detect vertex consistency failure"
    );
}

#[test]
fn extract_entities_covers_remaining_variants() {
    // BrokenCycle
    let e = ValidationError::BrokenCycle {
        half_edge: Default::default(),
        loop_: Default::default(),
    };
    let ents = extract_entities(&e);
    assert_eq!(ents.len(), 2);
    assert!(matches!(ents[0], EntityRef::HalfEdge(_)));
    assert!(matches!(ents[1], EntityRef::Loop(_)));

    // InconsistentVertex
    let e = ValidationError::InconsistentVertex {
        vertex: Default::default(),
        half_edge: Default::default(),
    };
    let ents = extract_entities(&e);
    assert_eq!(ents.len(), 2);
    assert!(matches!(ents[0], EntityRef::Vertex(_)));

    // LowValence
    let e = ValidationError::LowValence {
        vertex: Default::default(),
        valence: 1,
    };
    let ents = extract_entities(&e);
    assert_eq!(ents.len(), 1);
    assert!(matches!(ents[0], EntityRef::Vertex(_)));

    // OpenLoop
    let e = ValidationError::OpenLoop {
        loop_: Default::default(),
        gap: 0.1,
    };
    let ents = extract_entities(&e);
    assert_eq!(ents.len(), 1);
    assert!(matches!(ents[0], EntityRef::Loop(_)));

    // CurveSurfaceDeviation
    let e = ValidationError::CurveSurfaceDeviation {
        edge: Default::default(),
        max_deviation: 0.5,
    };
    let ents = extract_entities(&e);
    assert_eq!(ents.len(), 1);
    assert!(matches!(ents[0], EntityRef::Edge(_)));

    // SelfIntersection
    let e = ValidationError::SelfIntersection {
        face_a: Default::default(),
        face_b: Default::default(),
    };
    let ents = extract_entities(&e);
    assert_eq!(ents.len(), 1);
    assert!(matches!(ents[0], EntityRef::FacePair(_, _)));

    // TwistedSurface
    let e = ValidationError::TwistedSurface {
        face: Default::default(),
    };
    let ents = extract_entities(&e);
    assert_eq!(ents.len(), 1);
    assert!(matches!(ents[0], EntityRef::Face(_)));

    // Precision (no entities)
    let e = ValidationError::Precision {
        diagnostic: crate::precision::PrecisionDiagnostic {
            severity: crate::precision::DiagnosticSeverity::Info,
            kind: crate::precision::PrecisionDiagnosticKind::ExceedsMaxBound,
            entities: Vec::new(),
            measured: 0.0,
            threshold: 1.0,
            remediation: crate::precision::Remediation::None,
        },
    };
    let ents = extract_entities(&e);
    assert!(ents.is_empty());
}

#[test]
fn extract_measurements_covers_all_variants() {
    // OpenLoop carries gap
    let e = ValidationError::OpenLoop {
        loop_: Default::default(),
        gap: 0.123,
    };
    let (m, t) = extract_measurements(&e);
    assert_eq!(m, Some(0.123));
    assert!(t.is_none());

    // CurveSurfaceDeviation carries max_deviation
    let e = ValidationError::CurveSurfaceDeviation {
        edge: Default::default(),
        max_deviation: 0.456,
    };
    let (m, _) = extract_measurements(&e);
    assert_eq!(m, Some(0.456));

    // Non-measurement variant returns None
    let e = ValidationError::BrokenTwin {
        half_edge: Default::default(),
    };
    let (m, t) = extract_measurements(&e);
    assert!(m.is_none());
    assert!(t.is_none());
}

#[test]
fn validation_error_to_defect_preserves_info() {
    let e = ValidationError::ShellGap {
        edge: Default::default(),
        gap: 0.05,
    };
    let defect = validation_error_to_defect(&e, Severity::Error);
    assert_eq!(defect.category, DefectCategory::Topological);
    assert_eq!(defect.severity, Severity::Error);
    assert!(defect.description.contains("gap") || defect.description.contains("Gap"));
    assert_eq!(defect.measured_value, Some(0.05));
    assert!(!defect.entities.is_empty());
}

#[test]
fn check_face_orientation_on_degenerate_surface() {
    use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let fh = brep.faces.keys().next().unwrap();
    // Zero-radius sphere has no valid normal -> orientation failure.
    brep.faces[fh].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 0.0,
        },
        orientation: true,
    });
    let mut report = HealingReport::new();
    check_face_orientation(&brep, &mut report);
    assert!(!report.is_empty());
    assert_eq!(report.defects()[0].severity, Severity::Warning);
}

#[test]
fn check_shell_closure_detects_gap_from_moved_vertex() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Move a vertex far away to create a watertight gap.
    let vh = brep.vertices.keys().next().unwrap();
    brep.vertices[vh].position = Point3::new(999.0, 999.0, 999.0);
    let mut report = HealingReport::new();
    check_shell_closure(&brep, &mut report);
    assert!(
        !report.is_empty(),
        "Moving a vertex should cause a shell closure gap"
    );
}

#[test]
fn check_all_on_corrupted_brep_collects_multiple_defects() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Corrupt multiple things: move a vertex and break a twin.
    let vh = brep.vertices.keys().next().unwrap();
    brep.vertices[vh].position = Point3::new(100.0, 100.0, 100.0);
    let he = brep.half_edges.keys().next().unwrap();
    brep.half_edges[he].twin = Default::default();
    let tol = Tolerance::default();
    let mut report = HealingReport::new();
    check_all(&brep, &tol, &mut report);
    assert!(
        report.len() >= 2,
        "Expected at least 2 defects from multiple corruptions, got {}",
        report.len()
    );
}

#[test]
fn check_inner_loop_orientation_on_face_with_no_inner_loops() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut report = HealingReport::new();
    check_inner_loop_orientation(&brep, &mut report);
    assert!(
        report.is_empty(),
        "Cube has no inner loops, should be clean"
    );
}

#[test]
fn check_manifold_edges_on_valid_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut report = HealingReport::new();
    check_manifold_edges(&brep, &mut report);
    assert!(report.is_empty(), "Valid cube should have manifold edges");
}

#[test]
fn check_vertex_consistency_on_valid_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut report = HealingReport::new();
    check_vertex_consistency(&brep, &mut report);
    assert!(
        report.is_empty(),
        "Valid cube should pass vertex consistency"
    );
}

#[test]
fn check_half_edge_cycles_on_valid_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let mut report = HealingReport::new();
    check_half_edge_cycles(&brep, &mut report);
    assert!(
        report.is_empty(),
        "Valid cube should have valid half-edge cycles"
    );
}

#[test]
fn validation_error_to_defect_for_broken_cycle() {
    let e = ValidationError::BrokenCycle {
        half_edge: Default::default(),
        loop_: Default::default(),
    };
    let defect = validation_error_to_defect(&e, Severity::Error);
    assert_eq!(defect.category, DefectCategory::Topological);
    assert_eq!(defect.severity, Severity::Error);
    assert_eq!(defect.entities.len(), 2);
}

#[test]
fn validation_error_to_defect_for_inconsistent_vertex() {
    let e = ValidationError::InconsistentVertex {
        vertex: Default::default(),
        half_edge: Default::default(),
    };
    let defect = validation_error_to_defect(&e, Severity::Warning);
    assert_eq!(defect.severity, Severity::Warning);
    assert_eq!(defect.entities.len(), 2);
    assert!(matches!(defect.entities[0], EntityRef::Vertex(_)));
}

#[test]
fn validation_error_to_defect_for_curve_surface_deviation() {
    let e = ValidationError::CurveSurfaceDeviation {
        edge: Default::default(),
        max_deviation: 0.789,
    };
    let defect = validation_error_to_defect(&e, Severity::Error);
    assert_eq!(defect.measured_value, Some(0.789));
    assert!(defect.entities.len() == 1);
    assert!(matches!(defect.entities[0], EntityRef::Edge(_)));
}
