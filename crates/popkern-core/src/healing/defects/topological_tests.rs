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
            axis: crate::geometry::Vec3::Z,
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

// -----------------------------------------------------------------------
// Tests targeting uncovered branches
// -----------------------------------------------------------------------

use crate::constructors::builder::BRepBuilder;
use crate::geometry::Vec3;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

#[test]
fn inner_loop_wrong_winding_detected() {
    let mut b = BRepBuilder::new(Tolerance::default());
    // Outer loop: CCW square in XY plane.
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(4.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(4.0, 4.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 4.0, 0.0));
    // Inner loop: same CCW winding (wrong -- should be CW).
    let h0 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let h1 = b.add_vertex(Point3::new(3.0, 1.0, 0.0));
    let h2 = b.add_vertex(Point3::new(3.0, 3.0, 0.0));
    let h3 = b.add_vertex(Point3::new(1.0, 3.0, 0.0));
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let (_, ho01, _) = b.add_edge(v0, v1);
    let (_, ho12, _) = b.add_edge(v1, v2);
    let (_, ho23, _) = b.add_edge(v2, v3);
    let (_, ho30, _) = b.add_edge(v3, v0);
    let (_, hi01, _) = b.add_edge(h0, h1);
    let (_, hi12, _) = b.add_edge(h1, h2);
    let (_, hi23, _) = b.add_edge(h2, h3);
    let (_, hi30, _) = b.add_edge(h3, h0);
    let face = b.add_face(shell);
    b.add_loop(&[ho01, ho12, ho23, ho30], face).unwrap();
    b.add_inner_loop(&[hi01, hi12, hi23, hi30], face).unwrap();
    let brep = b.build_unchecked();
    let mut report = HealingReport::new();
    check_inner_loop_orientation(&brep, &mut report);
    assert!(
        !report.is_empty(),
        "same-winding inner loop should be flagged"
    );
    assert_eq!(report.defects()[0].severity, Severity::Warning);
    assert!(report.defects()[0].description.contains("same winding"));
}

#[test]
fn boundary_simplicity_detects_crossing_edges() {
    let mut b = BRepBuilder::new(Tolerance::default());
    // Bowtie shape: edges v0->v1 and v2->v3 cross at (1,1).
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(2.0, 2.0, 0.0));
    let v2 = b.add_vertex(Point3::new(2.0, 0.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 2.0, 0.0));
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let (e01, h01, _) = b.add_edge(v0, v1);
    let (e12, h12, _) = b.add_edge(v1, v2);
    let (e23, h23, _) = b.add_edge(v2, v3);
    let (e30, h30, _) = b.add_edge(v3, v0);
    for (eh, s, e) in [
        (e01, [0., 0., 0.], [2., 2., 0.]),
        (e12, [2., 2., 0.], [2., 0., 0.]),
        (e23, [2., 0., 0.], [0., 2., 0.]),
        (e30, [0., 2., 0.], [0., 0., 0.]),
    ] {
        let dir = Vec3::new(e[0] - s[0], e[1] - s[1], e[2] - s[2]);
        b.set_edge_curve(
            eh,
            CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: Point3::new(s[0], s[1], s[2]),
                    direction: dir,
                },
                t_start: 0.0,
                t_end: 1.0,
            },
        );
    }
    let face = b.add_face(shell);
    b.add_loop(&[h01, h12, h23, h30], face).unwrap();
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
    let brep = b.build_unchecked();
    // Use a very large model scale to get a linear tolerance > 1.0
    // so that crossing sample points (which are ~0.0 apart at the crossing)
    // are detected within tolerance.
    let loose_tol = Tolerance::from_model_scale(2_000_000.0);
    let mut report = HealingReport::new();
    check_boundary_simplicity(&brep, &loose_tol, &mut report);
    assert!(
        !report.is_empty(),
        "crossing edges should trigger self-intersection"
    );
    assert_eq!(report.defects()[0].severity, Severity::Error);
    assert!(
        report.defects()[0]
            .description
            .contains("self-intersection")
    );
}
