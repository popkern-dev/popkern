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
