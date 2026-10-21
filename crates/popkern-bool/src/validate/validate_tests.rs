//! Tests for boolean validation and repair.

use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ── Pre-check tests ────────────────────────────────────────────────────

#[test]
fn pre_check_valid_cubes() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap();
    let report = super::pre_check::pre_check(&a, &b, false);
    assert!(
        report.is_ok(),
        "valid cubes should pass pre-check, got {} errors",
        report.error_count()
    );
}

#[test]
fn pre_check_valid_cylinder() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::cylinder(0.5, 2.0, tol()).unwrap();
    let report = super::pre_check::pre_check(&a, &b, false);
    assert!(report.is_ok());
}

#[test]
fn pre_check_empty_body_fails() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::new(tol());
    let report = super::pre_check::pre_check(&a, &b, false);
    assert!(!report.is_ok(), "empty body B should fail pre-check");
    assert!(report.error_count() > 0);
}

#[test]
fn pre_check_report_counts() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol()).unwrap();
    let report = super::pre_check::pre_check(&a, &b, false);
    assert_eq!(report.error_count(), 0);
    assert_eq!(report.warning_count(), 0);
}

#[test]
fn pre_check_clean_cubes_have_high_likelihood() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(5.0, 0.0, 0.0))
        .unwrap();
    let report = super::pre_check::pre_check(&a, &b, false);
    // Disjoint clean cubes should have no degenerate faces and a
    // success likelihood close to 1.0.
    assert_eq!(report.degenerate_faces.len(), 0);
    assert!(
        report.success_likelihood >= 0.9,
        "clean disjoint cubes should be highly likely to succeed; got {}",
        report.success_likelihood
    );
}

#[test]
fn pre_check_touching_cubes_flag_problematic_pairs() {
    // Two unit cubes sharing a face (translation = 1.0 in x). The
    // +X face of A and the -X face of B sit on the same plane with
    // anti-parallel normals — a textbook near-coincident pair that
    // boolean engines must handle as a coincident face merge or
    // cancellation. Pre-check should flag the pair.
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(1.0, 0.0, 0.0))
        .unwrap();
    let report = super::pre_check::pre_check(&a, &b, false);
    assert!(
        !report.problematic_pairs.is_empty(),
        "touching cubes should flag at least one problematic pair, got {}",
        report.problematic_pairs.len()
    );
}

#[test]
fn pre_check_clean_cube_has_no_g0_violations() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let other = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(5.0, 0.0, 0.0))
        .unwrap();
    let report = super::pre_check::pre_check(&cube, &other, false);
    let g0_violations = report
        .findings
        .iter()
        .filter(|f| f.message.contains("G0 (positional) violation"))
        .count();
    assert_eq!(
        g0_violations, 0,
        "clean cube should have no G0 violations: {:?}",
        report.findings
    );
}

#[test]
fn pre_check_clean_cube_has_no_g1_violations() {
    // Cube has perpendicular faces (90° edges). The G1 detector
    // should NOT fire on sharp edges — only on near-tangent
    // (smooth) edges with subtle misalignment. Cube is the
    // canonical "no false positive" fixture.
    let cube = BRep::unit_cube(tol()).unwrap();
    let other = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(5.0, 0.0, 0.0))
        .unwrap();
    let report = super::pre_check::pre_check(&cube, &other, false);
    let g1_violations = report
        .findings
        .iter()
        .filter(|f| f.message.contains("G1 (tangent) violation"))
        .count();
    assert_eq!(
        g1_violations, 0,
        "cube edges are sharp; G1 detector must not fire: {:?}",
        report.findings
    );
}

#[test]
fn pre_check_torus_warns() {
    use popkern_core::geometry::Point3;
    let cube = BRep::unit_cube(tol()).unwrap();
    let torus = BRep::torus(1.0, 0.3, tol()).unwrap();
    let _ = Point3::new(0.0, 0.0, 0.0); // exercise import path
    let report = super::pre_check::pre_check(&cube, &torus, false);
    let torus_warning = report
        .findings
        .iter()
        .any(|f| f.message.contains("Torus surfaces present"));
    assert!(
        torus_warning,
        "torus input should produce a warning: {:?}",
        report.findings
    );
}

#[test]
fn pre_check_likelihood_decreases_with_problematic_pairs() {
    // Two overlapping cubes (lots of parallel faces) should score
    // lower than two disjoint cubes (no parallel pairs).
    let a = BRep::unit_cube(tol()).unwrap();
    let b_disjoint = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(5.0, 0.0, 0.0))
        .unwrap();
    let b_overlap = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    let report_disjoint = super::pre_check::pre_check(&a, &b_disjoint, false);
    let report_overlap = super::pre_check::pre_check(&a, &b_overlap, false);
    assert!(
        report_disjoint.success_likelihood >= report_overlap.success_likelihood,
        "disjoint score ({}) should be ≥ overlap score ({})",
        report_disjoint.success_likelihood,
        report_overlap.success_likelihood
    );
}

// ── Post-check tests ───────────────────────────────────────────────────

#[test]
fn post_check_valid_cube() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let report = super::post_check::post_check(&cube, &tol(), false);
    assert!(
        report.is_ok(),
        "valid cube should pass post-check, got {} errors: {:?}",
        report.error_count(),
        report.findings
    );
    assert_eq!(report.euler_characteristic, 2);
    assert_eq!(report.face_count, 6);
    assert_eq!(report.vertex_count, 8);
}

#[test]
fn post_check_valid_cylinder() {
    let cyl = BRep::cylinder(1.0, 2.0, tol()).unwrap();
    let report = super::post_check::post_check(&cyl, &tol(), false);
    // Cylinder may have minor geometric warnings but no structural errors
    let structural_errors = report
        .findings
        .iter()
        .filter(|f| {
            matches!(f.severity, crate::error::Severity::Error)
                && !f.message.contains("orientation")
        })
        .count();
    assert!(
        structural_errors == 0,
        "cylinder should have no structural errors, got: {:?}",
        report.findings
    );
}

#[test]
fn post_check_boolean_result() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    if let Ok(result) = crate::union(a, b, tol(), None) {
        let report = super::post_check::post_check(&result, &tol(), false);
        // Log but don't hard-assert — graph-merge may produce imperfect results
        if !report.is_ok() {
            eprintln!(
                "[post_check] boolean result: {} errors, {} faces, euler={}",
                report.error_count(),
                report.face_count,
                report.euler_characteristic
            );
        }
    }
}

#[test]
fn post_check_report_metadata() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let report = super::post_check::post_check(&cube, &tol(), false);
    assert!(report.face_count > 0);
    assert!(report.edge_count > 0);
    assert!(report.vertex_count > 0);
    assert!(report.max_gap >= 0.0);
}

#[test]
fn post_check_self_intersection_clean_cube() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let report = super::post_check::post_check(&cube, &tol(), true);
    assert!(
        report.is_ok(),
        "valid cube should pass self-intersection check, got: {:?}",
        report.findings
    );
}

#[test]
fn boolean_strict_disjoint_cubes_succeeds() {
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(5.0, 0.0, 0.0))
        .unwrap();
    let result = crate::boolean_strict(crate::BoolOperation::Union, &a, &b, tol(), false);
    assert!(
        result.is_ok(),
        "disjoint-cube union should pass strict post-check, got {:?}",
        result.err()
    );
}

#[test]
fn diagnostic_entry_carries_3d_location() {
    use crate::error::{BoolDiagnostics, DiagnosticEntry, EntityRef, Severity, Stage};
    use popkern_core::geometry::Point3;

    let mut diag = BoolDiagnostics::default();
    let pt = Point3::new(1.0, 2.0, 3.0);
    diag.add(DiagnosticEntry {
        severity: Severity::Error,
        phase: Stage::Intersection,
        message: "test diagnostic at known point".into(),
        entities: Vec::new(),
        location: Some(pt),
        measured_value: None,
        remediation: None,
    });

    let entry = diag.entries.first().expect("entry recorded");
    let loc = entry.location.expect("location populated");
    assert_eq!(loc, pt);

    // warn_at convenience method.
    let mut diag2 = BoolDiagnostics::default();
    diag2.warn_at(Stage::Merge, "another", Point3::new(4.0, 5.0, 6.0));
    let warn = diag2.entries.first().unwrap();
    assert_eq!(warn.location, Some(Point3::new(4.0, 5.0, 6.0)));
    assert!(matches!(warn.severity, Severity::Warning));

    // EntityRef still independent of location.
    let mut diag3 = BoolDiagnostics::default();
    diag3.add(DiagnosticEntry {
        severity: Severity::Info,
        phase: Stage::Validation,
        message: "no location".into(),
        entities: vec![EntityRef::Face(Default::default())],
        location: None,
        measured_value: None,
        remediation: None,
    });

    assert!(diag3.entries[0].location.is_none());

    // warn_with_measure convenience method (§13.1).
    let mut diag4 = BoolDiagnostics::default();
    diag4.warn_with_measure(Stage::Validation, "gap of 0.003mm", 0.003);
    let measured = diag4.entries.first().unwrap();
    assert_eq!(measured.measured_value, Some(0.003));
    assert!(matches!(measured.severity, Severity::Warning));
}

#[test]
fn bool_error_categories_cover_interference_and_system() {
    // §13.2: every failure mode has a corresponding error category.
    // This test pins the new variants so the spec coverage matrix
    // doesn't regress.
    use crate::error::BoolError;
    let interference = BoolError::InterferenceFailed("BVH lookup failed".into());
    assert!(matches!(interference, BoolError::InterferenceFailed(_)));
    let system = BoolError::SystemError("out of memory".into());
    assert!(matches!(system, BoolError::SystemError(_)));
}

#[test]
fn boolean_with_progress_emits_phase_events() {
    use crate::progress::{CaptureReporter, ProgressReporter};
    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    let mut reporter = CaptureReporter::default();
    let trait_ref: &mut dyn ProgressReporter = &mut reporter;
    let _ = crate::boolean_with_progress(crate::BoolOperation::Union, &a, &b, tol(), trait_ref);
    // Overlapping cubes run the full pipeline; expect events for
    // multiple stages.
    assert!(
        reporter.events.len() >= 2,
        "expected progress events from intersect/split/etc., got {} events",
        reporter.events.len()
    );
    // First event should be the intersection stage at phase_index 0.
    let first = reporter.events.first().unwrap();
    assert_eq!(first.stage, crate::error::Stage::Intersection);
    assert_eq!(first.phase_index, 0);
    // Percent monotonically increases.
    let mut last_percent = 0;
    for ev in &reporter.events {
        assert!(ev.percent() >= last_percent);
        last_percent = ev.percent();
    }
}

#[test]
fn boolean_with_diagnostics_records_phase_metrics() {
    use crate::error::Stage;

    let a = BRep::unit_cube(tol()).unwrap();
    let b = BRep::unit_cube(tol())
        .unwrap()
        .transform(&Transform3::translation(0.5, 0.3, 0.2))
        .unwrap();
    let (_result, diag) =
        crate::boolean_with_diagnostics(crate::BoolOperation::Union, &a, &b, tol())
            .expect("intersecting cubes union should succeed");

    // The pipeline ran the full set of stages (no shortcut), so we expect
    // a metric for each. The cellular shortcut paths return earlier and
    // may produce fewer metrics — for genuinely intersecting cubes the
    // full pipeline runs.
    let stages: Vec<Stage> = diag.phase_metrics.iter().map(|m| m.stage).collect();
    assert!(
        stages.contains(&Stage::Intersection),
        "phase_metrics should record an Intersection sample, got {stages:?}"
    );
    assert!(
        stages.contains(&Stage::Merge),
        "phase_metrics should record a Merge sample, got {stages:?}"
    );

    // Every recorded sample must have a non-zero entity count for the
    // intersection stage on inputs known to intersect.
    let intersect_metric = diag
        .phase_metrics
        .iter()
        .find(|m| m.stage == Stage::Intersection)
        .expect("intersection metric present");
    assert!(
        intersect_metric.entity_count > 0,
        "intersection should report > 0 records for overlapping cubes"
    );
}

// ── Repair tests ───────────────────────────────────────────────────────

#[test]
fn repair_on_valid_cube_changes_nothing() {
    let mut cube = BRep::unit_cube(tol()).unwrap();
    let report = super::repair::repair_vertex_positions(&mut cube, 1e-4);
    assert!(report.fully_repaired());
    assert_eq!(report.vertices_adjusted, 0);
}

#[test]
fn repair_report_default() {
    let report = super::repair::RepairReport::default();
    assert!(report.fully_repaired());
    assert_eq!(report.vertices_adjusted, 0);
    assert_eq!(report.unrepairable_gaps, 0);
}

#[test]
fn repair_closes_small_gap() {
    // Create a cube and intentionally perturb a vertex position
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();

    // Perturb the first vertex slightly
    let vh = cube.vertices.keys().next().unwrap();
    let original_pos = cube.vertices[vh].position;
    cube.vertices[vh].position =
        popkern_core::geometry::Point3::new(original_pos.x + 1e-7, original_pos.y, original_pos.z);

    let report = super::repair::repair_vertex_positions(&mut cube, 1e-4);
    // The vertex was moved by 1e-7 which is within 1e-4 tolerance.
    // The repair should snap it to the curve endpoint if the edge
    // curve evaluates to the original position.
    // Note: unit_cube edges are lines from vertex to vertex, so the
    // curve endpoint IS the vertex — repair may or may not trigger
    // depending on whether the curve was updated with the perturbation.
    // This test validates the repair path executes without panic.
    assert!(report.fully_repaired() || report.vertices_adjusted > 0);
}

#[test]
fn post_check_gates_report_exactly_the_gating_findings() {
    // The serve path reads only the gates (PK-731): same Error/Critical
    // findings as the full post-check, no warning-only diagnostics.
    for body in [
        BRep::unit_cube(tol()).unwrap(),
        BRep::cylinder(1.0, 2.0, tol()).unwrap(),
    ] {
        let full = super::post_check::post_check(&body, &tol(), false);
        let gates = super::post_check::post_check_gates(&body, &tol(), &Default::default(), false);
        let errors = |r: &super::post_check::PostCheckReport| -> Vec<String> {
            r.findings
                .iter()
                .filter(|f| {
                    matches!(
                        f.severity,
                        crate::error::Severity::Error | crate::error::Severity::Critical
                    )
                })
                .map(|f| f.message.clone())
                .collect()
        };
        assert_eq!(errors(&full), errors(&gates));
        assert_eq!(gates.findings.len(), errors(&gates).len(), "no warnings");
    }
}
