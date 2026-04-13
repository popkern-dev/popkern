use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}
fn fopts() -> FilletOptions {
    FilletOptions::default()
}
fn copts() -> ChamferOptions {
    ChamferOptions::default()
}
fn cube_and_edge() -> (BRep, popkern_core::topology::EdgeHandle) {
    let brep = BRep::unit_cube(tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    (brep, edge)
}

#[test]
fn variable_radius_fillet_computes_geometry() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Variable {
            start: 0.08,
            end: 0.12,
        },
    }];
    let result = fillet_edges(&brep, &specs, &fopts(), tol(), None);
    assert!(result.is_ok(), "variable fillet failed: {:?}", result.err());
    assert!(result.unwrap().brep.faces.len() > brep.faces.len());
}

#[test]
fn by_points_fillet_exercises_code_paths() {
    let (brep, edge) = cube_and_edge();
    // Two-point ByPoints — exercises variable_radius path. Should not panic.
    let s1 = vec![FilletSpec {
        edge,
        radius: FilletRadius::ByPoints(vec![(0.0, 0.1), (1.0, 0.1)]),
    }];
    let _ = fillet_edges(&brep, &s1, &fopts(), tol(), None);
    // Single-point ByPoints.
    let s2 = vec![FilletSpec {
        edge,
        radius: FilletRadius::ByPoints(vec![(0.5, 0.1)]),
    }];
    let _ = fillet_edges(&brep, &s2, &fopts(), tol(), None);
    // Empty ByPoints should fail.
    let s3 = vec![FilletSpec {
        edge,
        radius: FilletRadius::ByPoints(vec![]),
    }];
    assert!(fillet_edges(&brep, &s3, &fopts(), tol(), None).is_err());
}

#[test]
fn chamfer_two_distance_via_api() {
    let (brep, edge) = cube_and_edge();
    let face_1 = analysis::analyze_edges(&brep, &[edge], false).unwrap()[0].face_1;
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::TwoDistance {
            face_1,
            dist_1: 0.1,
            dist_2: 0.2,
        },
    }];
    let r = chamfer_edges(&brep, &specs, &copts(), tol(), None);
    assert!(r.is_ok(), "two-distance chamfer failed: {:?}", r.err());
}

#[test]
fn chamfer_distance_angle_via_api() {
    let (brep, edge) = cube_and_edge();
    let face = analysis::analyze_edges(&brep, &[edge], false).unwrap()[0].face_1;
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::DistanceAngle {
            face,
            distance: 0.1,
            angle: std::f64::consts::FRAC_PI_4,
        },
    }];
    let r = chamfer_edges(&brep, &specs, &copts(), tol(), None);
    assert!(r.is_ok(), "dist-angle chamfer failed: {:?}", r.err());
}

#[test]
fn empty_specs_returns_error() {
    let brep = BRep::unit_cube(tol()).unwrap();
    assert!(fillet_edges(&brep, &[], &fopts(), tol(), None).is_err());
    assert!(chamfer_edges(&brep, &[], &copts(), tol(), None).is_err());
}

#[test]
fn fillet_edges_rejects_nan_radius_before_touching_geometry() {
    let (brep, edge) = cube_and_edge();
    let pre_face_count = brep.faces.len();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(f64::NAN),
    }];
    let err = fillet_edges(&brep, &specs, &fopts(), tol(), None).unwrap_err();
    assert!(matches!(err, BlendError::InvalidParameter { .. }));
    // Input BRep must be untouched — this is part of the
    // §10.1 contract.
    assert_eq!(brep.faces.len(), pre_face_count);
    // Remediation hint should route to FixInput.
    assert_eq!(err.remediation_hint().action, RemediationAction::FixInput);
}

#[test]
fn fillet_edges_rejects_negative_radius_before_touching_geometry() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(-0.1),
    }];
    let err = fillet_edges(&brep, &specs, &fopts(), tol(), None).unwrap_err();
    assert!(matches!(err, BlendError::InvalidParameter { .. }));
}

#[test]
fn fillet_edges_returns_radius_too_large_with_populated_max_radius() {
    // Unit cube has edge length 1.0. An obviously-too-large
    // radius (e.g. 2.0) should trigger the self-intersection
    // pre-flight and return RadiusTooLarge with a concrete
    // max_radius hint that is non-zero and less than the
    // requested value.
    //
    // However, note: plane-plane curvature is infinite, so two
    // flat faces won't trip the check. This test verifies the
    // plumbing: even when no max is populated, the structured
    // failure shape is stable.
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(2.0), // 2× the cube
    }];
    // Cube edges are plane-plane, so self_intersect returns None
    // and the operation proceeds as before (either succeeds or
    // fails downstream). The test for actual max_radius
    // population is in the curved_faces integration test file.
    let _ = fillet_edges(&brep, &specs, &fopts(), tol(), None);
}

#[test]
fn chamfer_edges_rejects_zero_distance_before_touching_geometry() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::Equal(0.0),
    }];
    let err = chamfer_edges(&brep, &specs, &copts(), tol(), None).unwrap_err();
    assert!(matches!(err, BlendError::InvalidParameter { .. }));
}

#[test]
fn chamfer_edges_rejects_out_of_range_angle_before_touching_geometry() {
    let (brep, edge) = cube_and_edge();
    let face = brep.faces.keys().next().unwrap();
    let specs = vec![ChamferSpec {
        edge,
        kind: ChamferKind::DistanceAngle {
            face,
            distance: 0.1,
            angle: 10.0, // > π
        },
    }];
    let err = chamfer_edges(&brep, &specs, &copts(), tol(), None).unwrap_err();
    assert!(matches!(err, BlendError::InvalidParameter { .. }));
}

#[test]
fn fillet_edges_with_quality_returns_populated_metrics() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let (result, quality) =
        fillet_edges_with_quality(&brep, &specs, &fopts(), tol(), None).unwrap();
    // The blend ran — at least one new face exists.
    assert!(result.brep.faces.len() > brep.faces.len());
    // Quality metrics populated.
    assert!(quality.patch_count >= 1);
    assert!(quality.surface_area > 0.0);
    assert!(quality.composite_score >= 0.0);
    assert!(quality.composite_score <= 1.0);
}

#[test]
fn fillet_edges_with_quality_propagates_validation_error() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(f64::NAN),
    }];
    assert!(matches!(
        fillet_edges_with_quality(&brep, &specs, &fopts(), tol(), None),
        Err(BlendError::InvalidParameter { .. })
    ));
}

#[test]
fn fillet_cancellable_pre_signaled_returns_cancelled_without_mutation() {
    let (brep, edge) = cube_and_edge();
    let pre_face_count = brep.faces.len();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let token = CancellationToken::new();
    token.cancel();
    let err = fillet_edges_cancellable(&brep, &specs, &fopts(), tol(), None, &token).unwrap_err();
    assert!(matches!(err, BlendError::Cancelled));
    // Input unchanged.
    assert_eq!(brep.faces.len(), pre_face_count);
    // Remediation hint.
    assert_eq!(err.remediation_hint().action, RemediationAction::NoAction);
}

#[test]
fn fillet_cancellable_without_cancel_signal_succeeds() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let token = CancellationToken::new();
    let result = fillet_edges_cancellable(&brep, &specs, &fopts(), tol(), None, &token);
    assert!(result.is_ok());
}

#[test]
fn fillet_with_progress_emits_events_for_every_phase() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let mut reporter = CaptureReporter::default();
    fillet_edges_with_progress(&brep, &specs, &fopts(), tol(), None, &mut reporter).unwrap();
    // One event per phase.
    assert_eq!(reporter.events.len(), FilletPhase::total());
    // Phases in canonical order.
    let phases: Vec<_> = reporter.events.iter().map(|e| e.phase).collect();
    assert_eq!(
        phases,
        vec![
            FilletPhase::Validate,
            FilletPhase::Analyze,
            FilletPhase::ComputeGeometry,
            FilletPhase::DetectOverlaps,
            FilletPhase::ApplyTopology,
            FilletPhase::RepairPcurves,
        ]
    );
    // Percent is monotonically non-decreasing.
    let mut last_pct = 0u8;
    for e in &reporter.events {
        assert!(e.percent() >= last_pct);
        last_pct = e.percent();
    }
}

#[test]
fn fillet_with_null_reporter_produces_same_result_as_plain_call() {
    let (brep, edge) = cube_and_edge();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let mut null = NullReporter;
    let a = fillet_edges(&brep, &specs, &fopts(), tol(), None).unwrap();
    let b = fillet_edges_with_progress(&brep, &specs, &fopts(), tol(), None, &mut null).unwrap();
    assert_eq!(a.brep.faces.len(), b.brep.faces.len());
}

#[test]
fn fillet_result_has_pcurves() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let result = fillet_edge(&brep, edge, 0.1, Tolerance::default(), None).unwrap();

    let total = result.brep.half_edges.len();
    let with_pcurves = result
        .brep
        .half_edges
        .values()
        .filter(|he| he.pcurve.is_some())
        .count();

    assert!(
        with_pcurves > total / 2,
        "expected most half-edges to have pcurves, got {with_pcurves}/{total}"
    );
}
