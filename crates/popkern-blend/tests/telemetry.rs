//! §18 telemetry end-to-end tests.
//!
//! Verifies that `fillet_edges_with_progress` + `CaptureReporter`
//! produces a queryable [`popkern_blend::FilletTelemetry`] report
//! with per-edge timing, Newton iteration counts, and aggregate
//! accessors (mean / max / histogram).

use popkern_blend::{
    CaptureReporter, FilletOptions, FilletPhase, FilletRadius, FilletSpec,
    fillet_edges_with_progress,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------------------
// Per-edge event emission
// ---------------------------------------------------------------------------

#[test]
fn capture_reporter_records_one_edge_event_per_spec() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .take(3)
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let original_count = specs.len();

    let mut reporter = CaptureReporter::default();
    let _result = fillet_edges_with_progress(
        &cube,
        &specs,
        &FilletOptions::default(),
        tol(),
        None,
        &mut reporter,
    );
    let tel = reporter.telemetry();

    assert_eq!(
        tel.edge_events.len(),
        original_count,
        "expected one edge event per spec, got {} events for {} specs",
        tel.edge_events.len(),
        original_count
    );
    // Every event should have a non-zero duration.
    for ev in &tel.edge_events {
        assert!(
            ev.duration.as_nanos() > 0,
            "edge event should have measurable duration, got {:?}",
            ev.duration
        );
    }
}

// ---------------------------------------------------------------------------
// Per-corner event emission (PK-276, acceptance row 413)
// ---------------------------------------------------------------------------

/// The 12-edge cube fillet constructs 8 vertex-blend corners inside
/// ApplyTopology; the reporter must tick once per corner with
/// monotone 0-based indices and an honest total. Before PK-276 the
/// dominant apply phase emitted exactly one event, at completion.
#[test]
fn capture_reporter_records_one_corner_event_per_vertex_blend() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    assert_eq!(specs.len(), 12);

    let mut reporter = CaptureReporter::default();
    let result = fillet_edges_with_progress(
        &cube,
        &specs,
        &FilletOptions::default(),
        tol(),
        None,
        &mut reporter,
    );
    assert!(result.is_ok(), "12-edge cube fillet serves: {result:?}");

    assert_eq!(
        reporter.corner_events.len(),
        8,
        "cube has 8 corners; expected one corner event each, got {}",
        reporter.corner_events.len()
    );
    let mut seen_vertices = std::collections::BTreeSet::new();
    for (i, ev) in reporter.corner_events.iter().enumerate() {
        assert_eq!(
            ev.corner_index, i,
            "corner indices must be 0-based monotone"
        );
        assert_eq!(ev.total_corners, 8, "total must be the honest corner count");
        assert!(
            ev.duration.as_nanos() > 0,
            "corner event should have measurable duration, got {:?}",
            ev.duration
        );
        assert!(
            seen_vertices.insert(ev.vertex),
            "corner vertices must be distinct, {:?} repeated",
            ev.vertex
        );
    }
}

// ---------------------------------------------------------------------------
// Phase + edge coverage
// ---------------------------------------------------------------------------

#[test]
fn capture_reporter_records_all_pipeline_phases() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];

    let mut reporter = CaptureReporter::default();
    let _ = fillet_edges_with_progress(
        &cube,
        &specs,
        &FilletOptions::default(),
        tol(),
        None,
        &mut reporter,
    );
    let tel = reporter.telemetry();

    // Phase events should cover the standard pipeline sequence.
    // We accept phase events in any subset (e.g. if the pipeline
    // short-circuits on a failure) but every phase that DID emit
    // must carry a nonzero duration.
    for e in &tel.phase_events {
        assert!(
            e.duration.as_nanos() > 0,
            "phase {:?} event should have nonzero duration, got {:?}",
            e.phase,
            e.duration
        );
    }
    // The Validate phase always runs regardless of downstream
    // outcome, so we expect at least one phase event in total.
    assert!(!tel.phase_events.is_empty(), "expected ≥1 phase event");
    assert!(
        tel.phase_events
            .iter()
            .any(|e| e.phase == FilletPhase::Validate),
        "Validate phase event should always be emitted"
    );
}

// ---------------------------------------------------------------------------
// Newton iteration telemetry on the walker path (cone base)
// ---------------------------------------------------------------------------

#[test]
fn walker_path_populates_newton_iterations() {
    // The cone base fillet routes through the phase-1b walker,
    // which populates `center_result.contact_solutions` and thus
    // `EdgeGeometry::iterations_used`. The telemetry should
    // report `Some(n)` for that edge, not `None`.
    let cone = BRep::cone(1.0, 2.0, tol()).expect("cone BRep should build");
    let chosen = cone
        .edges
        .keys()
        .find(|&e| {
            let he = cone.edges[e].half_edge;
            let twin = cone.half_edges[he].twin;
            let face_a = cone.loops[cone.half_edges[he].loop_].face;
            let face_b = cone.loops[cone.half_edges[twin].loop_].face;
            face_a != face_b
        })
        .unwrap_or_else(|| cone.edges.keys().next().unwrap());
    let specs = vec![FilletSpec {
        edge: chosen,
        radius: FilletRadius::Constant(0.005),
    }];

    let mut reporter = CaptureReporter::default();
    let _ = fillet_edges_with_progress(
        &cone,
        &specs,
        &FilletOptions::default(),
        tol(),
        None,
        &mut reporter,
    );
    let tel = reporter.telemetry();

    // Acceptable outcomes:
    //  - edge succeeded with Some(n) Newton count
    //  - edge failed (e.g. cone-corner topology phase) with None
    // We only assert the walker populates iterations when it
    // actually ran, not that it always runs.
    if let Some(ev) = tel.edge_events.first()
        && ev.succeeded
    {
        assert!(
            ev.newton_iterations.is_some(),
            "successful walker-path edge should report Newton iterations, got None"
        );
        let count = ev.newton_iterations.unwrap();
        assert!(
            count > 0,
            "walker ran but reported 0 total Newton iterations"
        );
        // And the aggregates should reflect this.
        assert_eq!(tel.total_newton_iterations(), count);
        assert_eq!(tel.max_newton_iterations(), Some(count));
        assert!(tel.mean_newton_iterations().is_some());
        let histogram = tel.newton_histogram();
        assert!(
            !histogram.is_empty(),
            "histogram should include at least the walker edge's bucket"
        );

        // Step rejections: the walker populates this whenever it
        // populates iterations — expected `Some(n)` with n ≥ 0.
        // On a healthy cone-base solve the line search usually
        // accepts on the first try, so 0 is a valid count.
        assert!(
            ev.step_rejections.is_some(),
            "successful walker-path edge should report step rejections, got None"
        );
        // Ratio accessor should also return Some here.
        assert!(tel.step_rejection_ratio().is_some());
    }
}

// ---------------------------------------------------------------------------
// Offset-intersect path leaves Newton iterations = None
// ---------------------------------------------------------------------------

#[test]
fn offset_intersect_path_leaves_newton_iterations_none() {
    // A plain cube fillet at a small radius takes the offset-
    // intersect vertex-blend path (no walker, no contact solutions),
    // so `iterations_used` stays None.
    let cube = BRep::unit_cube(tol()).unwrap();
    let edge = cube.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];

    let mut reporter = CaptureReporter::default();
    let _ = fillet_edges_with_progress(
        &cube,
        &specs,
        &FilletOptions::default(),
        tol(),
        None,
        &mut reporter,
    );
    let tel = reporter.telemetry();

    // Only check successful edges — failures naturally have None.
    for ev in &tel.edge_events {
        if ev.succeeded {
            assert!(
                ev.newton_iterations.is_none(),
                "offset-intersect path should leave iterations None, got {:?}",
                ev.newton_iterations
            );
            assert!(
                ev.step_rejections.is_none(),
                "offset-intersect path should leave step rejections None, got {:?}",
                ev.step_rejections
            );
        }
    }
    // Aggregates should report nothing when no walker edges ran.
    assert_eq!(tel.total_newton_iterations(), 0);
    assert_eq!(tel.mean_newton_iterations(), None);
    assert_eq!(tel.max_newton_iterations(), None);
    assert!(tel.newton_histogram().is_empty());
    // Step-rejection aggregates mirror Newton aggregates.
    assert_eq!(tel.total_step_rejections(), 0);
    assert_eq!(tel.mean_step_rejections(), None);
    assert_eq!(tel.max_step_rejections(), None);
    assert_eq!(tel.step_rejection_ratio(), None);
}

// ---------------------------------------------------------------------------
// Aggregate accessors
// ---------------------------------------------------------------------------

#[test]
fn telemetry_aggregate_accessors_are_consistent() {
    use popkern_blend::{EdgeTelemetry, FilletTelemetry};
    use std::time::Duration;

    let mut tel = FilletTelemetry::default();
    // Synthesize three edge events with known iteration counts.
    let dummy = popkern_core::topology::EdgeHandle::default();
    tel.edge_events = vec![
        EdgeTelemetry {
            edge: dummy,
            duration: Duration::from_micros(100),
            newton_iterations: Some(3),
            step_rejections: Some(0),
            succeeded: true,
        },
        EdgeTelemetry {
            edge: dummy,
            duration: Duration::from_micros(200),
            newton_iterations: Some(5),
            step_rejections: Some(2),
            succeeded: true,
        },
        EdgeTelemetry {
            edge: dummy,
            duration: Duration::from_micros(150),
            newton_iterations: Some(3),
            step_rejections: Some(1),
            succeeded: false,
        },
    ];

    assert_eq!(tel.edges_succeeded(), 2);
    assert_eq!(tel.edges_failed(), 1);
    assert_eq!(tel.total_newton_iterations(), 11);
    assert!((tel.mean_newton_iterations().unwrap() - (11.0 / 3.0)).abs() < 1e-12);
    assert_eq!(tel.max_newton_iterations(), Some(5));
    assert_eq!(
        tel.total_edge_duration(),
        Duration::from_micros(100 + 200 + 150)
    );

    // Histogram: two edges at 3 iter, one edge at 5 iter.
    let histogram = tel.newton_histogram();
    assert_eq!(histogram, vec![(3, 2), (5, 1)]);

    // Step rejection aggregates: 0 + 2 + 1 = 3 total across 3
    // edges, so mean = 1.0 and max = 2. The rejection ratio is
    // 3 rejections / 11 Newton iterations ≈ 0.2727.
    assert_eq!(tel.total_step_rejections(), 3);
    assert!((tel.mean_step_rejections().unwrap() - 1.0).abs() < 1e-12);
    assert_eq!(tel.max_step_rejections(), Some(2));
    assert!((tel.step_rejection_ratio().unwrap() - (3.0 / 11.0)).abs() < 1e-12);
}
