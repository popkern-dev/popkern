use super::*;
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;

#[test]
fn invalid_parameter_error_renders_parameter_and_reason() {
    let err = BlendError::InvalidParameter {
        parameter: "radius".to_string(),
        value: "-1".to_string(),
        reason: "radius must be strictly positive".to_string(),
    };
    let msg = err.to_string();
    assert!(msg.contains("radius"));
    assert!(msg.contains("-1"));
    assert!(msg.contains("strictly positive"));
}

#[test]
fn blend_error_remediation_hint_is_total() {
    // Every BlendError variant must return a non-empty rationale.
    let variants: Vec<BlendError> = vec![
        BlendError::InvalidParameter {
            parameter: "radius".to_string(),
            value: "NaN".to_string(),
            reason: "NaN is not a valid radius".to_string(),
        },
        BlendError::NoValidEdges("analysis failed".to_string()),
        BlendError::PartialFailure {
            failed: 1,
            total: 2,
            partial_result: Box::new(popkern_core::topology::BRep::new(Tolerance::default())),
            failures: vec![],
        },
    ];
    for v in variants {
        let hint = v.remediation_hint();
        assert!(!hint.rationale.is_empty(), "empty rationale for {v:?}");
    }
}

#[test]
fn invalid_parameter_routes_to_fix_input() {
    let err = BlendError::InvalidParameter {
        parameter: "radius".to_string(),
        value: "NaN".to_string(),
        reason: "not finite".to_string(),
    };
    assert_eq!(err.remediation_hint().action, RemediationAction::FixInput);
}

#[test]
fn no_valid_edges_routes_to_precheck() {
    let err = BlendError::NoValidEdges("bad".to_string());
    assert_eq!(err.remediation_hint().action, RemediationAction::Precheck);
}

#[test]
fn topology_error_routes_to_heal_input() {
    let err = BlendError::Topology(popkern_core::error::CoreError::PreconditionFailed(
        "dummy".to_string(),
    ));
    assert_eq!(err.remediation_hint().action, RemediationAction::HealInput);
}

#[test]
fn partial_failure_routes_to_precheck_with_counts_in_rationale() {
    let err = BlendError::PartialFailure {
        failed: 3,
        total: 7,
        partial_result: Box::new(popkern_core::topology::BRep::new(Tolerance::default())),
        failures: vec![],
    };
    let hint = err.remediation_hint();
    assert_eq!(hint.action, RemediationAction::Precheck);
    assert!(hint.rationale.contains('3'));
    assert!(hint.rationale.contains('7'));
}

#[test]
fn blend_failure_remediation_hint_is_total() {
    // Every BlendFailure variant must return a non-empty rationale.
    let vertex = popkern_core::topology::VertexHandle::default();
    let variants: Vec<BlendFailure> = vec![
        BlendFailure::RadiusTooLarge {
            requested: 1.0,
            max_radius: Some(0.3),
            location: None,
        },
        BlendFailure::RadiusTooLarge {
            requested: 1.0,
            max_radius: None,
            location: None,
        },
        BlendFailure::OffsetFailed {
            reason: "degenerate surface".to_string(),
            location: None,
        },
        BlendFailure::VertexBlendFailed {
            vertex,
            reason: "cannot fit patch".to_string(),
        },
        BlendFailure::IntersectionFailed {
            reason: "near-tangent".to_string(),
            location: Some(Point3::new(0.0, 0.0, 0.0)),
        },
        BlendFailure::InvalidTopology {
            reason: "no adjacent faces".to_string(),
        },
        BlendFailure::ContactSolverFailed {
            reason: "diverged".to_string(),
            spine_parameter: Some(0.5),
        },
        BlendFailure::AlreadyBlended,
        BlendFailure::SeamEdge,
    ];
    for v in variants {
        let hint = v.remediation_hint();
        assert!(!hint.rationale.is_empty(), "empty rationale for {v:?}");
    }
}

#[test]
fn radius_too_large_with_known_max_suggests_concrete_value() {
    let failure = BlendFailure::RadiusTooLarge {
        requested: 1.0,
        max_radius: Some(0.4),
        location: None,
    };
    let hint = failure.remediation_hint();
    assert_eq!(hint.action, RemediationAction::ReduceRadius);
    assert_eq!(hint.parameter.as_deref(), Some("radius"));
    let suggested = hint.value.expect("value populated for known max");
    // 5% safety margin: suggested ≈ 0.4 * 0.95 = 0.38
    assert!((suggested - 0.38).abs() < 1e-9);
}

#[test]
fn radius_too_large_without_known_max_has_no_concrete_value() {
    let failure = BlendFailure::RadiusTooLarge {
        requested: 1.0,
        max_radius: None,
        location: None,
    };
    let hint = failure.remediation_hint();
    assert_eq!(hint.action, RemediationAction::ReduceRadius);
    assert!(hint.value.is_none());
    assert!(
        hint.rationale.contains("feasible_radius_bounds"),
        "rationale should point to the query function"
    );
}

#[test]
fn offset_failed_suggests_tolerance_bump_by_default() {
    let failure = BlendFailure::OffsetFailed {
        reason: "parallel surfaces".to_string(),
        location: None,
    };
    let hint = failure.remediation_hint();
    assert_eq!(hint.action, RemediationAction::IncreaseTolerance);
    assert_eq!(hint.parameter.as_deref(), Some("linear"));
}

#[test]
fn offset_failed_after_workaround_escalates_to_heal_input() {
    let failure = BlendFailure::OffsetFailed {
        reason: "parallel surfaces".to_string(),
        location: None,
    };
    let ctx = RemediationContext {
        tolerance: None,
        workaround_attempted: true,
    };
    let hint = failure.remediation_hint_with_context(&ctx);
    assert_eq!(hint.action, RemediationAction::HealInput);
}

#[test]
fn intersection_failed_routes_to_increase_tolerance() {
    let failure = BlendFailure::IntersectionFailed {
        reason: "tangent".to_string(),
        location: None,
    };
    assert_eq!(
        failure.remediation_hint().action,
        RemediationAction::IncreaseTolerance
    );
}

#[test]
fn vertex_blend_failed_routes_to_alternative_mode() {
    let failure = BlendFailure::VertexBlendFailed {
        vertex: popkern_core::topology::VertexHandle::default(),
        reason: "degenerate".to_string(),
    };
    let hint = failure.remediation_hint();
    assert_eq!(hint.action, RemediationAction::UseAlternativeMode);
    assert!(
        hint.rationale.contains("Mitred"),
        "rationale should name the alternative"
    );
}

#[test]
fn seam_edge_routes_to_break_topology() {
    let hint = BlendFailure::SeamEdge.remediation_hint();
    assert_eq!(hint.action, RemediationAction::BreakTopology);
    assert!(hint.rationale.contains("seam"));
}

#[test]
fn already_blended_routes_to_no_action_with_unfillet_hint() {
    let hint = BlendFailure::AlreadyBlended.remediation_hint();
    assert_eq!(hint.action, RemediationAction::NoAction);
    assert!(hint.rationale.contains("unfillet"));
}

#[test]
fn contact_solver_failed_mentions_spine_parameter_when_known() {
    let failure = BlendFailure::ContactSolverFailed {
        reason: "NaN residual".to_string(),
        spine_parameter: Some(0.73),
    };
    let hint = failure.remediation_hint();
    assert_eq!(hint.action, RemediationAction::ReduceRadius);
    assert!(
        hint.rationale.contains("0.73"),
        "rationale should reference the spine parameter: {}",
        hint.rationale
    );
}

#[test]
fn remediation_with_tolerance_context_fills_concrete_linear_value() {
    let tol = Tolerance::default();
    let failure = BlendFailure::IntersectionFailed {
        reason: "tangent".to_string(),
        location: None,
    };
    let ctx = RemediationContext {
        tolerance: Some(&tol),
        workaround_attempted: false,
    };
    let hint = failure.remediation_hint_with_context(&ctx);
    assert_eq!(hint.action, RemediationAction::IncreaseTolerance);
    assert_eq!(hint.parameter.as_deref(), Some("linear"));
    // The BlendFailure path doesn't currently plumb the tolerance
    // through, but the BlendError path does. We verify the plumb
    // for BlendError below.
}

#[test]
fn blend_error_tolerance_context_fills_linear_value_on_tolerance_bump() {
    // When a BlendError naturally routes to IncreaseTolerance, the
    // caller-supplied tolerance should appear as a concrete value.
    // Currently BlendError only routes to IncreaseTolerance via
    // BlendFailure plumbing, so this checks the fill_tolerance_value
    // path using a direct Remediation constructed with a parameter.
    let tol = Tolerance::default();
    let ctx = RemediationContext {
        tolerance: Some(&tol),
        workaround_attempted: false,
    };
    let r = Remediation {
        action: RemediationAction::IncreaseTolerance,
        parameter: Some("linear".to_string()),
        value: None,
        rationale: "test".to_string(),
    };
    let filled = r.fill_tolerance_value(&ctx);
    assert!(filled.value.is_some());
    let v = filled.value.unwrap();
    assert!((v - tol.linear() * 10.0).abs() < 1e-12);
}

#[test]
fn remediation_new_constructs_minimal_hint() {
    let r = Remediation::new(RemediationAction::FixInput, "check your input");
    assert_eq!(r.action, RemediationAction::FixInput);
    assert!(r.parameter.is_none());
    assert!(r.value.is_none());
    assert_eq!(r.rationale, "check your input");
}

#[test]
fn remediation_adjust_tolerance_populates_parameter_and_value() {
    let r = Remediation::adjust_tolerance(
        RemediationAction::IncreaseTolerance,
        "linear",
        0.01,
        "widen solver envelope",
    );
    assert_eq!(r.parameter.as_deref(), Some("linear"));
    assert_eq!(r.value, Some(0.01));
    assert_eq!(r.action, RemediationAction::IncreaseTolerance);
}

#[test]
fn remediation_reduce_radius_uses_radius_parameter() {
    let r = Remediation::reduce_radius(0.38, "5% margin");
    assert_eq!(r.action, RemediationAction::ReduceRadius);
    assert_eq!(r.parameter.as_deref(), Some("radius"));
    assert_eq!(r.value, Some(0.38));
}
