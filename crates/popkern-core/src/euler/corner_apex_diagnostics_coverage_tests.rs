use super::coverage::{
    CornerApexCoverageFacetPayload, CornerApexCoverageRow, MutatorReadinessRow,
    validate_corner_apex_coverage_row, validate_mutator_readiness_row,
};

fn complete_row() -> CornerApexCoverageRow {
    CornerApexCoverageRow {
        row_id: "corner-apex:three-face-tangency:non-planar-support",
        terminal_label: "three_face_tangency_non_planar_support",
        primitive_family: Some("three-face-tangent-patch"),
        invalidity_reason: Some("not-applicable-ready-row"),
        tolerance_source: Some("local-geometry-predicate"),
        candidate_family: Some("corner-apex-three-face"),
        facet_payload: Some(CornerApexCoverageFacetPayload::complete("covered")),
    }
}

fn remove_primitive_family(row: &mut CornerApexCoverageRow) {
    row.primitive_family = None;
}

fn remove_invalidity_reason(row: &mut CornerApexCoverageRow) {
    row.invalidity_reason = None;
}

fn remove_tolerance_source(row: &mut CornerApexCoverageRow) {
    row.tolerance_source = None;
}

fn remove_candidate_family(row: &mut CornerApexCoverageRow) {
    row.candidate_family = None;
}

fn remove_facet_payload(row: &mut CornerApexCoverageRow) {
    row.facet_payload = None;
}

fn remove_topology(payload: &mut CornerApexCoverageFacetPayload) {
    payload.topology = None;
}

fn remove_geometry(payload: &mut CornerApexCoverageFacetPayload) {
    payload.geometry = None;
}

fn remove_dependency(payload: &mut CornerApexCoverageFacetPayload) {
    payload.dependency = None;
}

fn remove_repair(payload: &mut CornerApexCoverageFacetPayload) {
    payload.repair = None;
}

fn remove_semantic_history(payload: &mut CornerApexCoverageFacetPayload) {
    payload.semantic_history = None;
}

fn remove_downstream(payload: &mut CornerApexCoverageFacetPayload) {
    payload.downstream = None;
}

#[test]
fn coverage_gate_rejects_generic_terminal_labels() {
    for terminal_label in [
        "missing_operator",
        "geometry_failed",
        "invalid_topology",
        "unsupported_case",
        "three-face-tangency-nonplanar-support",
    ] {
        let mut row = complete_row();
        row.terminal_label = terminal_label;

        let error = validate_corner_apex_coverage_row(&row).unwrap_err();

        assert_eq!(error.field, "terminal_label");
    }
}

#[test]
fn coverage_gate_names_missing_top_level_facets() {
    let missing_fields: [(&str, fn(&mut CornerApexCoverageRow)); 5] = [
        ("primitive_family", remove_primitive_family),
        ("invalidity_reason", remove_invalidity_reason),
        ("tolerance_source", remove_tolerance_source),
        ("candidate_family", remove_candidate_family),
        ("facet_payload", remove_facet_payload),
    ];

    for (expected_field, remove_field) in missing_fields {
        let mut row = complete_row();
        remove_field(&mut row);

        let error = validate_corner_apex_coverage_row(&row).unwrap_err();

        assert_eq!(error.field, expected_field);
    }
}

#[test]
fn coverage_gate_names_missing_six_facet_payload_members() {
    let missing_fields: [(&str, fn(&mut CornerApexCoverageFacetPayload)); 6] = [
        ("facet_payload.topology", remove_topology),
        ("facet_payload.geometry", remove_geometry),
        ("facet_payload.dependency", remove_dependency),
        ("facet_payload.repair", remove_repair),
        ("facet_payload.semantic_history", remove_semantic_history),
        ("facet_payload.downstream", remove_downstream),
    ];

    for (expected_field, remove_field) in missing_fields {
        let mut row = complete_row();
        let mut payload = row.facet_payload.expect("complete fixture has payload");
        remove_field(&mut payload);
        row.facet_payload = Some(payload);

        let error = validate_corner_apex_coverage_row(&row).unwrap_err();

        assert_eq!(error.field, expected_field);
    }
}

#[test]
fn mutator_readiness_requires_diagnostic_reference_or_closure_task() {
    assert!(
        validate_mutator_readiness_row(&MutatorReadinessRow {
            row_id: "mutator:split-cylinder-resource-merge",
            diagnostic_row_ref: Some("corner-apex:split-cylinder-resource-merge"),
            closure_task_id: None,
        })
        .is_ok()
    );
    assert!(
        validate_mutator_readiness_row(&MutatorReadinessRow {
            row_id: "mutator:boundary-side-separator",
            diagnostic_row_ref: None,
            closure_task_id: Some("PK-005"),
        })
        .is_ok()
    );

    let error = validate_mutator_readiness_row(&MutatorReadinessRow {
        row_id: "mutator:untracked",
        diagnostic_row_ref: None,
        closure_task_id: None,
    })
    .unwrap_err();

    assert_eq!(error.field, "diagnostic_row_ref_or_closure_task_id");
}
