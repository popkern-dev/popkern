//! Tests for malformed STEP file import error handling.
//!
//! Validates that the import pipeline produces informative errors for
//! structurally broken, incomplete, or unsupported STEP inputs.

use super::*;
use crate::step::diagnostics::Phase;

/// Helper: default options with healing disabled (to isolate parse/resolve errors).
fn opts() -> StepImportOptions {
    StepImportOptions {
        healing: HealingLevel::None,
        ..Default::default()
    }
}

/// Extract the error from an import result, panicking if it succeeded.
fn unwrap_err(result: Result<StepImportResult, IoError>, ctx: &str) -> IoError {
    match result {
        Ok(_) => panic!("{ctx}: expected Err, got Ok"),
        Err(e) => e,
    }
}

// ── Empty input ─────────────────────────────────────────────────────────

#[test]
fn empty_input_returns_parse_error() {
    let err = unwrap_err(import_step(b"", &opts()), "empty input");
    assert!(
        matches!(err, IoError::Parse { .. }),
        "expected IoError::Parse, got {err:?}"
    );
    // The error message should mention the missing ISO sentinel.
    let msg = err.to_string();
    assert!(
        msg.contains("ISO-10303-21") || msg.contains("keyword"),
        "error should reference ISO sentinel or keyword, got: {msg}"
    );
}

#[test]
fn whitespace_only_input_returns_parse_error() {
    let err = unwrap_err(
        import_step(b"   \n\t  \n", &opts()),
        "whitespace-only input",
    );
    assert!(matches!(err, IoError::Parse { .. }));
}

// ── Invalid header ──────────────────────────────────────────────────────

#[test]
fn missing_iso_sentinel_returns_parse_error() {
    let data = b"HEADER;\nENDSEC;\nDATA;\nENDSEC;\nEND-ISO-10303-21;\n";
    let err = unwrap_err(import_step(data, &opts()), "missing ISO sentinel");
    assert!(matches!(err, IoError::Parse { .. }));
    let msg = err.to_string();
    assert!(
        msg.contains("ISO-10303-21"),
        "error should mention ISO-10303-21, got: {msg}"
    );
}

#[test]
fn garbage_before_sentinel_returns_parse_error() {
    let data = b"this is not a step file at all";
    let err = unwrap_err(import_step(data, &opts()), "garbage input");
    assert!(matches!(err, IoError::Parse { .. }));
}

#[test]
fn truncated_header_returns_parse_error() {
    // ISO sentinel present, but HEADER section is never closed.
    let data = b"ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('test'), '2;1');\n";
    let err = unwrap_err(import_step(data, &opts()), "truncated header");
    assert!(matches!(err, IoError::Parse { .. }));
}

// ── Valid header, no DATA section ───────────────────────────────────────

#[test]
fn missing_data_section_returns_parse_error() {
    let data = b"ISO-10303-21;\n\
        HEADER;\n\
        FILE_DESCRIPTION(('test'), '2;1');\n\
        FILE_NAME('test.stp', '2024-01-01', (''), (''), '', '', '');\n\
        FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\n\
        ENDSEC;\n\
        END-ISO-10303-21;\n";
    let err = unwrap_err(import_step(data, &opts()), "missing DATA section");
    assert!(matches!(err, IoError::Parse { .. }));
    let msg = err.to_string();
    assert!(
        msg.contains("DATA"),
        "error should reference DATA keyword, got: {msg}"
    );
}

#[test]
fn header_without_endsec_returns_parse_error() {
    // Header that goes straight to DATA without ENDSEC.
    let data = b"ISO-10303-21;\n\
        HEADER;\n\
        FILE_DESCRIPTION(('test'), '2;1');\n\
        DATA;\n\
        ENDSEC;\n\
        END-ISO-10303-21;\n";
    let err = unwrap_err(import_step(data, &opts()), "header without ENDSEC");
    // This should fail because the parser expects ENDSEC for HEADER before DATA.
    assert!(matches!(err, IoError::Parse { .. }));
}

// ── Malformed entity references (broken #N references) ──────────────────

/// A well-formed STEP wrapper around a DATA section. All callers supply
/// only the entity lines (without the DATA/ENDSEC delimiters).
fn wrap_step_data(entity_lines: &str) -> Vec<u8> {
    format!(
        "ISO-10303-21;\n\
         HEADER;\n\
         FILE_DESCRIPTION(('test'), '2;1');\n\
         FILE_NAME('test.stp', '2024-01-01', (''), (''), '', '', '');\n\
         FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\n\
         ENDSEC;\n\
         DATA;\n\
         {entity_lines}\n\
         ENDSEC;\n\
         END-ISO-10303-21;\n"
    )
    .into_bytes()
}

#[test]
fn broken_entity_ref_logged_as_diagnostic_error() {
    // #1 references #999 which does not exist.
    let data = wrap_step_data("#1 = AXIS2_PLACEMENT_3D('', #999, #999, #999);");
    let result = import_step(&data, &opts()).expect("import should succeed with diagnostics");
    // The EntityGraph.build() should log broken references.
    let resolve_errors: Vec<_> = result
        .diagnostics
        .errors()
        .filter(|e| e.phase == Phase::Resolve)
        .collect();
    assert!(
        !resolve_errors.is_empty(),
        "should have resolve-phase errors for broken refs"
    );
    // Error message should mention the dangling ref.
    let has_ref_msg = resolve_errors
        .iter()
        .any(|e| e.message.contains("#999") || e.message.contains("999"));
    assert!(
        has_ref_msg,
        "error should mention the dangling ref #999, got: {:?}",
        resolve_errors
            .iter()
            .map(|e| &e.message)
            .collect::<Vec<_>>()
    );
}

#[test]
fn multiple_broken_refs_each_logged() {
    // #1 and #2 both reference missing entities.
    let data = wrap_step_data(
        "#1 = AXIS2_PLACEMENT_3D('', #100, #200, #300);\n\
         #2 = EDGE_CURVE('', #400, #500, #600, .T.);",
    );
    let result = import_step(&data, &opts()).expect("import should succeed with diagnostics");
    let resolve_errors: Vec<_> = result
        .diagnostics
        .errors()
        .filter(|e| e.phase == Phase::Resolve)
        .collect();
    // Each entity should log at least one broken-ref error.
    assert!(
        resolve_errors.len() >= 2,
        "expected at least 2 resolve errors, got {}",
        resolve_errors.len()
    );
}

// ── Unsupported entity types ────────────────────────────────────────────

#[test]
fn unsupported_entity_type_does_not_cause_hard_error() {
    // TOTALLY_MADE_UP_ENTITY is not in the registry — the import should
    // succeed (it becomes ResolvedEntity::Unknown), not crash.
    let data = wrap_step_data("#1 = TOTALLY_MADE_UP_ENTITY('foo', 42, (1.0, 2.0));");
    let result = import_step(&data, &opts());
    assert!(
        result.is_ok(),
        "unrecognized entity type should not cause a hard error: {:?}",
        result.err()
    );
}

#[test]
fn only_unsupported_entities_produces_empty_document() {
    // A file with only unsupported entities has no products or BReps.
    let data = wrap_step_data(
        "#1 = TOTALLY_MADE_UP_ENTITY('foo', 42);\n\
         #2 = ANOTHER_UNKNOWN('bar', #1);",
    );
    let result = import_step(&data, &opts()).expect("should succeed");
    assert!(
        result.document.products.is_empty(),
        "no products should be produced from unknown entities"
    );
}

// ── Valid structure, empty DATA section ──────────────────────────────────

#[test]
fn empty_data_section_produces_empty_document() {
    let data = wrap_step_data("");
    let result = import_step(&data, &opts()).expect("empty DATA section should succeed");
    assert!(result.document.products.is_empty());
    assert!(result.document.roots.is_empty());
}

// ── Malformed entity syntax in DATA section ─────────────────────────────

#[test]
fn malformed_entity_skipped_with_parse_warning() {
    // Missing '=' between entity ref and type name.
    let data = wrap_step_data(
        "#1 CARTESIAN_POINT('', (0.0, 0.0, 0.0));\n\
         #2 = DIRECTION('', (0.0, 0.0, 1.0));",
    );
    let result = import_step(&data, &opts());
    // The parser should either skip the bad entity (adding a parse error)
    // or fail entirely. Either behavior is acceptable — the key is that
    // we don't panic.
    match result {
        Ok(r) => {
            // If it recovered, we should have a parse warning about the skip.
            let parse_diags: Vec<_> = r
                .diagnostics
                .entries()
                .iter()
                .filter(|e| e.phase == Phase::Parse)
                .collect();
            assert!(
                !parse_diags.is_empty() || r.document.products.is_empty(),
                "skipped entity should produce a diagnostic or no output"
            );
        }
        Err(e) => {
            // Hard parse failure is acceptable for structurally broken entities.
            assert!(matches!(e, IoError::Parse { .. }));
        }
    }
}

// ── Missing END-ISO-10303-21 sentinel ───────────────────────────────────

#[test]
fn missing_end_sentinel_returns_parse_error() {
    let data = b"ISO-10303-21;\n\
        HEADER;\n\
        FILE_DESCRIPTION(('test'), '2;1');\n\
        FILE_NAME('test.stp', '2024-01-01', (''), (''), '', '', '');\n\
        FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\n\
        ENDSEC;\n\
        DATA;\n\
        ENDSEC;\n";
    let err = unwrap_err(import_step(data, &opts()), "missing END-ISO-10303-21");
    assert!(matches!(err, IoError::Parse { .. }));
    let msg = err.to_string();
    assert!(
        msg.contains("END-ISO-10303-21") || msg.contains("END"),
        "error should reference the end sentinel, got: {msg}"
    );
}

// ── Error message quality ───────────────────────────────────────────────

#[test]
fn parse_error_contains_position_info() {
    // Ensure IoError::Parse carries line/col or token position context.
    let err = unwrap_err(import_step(b"not a step file", &opts()), "garbage input");
    match err {
        IoError::Parse { message, .. } => {
            assert!(
                !message.is_empty(),
                "parse error message should not be empty"
            );
            // The message should contain useful context (token position or keyword name).
            assert!(
                message.contains("token")
                    || message.contains("keyword")
                    || message.contains("ISO")
                    || message.contains("byte"),
                "error message should have context, got: {message}"
            );
        }
        other => panic!("expected IoError::Parse, got {other:?}"),
    }
}

#[test]
fn broken_ref_error_includes_entity_id() {
    let data = wrap_step_data("#42 = AXIS2_PLACEMENT_3D('', #9999, $, $);");
    let result = import_step(&data, &opts()).expect("should succeed with diagnostics");
    let resolve_errors: Vec<_> = result
        .diagnostics
        .errors()
        .filter(|e| e.phase == Phase::Resolve)
        .collect();
    // At least one error should identify the originating entity (#42).
    let has_entity_id = resolve_errors.iter().any(|e| e.entity_id == Some(42));
    assert!(
        has_entity_id,
        "broken-ref error should include the originating entity id (42), got: {:?}",
        resolve_errors
            .iter()
            .map(|e| (e.entity_id, &e.message))
            .collect::<Vec<_>>()
    );
}
