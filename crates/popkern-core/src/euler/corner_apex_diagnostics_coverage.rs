//! Coverage-gate validation for corner/apex diagnostic rows.
#![allow(missing_docs)]

const GENERIC_TERMINAL_LABELS: [&str; 5] = [
    "missing_operator",
    "geometry_failed",
    "invalid_topology",
    "unsupported_case",
    "three-face-tangency-nonplanar-support",
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CornerApexCoverageFacetPayload {
    pub topology: Option<&'static str>,
    pub geometry: Option<&'static str>,
    pub dependency: Option<&'static str>,
    pub repair: Option<&'static str>,
    pub semantic_history: Option<&'static str>,
    pub downstream: Option<&'static str>,
}

impl CornerApexCoverageFacetPayload {
    pub const fn complete(value: &'static str) -> Self {
        Self {
            topology: Some(value),
            geometry: Some(value),
            dependency: Some(value),
            repair: Some(value),
            semantic_history: Some(value),
            downstream: Some(value),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CornerApexCoverageRow {
    pub row_id: &'static str,
    pub terminal_label: &'static str,
    pub primitive_family: Option<&'static str>,
    pub invalidity_reason: Option<&'static str>,
    pub tolerance_source: Option<&'static str>,
    pub candidate_family: Option<&'static str>,
    pub facet_payload: Option<CornerApexCoverageFacetPayload>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MutatorReadinessRow {
    pub row_id: &'static str,
    pub diagnostic_row_ref: Option<&'static str>,
    pub closure_task_id: Option<&'static str>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CornerApexCoverageValidationError {
    pub row_id: &'static str,
    pub field: &'static str,
}

pub fn validate_corner_apex_coverage_row(
    row: &CornerApexCoverageRow,
) -> Result<(), CornerApexCoverageValidationError> {
    if GENERIC_TERMINAL_LABELS.contains(&row.terminal_label) {
        return Err(missing(row.row_id, "terminal_label"));
    }

    require(row.row_id, row.primitive_family, "primitive_family")?;
    require(row.row_id, row.invalidity_reason, "invalidity_reason")?;
    require(row.row_id, row.tolerance_source, "tolerance_source")?;
    require(row.row_id, row.candidate_family, "candidate_family")?;

    let payload = row
        .facet_payload
        .ok_or_else(|| missing(row.row_id, "facet_payload"))?;
    require(row.row_id, payload.topology, "facet_payload.topology")?;
    require(row.row_id, payload.geometry, "facet_payload.geometry")?;
    require(row.row_id, payload.dependency, "facet_payload.dependency")?;
    require(row.row_id, payload.repair, "facet_payload.repair")?;
    require(
        row.row_id,
        payload.semantic_history,
        "facet_payload.semantic_history",
    )?;
    require(row.row_id, payload.downstream, "facet_payload.downstream")?;

    Ok(())
}

pub fn validate_mutator_readiness_row(
    row: &MutatorReadinessRow,
) -> Result<(), CornerApexCoverageValidationError> {
    if row.diagnostic_row_ref.is_some() || row.closure_task_id.is_some() {
        return Ok(());
    }

    Err(missing(row.row_id, "diagnostic_row_ref_or_closure_task_id"))
}

fn require(
    row_id: &'static str,
    value: Option<&'static str>,
    field: &'static str,
) -> Result<(), CornerApexCoverageValidationError> {
    value.map(|_| ()).ok_or_else(|| missing(row_id, field))
}

const fn missing(row_id: &'static str, field: &'static str) -> CornerApexCoverageValidationError {
    CornerApexCoverageValidationError { row_id, field }
}
