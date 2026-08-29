use super::CornerApexMatrixRow;

#[path = "apex_fillet_diagnostic_matrix_corner_rows_a.rs"]
mod apex_fillet_diagnostic_matrix_corner_rows_a;
#[path = "apex_fillet_diagnostic_matrix_corner_rows_b.rs"]
mod apex_fillet_diagnostic_matrix_corner_rows_b;

pub(crate) fn corner_apex_rows() -> Vec<CornerApexMatrixRow> {
    let mut rows = Vec::with_capacity(44);
    rows.extend_from_slice(&apex_fillet_diagnostic_matrix_corner_rows_a::CORNER_APEX_ROWS_A);
    rows.extend_from_slice(&apex_fillet_diagnostic_matrix_corner_rows_b::CORNER_APEX_ROWS_B);
    rows
}
