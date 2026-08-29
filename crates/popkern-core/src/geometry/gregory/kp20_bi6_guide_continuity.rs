//! KP18 total-degree-6 guide sector continuity formulas.
//!
//! KP20 reuses KP18's degree-6 guide; this module keeps the KP18
//! equation (1) join weights. The full KP18 (2)-(5) join helpers
//! that once lived here were deleted: they addressed BB CPs with a
//! different indexing convention than
//! `solve_kp20_guide_lsq::constraint_rows` (the geometrically
//! correct, production encoding — it identifies CPs along the
//! shared spoke between adjacent sectors), no production caller
//! ever invoked them, and keeping a wrong transcription alive as
//! "informational" invited misuse. See
//! `docs/decisions/2026-06-23-kp20-bi6-paper-grade-guide-from-kp18.md`
//! Step 4 for the audit.

/// KP18 equation (1) join weights for adjacent sector wedges.
pub fn kp18_guide_join_weights(n: usize) -> Result<(f64, f64, f64), String> {
    if n < 3 {
        return Err(format!(
            "KP18/KP20 guide join weights require n >= 3; got {n}"
        ));
    }
    let c = (core::f64::consts::TAU / n as f64).cos();
    Ok((-1.0, 2.0 * c, 2.0 * (1.0 - c)))
}

#[cfg(test)]
#[path = "kp20_bi6_guide_continuity_tests.rs"]
mod tests;
