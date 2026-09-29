//! The setback cap's Gregory patches from the KP20 assembly (PK-599):
//! the constrained bi-9 resampling when it served, else the elevated
//! bi-6 cap (geometrically the paper's cap). The decline reason is on
//! the report; the `setback` trace channel names the outcome.

use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20PaperCapReport, RectangularGregoryBiNonic,
};

/// The cap's patches, sector-major, pieces in `Kp20PieceId::all` order.
pub(super) fn setback_cap_patches(report: &Kp20PaperCapReport) -> Vec<GregoryPatch> {
    if popkern_core::trace::trace_enabled("setback") {
        match &report.bi9 {
            Ok(r) => eprintln!(
                "  [bi9] resampled: relations {:?}, movement {:.3e}",
                r.relation_residual, r.deviation_after
            ),
            Err(e) => eprintln!("  [bi9] declined, bi-6 serves: {e}"),
        }
    }
    match &report.bi9 {
        Ok(resample) if !crate::dispatch_overrides::kp20_bi9_resample_disabled() => resample
            .pieces
            .iter()
            .map(|grid| {
                GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic {
                    grid: *grid,
                }))
            })
            .collect(),
        _ => super::super::kp20_polar_fold_summary::gregory_patches_from_macro_patches(
            &report.sectors,
        ),
    }
}
