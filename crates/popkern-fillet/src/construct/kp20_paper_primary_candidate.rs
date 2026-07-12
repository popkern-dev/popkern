//! Paper-faithful KP20 cap as the primary polar apex candidate.
//!
//! `assemble_kp20_paper_macro_cap` (popkern-core) builds the §4.4
//! cap exactly as printed — derived characteristic border, solved τ,
//! g∘χ̃-fitted guide, exact Hermite sampling, eq (14) border
//! overwrite, and the S-strip projection — and reports raw
//! eq (8)/(9) G¹/G² residuals plus C⁰/border/apex gates. This module
//! adapts it into the fillet candidate pipeline: geometric fold and
//! self-intersection checks reuse the polar fold sampler, while
//! acceptance gates on the paper report's own residuals (scaled to
//! the candidate's geometric size) instead of the legacy macro-cap
//! validator, whose piece chart and c/m identification encode the
//! superseded pre-recovery reading.

use super::super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, gregory_patches_from_macro_patches, summarize_candidate_intersections,
};
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20PaperCapReport, TensorBorderBi3, assemble_kp20_paper_macro_cap,
};

/// Relative residual gate for the paper cap's continuity residuals:
/// eq (8)/(9), C⁰, border, and apex rows must close to this fraction
/// of the candidate's geometric radius.
const PAPER_CAP_RELATIVE_RESIDUAL_GATE: f64 = 1.0e-6;

/// Relative gate on the INPUT border's corner-band disagreement at
/// the sector corners `c` (`input_corner_mismatch`), as a fraction of
/// the candidate's geometric radius — a G¹-class incompatibility of
/// the input that the assembly averages into the rim.
///
/// Why this is NOT the 1e-6 closure gate above: the mismatch is an
/// input-shape-class measurement (like `guide_fit_residual`), not an
/// assembly residual, and it is nonzero even on by-construction-valid
/// input — measured 1.7e-1 relative on the C² B-spline ring border
/// (fold-free, all continuity rows at machine zero) and 1.3e-1 on the
/// accepted pentagonal-pyramid G² setback candidate. Every folding
/// real-rolling-ball candidate measured to date sits at >= 2.3e-1
/// relative (hexagon6_variable setback/blend-boundary, pentagon5 and
/// hexagon6 midpoint-split adapters, PK-090 calibration). The gate is
/// placed in the measured separation band so corner-incompatible
/// input declines attributably (measured value in the detail string)
/// without regressing the accepted class.
const PAPER_CAP_RELATIVE_CORNER_MISMATCH_GATE: f64 = 2.0e-1;

/// Commercial-tier relative gate for the derivative rows (eq (8)/(9)
/// G¹/G² residuals and the Lemma-4 delta). Real rolling-ball apex data
/// is not in the paper-exact input class (its corner mismatch alone is
/// 1e-1-class), so a fold-free candidate whose derivative rows close
/// within this band is commercially far better than the mitred
/// fallback it would otherwise become. Calibration (PK-090 candidate
/// table): the accepted class measures g2 ≤ 4.2e-5 relative
/// (pentagon5_mixed setback: folds=0, g1 = 9.9e-14, g2 = 1.979e-6 abs
/// at scale 4.8e-2); every folding candidate measures ≥ 6.4 relative.
/// Positional rows (C⁰/border) keep the paper-exact gate — the tier
/// relaxes smoothness bands, never watertightness.
const PAPER_CAP_COMMERCIAL_RELATIVE_DERIVATIVE_GATE: f64 = 1.0e-3;

/// How a paper-border candidate cleared acceptance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum Kp20PaperAcceptanceTier {
    /// Every paper gate closed at the 1e-6 relative form.
    PaperExact,
    /// Fold-free and corner-compatible with positional rows paper-exact,
    /// derivative rows within the commercial band; carries the measured
    /// relative g2 residual for the relaxed-continuity record.
    Commercial { g2_relative: f64 },
}

/// Paper-cap candidate with its acceptance evidence.
pub(super) struct Kp20PaperPrimaryCandidate {
    /// Bi-9 Gregory patches, sector-major piece order.
    pub(super) patches: Vec<GregoryPatch>,
    /// Geometric fold / self-intersection summary.
    pub(super) summary: Kp20SamplerFoldSummary,
    /// Paper assembly report (raw G¹/G² gates).
    pub(super) report: Kp20PaperCapReport,
    /// Geometric radius used to scale the residual gates.
    pub(super) scale: f64,
}

impl Kp20PaperPrimaryCandidate {
    /// Accept when the paper gates close at the relative tolerance,
    /// the surface is regular, and the fold sampler agrees.
    pub(super) fn is_acceptable(&self) -> bool {
        let scale = self.scale.max(f64::MIN_POSITIVE);
        let gate = PAPER_CAP_RELATIVE_RESIDUAL_GATE * scale;
        // Corner-incompatible input must DECLINE (with the measured
        // `input_corner_mismatch` in the detail string), not ship a
        // kinked rim; see the corner gate constant for the measured
        // calibration and why its relative form differs from the
        // closure rows.
        let corner_gate = PAPER_CAP_RELATIVE_CORNER_MISMATCH_GATE * scale;
        self.report.max_g1_residual < gate
            && self.report.max_g2_residual < gate
            && self.report.max_c0_residual < gate
            && self.report.max_border_residual < gate
            && self.report.lemma4_delta < gate
            && self.report.input_corner_mismatch < corner_gate
            && self.report.min_regularity > 0.0
            && self.summary.fold_count == 0
    }

    /// Acceptance tier: paper-exact, commercial, or `None` (decline).
    /// The commercial tier keeps every positional/fold/corner
    /// requirement and relaxes only the derivative rows to the
    /// measured commercial band (see the gate constant).
    pub(super) fn acceptance_tier(&self) -> Option<Kp20PaperAcceptanceTier> {
        if self.is_acceptable() {
            return Some(Kp20PaperAcceptanceTier::PaperExact);
        }
        let scale = self.scale.max(f64::MIN_POSITIVE);
        let exact_gate = PAPER_CAP_RELATIVE_RESIDUAL_GATE * scale;
        let derivative_gate = PAPER_CAP_COMMERCIAL_RELATIVE_DERIVATIVE_GATE * scale;
        let corner_gate = PAPER_CAP_RELATIVE_CORNER_MISMATCH_GATE * scale;
        let commercial = self.summary.fold_count == 0
            && self.report.min_regularity > 0.0
            && self.report.input_corner_mismatch < corner_gate
            && self.report.max_c0_residual < exact_gate
            && self.report.max_border_residual < exact_gate
            && self.report.max_g1_residual < derivative_gate
            && self.report.max_g2_residual < derivative_gate
            && self.report.lemma4_delta < derivative_gate;
        commercial.then(|| Kp20PaperAcceptanceTier::Commercial {
            g2_relative: self.report.max_g2_residual / scale,
        })
    }

    /// Structured diagnostic detail (paper residuals + fold summary).
    /// `corner_mismatch` is the measured `input_corner_mismatch`, so
    /// corner-gate declines are attributable from the record alone.
    pub(super) fn detail(&self) -> String {
        format!(
            "map_family=paper_section_4_4;status=assembled;folds={};g1={:.3e};g2={:.3e};c0={:.3e};border={:.3e};lemma4={:.3e};corner_mismatch={:.3e};guide_fit={:.3e};strip_move={:.3e};min_regularity={:.3e};scale={:.3e}",
            self.summary.fold_count,
            self.report.max_g1_residual,
            self.report.max_g2_residual,
            self.report.max_c0_residual,
            self.report.max_border_residual,
            self.report.lemma4_delta,
            self.report.input_corner_mismatch,
            self.report.guide_fit_residual,
            self.report.max_strip_movement,
            self.report.min_regularity,
            self.scale,
        )
    }
}

/// Record a paper candidate's outcome under `label` and return the
/// accepted patches. Shared contract for every paper-border primary
/// candidate (setback, blend-boundary): an ACCEPTED candidate records
/// `label` with its full detail; a gate decline records
/// `{label}_declined` with the same detail (fold count and measured
/// `input_corner_mismatch` included); a construction error records
/// `{label}_error` with the error string.
pub(super) fn record_paper_candidate_outcome(
    diagnostics: &mut Vec<crate::types::ApexConstructionDiagnostic>,
    label: &str,
    n: usize,
    candidate: Result<Kp20PaperPrimaryCandidate, String>,
) -> Option<(Vec<GregoryPatch>, Kp20PaperAcceptanceTier)> {
    use super::super::kp20_polar_diagnostics::dispatch_fact;
    match candidate {
        Ok(candidate) => match candidate.acceptance_tier() {
            Some(tier @ Kp20PaperAcceptanceTier::PaperExact) => {
                diagnostics.push(dispatch_fact(
                    label,
                    n,
                    Some(candidate.detail()),
                    None,
                    None,
                ));
                Some((candidate.patches, tier))
            }
            Some(tier @ Kp20PaperAcceptanceTier::Commercial { .. }) => {
                diagnostics.push(dispatch_fact(
                    &format!("{label}_commercial"),
                    n,
                    Some(format!("tier=commercial;{}", candidate.detail())),
                    None,
                    None,
                ));
                Some((candidate.patches, tier))
            }
            None => {
                diagnostics.push(dispatch_fact(
                    &format!("{label}_declined"),
                    n,
                    Some(candidate.detail()),
                    None,
                    None,
                ));
                None
            }
        },
        Err(e) => {
            diagnostics.push(dispatch_fact(
                &format!("{label}_error"),
                n,
                Some(e),
                None,
                None,
            ));
            None
        }
    }
}

/// Assemble the paper cap for the apex and package it as a polar
/// candidate.
pub(super) fn try_kp20_paper_primary_candidate(
    star_point: Point3,
    tensor_border: &TensorBorderBi3,
) -> Result<Kp20PaperPrimaryCandidate, String> {
    let (output, report) =
        assemble_kp20_paper_macro_cap(star_point, tensor_border).map_err(|e| e.to_string())?;
    let patches = gregory_patches_from_macro_patches(&output.macro_patches);
    let intersections = patches
        .iter()
        .enumerate()
        .filter_map(|(idx, gp)| {
            let report = gp.self_intersection_report(SELF_INT_TOL);
            report.intersects().then_some((idx, report))
        })
        .collect::<Vec<_>>();
    let summary = summarize_candidate_intersections(&patches, &intersections);
    let scale = tensor_border
        .sectors
        .iter()
        .flat_map(|rows| rows.iter().flatten())
        .map(|cp| cp.distance_to(&star_point))
        .fold(0.0_f64, f64::max);
    Ok(Kp20PaperPrimaryCandidate {
        patches,
        summary,
        report,
        scale,
    })
}

/// Valid production-shaped border for tests: each depth row family
/// is a uniform periodic cubic B-spline ring (C² across the spokes by
/// construction), one cubic per sector, non-planar depth profile.
#[cfg(test)]
pub(crate) fn bspline_ring_border(n: usize) -> TensorBorderBi3 {
    let mut sectors = vec![[[Point3::new(0.0, 0.0, 0.0); 4]; 3]; n];
    for depth in 0..3_usize {
        let radius = 1.0 - 0.08 * depth as f64;
        let height = 0.12 * depth as f64;
        let control = |k: usize| -> Point3 {
            let a = core::f64::consts::TAU * (k % n) as f64 / n as f64;
            Point3::new(radius * a.cos(), radius * a.sin(), height)
        };
        for (k, sector) in sectors.iter_mut().enumerate() {
            let p = [
                control(k + n - 1),
                control(k),
                control(k + 1),
                control(k + 2),
            ];
            let combine = |w: [f64; 4]| {
                Point3::new(
                    w[0] * p[0].x + w[1] * p[1].x + w[2] * p[2].x + w[3] * p[3].x,
                    w[0] * p[0].y + w[1] * p[1].y + w[2] * p[2].y + w[3] * p[3].y,
                    w[0] * p[0].z + w[1] * p[1].z + w[2] * p[2].z + w[3] * p[3].z,
                )
            };
            sector[depth] = [
                combine([1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0, 0.0]),
                combine([0.0, 2.0 / 3.0, 1.0 / 3.0, 0.0]),
                combine([0.0, 1.0 / 3.0, 2.0 / 3.0, 0.0]),
                combine([0.0, 1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0]),
            ];
        }
    }
    TensorBorderBi3 { sectors }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// End-to-end production gate: the paper cap must be assembled
    /// AND accepted (G² at relative tolerance, fold-free) for a valid
    /// C² production border across the supported polar valences.
    #[test]
    fn paper_primary_candidate_is_accepted_on_valid_borders() {
        for n in [5, 6, 7, 8] {
            let border = bspline_ring_border(n);
            let star = Point3::new(0.0, 0.0, 0.3);
            let candidate =
                try_kp20_paper_primary_candidate(star, &border).expect("paper candidate");
            eprintln!("paper candidate n={n}: {}", candidate.detail());
            assert!(
                candidate.is_acceptable(),
                "n={n}: paper candidate rejected: {}",
                candidate.detail(),
            );
            assert_eq!(candidate.patches.len(), 4 * n);
        }
    }
}
