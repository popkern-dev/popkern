//! KP20 §4.4 paper-order cap assembly (WP4).
//!
//! Assembles the n-sector 2×2 bi-6 macro cap exactly in the paper's
//! order, from the WP2 τ and the WP3 guide:
//!
//! - step (iii): every coefficient from the exact Fig-9 Hermite
//!   sampling of `h := g ∘ τ` (`kp20_paper_cap_sample`); the ◦-marked
//!   apex coefficients of step (i) are automatically the Lemma-4
//!   cubic-expansion data because the sampling is exact (consistency
//!   is measured, not overlaid);
//! - step (ii): the gray rows along each sector's rim overwritten
//!   with the eq (14)-reparameterized split input tensor-border
//!   (the WP2-validated layout, shared placement builder);
//! - step (iv): per sector-separating curve, the strip set S is
//!   overwritten by the closest values (to the sampled f) on the
//!   affine manifold cut out by the raw §4.1/§4.2 equations — C⁰,
//!   the degree-5 rows (11), the seam C² rows (12), and the full
//!   eq (8)/(9) G¹/G² rows — which is exactly the manifold the paper
//!   parameterizes by S₀.
//!
//! Chart: the WP2 quadrant convention per sector (P1 owns `o`; sector
//! `s` spans spokes `s`, `s+1`; spoke `s` is sector `s`'s `a = 0`
//! edge and sector `s−1`'s `b = 0` edge). The curve between sectors
//! `s−1` and `s` is flanked by (P1, P2) of sector `s−1` (natural
//! frame for the eq (8)/(9) Jacobians) and (P1, P4) of sector `s`
//! (transposed indices).

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use super::kp20_bi6_cap::{Kp20MacroPatch, Kp20PieceId};
use super::kp20_paper_cap_sample::kp20_sample_paper_cap_sector;
use super::kp20_paper_cap_strip::{Kp20StripEnforcement, project_strip};
use super::kp20_paper_guide::Kp20PaperGuideBorder;
use super::kp20_paper_tau::{kp20_paper_border_band_placements, solve_kp20_paper_tau};
use crate::geometry::Point3;
use popkern_math::PortableMath;

/// Assembled paper cap with measured postconditions.
#[derive(Debug, Clone)]
pub struct Kp20PaperCapReport {
    /// Valence.
    pub n: usize,
    /// One 2×2 macro patch per sector (WP2 quadrant chart).
    pub sectors: Vec<Kp20MacroPatch>,
    /// Guide fit residual (KP18 LSQ misfit; a modeling measurement,
    /// not an error).
    pub guide_fit_residual: f64,
    /// Cross-sector C⁰ mismatch along the sector-separating curves.
    pub max_c0_residual: f64,
    /// Worst eq (8) G¹ residual BB coefficient across all curves.
    pub max_g1_residual: f64,
    /// Worst eq (9) G² residual BB coefficient across all curves.
    pub max_g2_residual: f64,
    /// Worst deviation of the gray rows from the eq (14) bands.
    pub max_border_residual: f64,
    /// Worst apex cross-sector consistency delta (orders 0..=3 of the
    /// Lemma-4 edge runs).
    pub lemma4_delta: f64,
    /// Largest movement step (iv) applied to any strip coefficient.
    pub max_strip_movement: f64,
    /// Corner-compatibility of the INPUT border at the sector
    /// corners `c`: the worst disagreement between each sector's two
    /// eq (14) corner bands before averaging. Near zero only for
    /// valid two-piece input (Hermite data of one cornered patch);
    /// large values mean the assembled rim will kink at `c`.
    pub input_corner_mismatch: f64,
    /// Minimum sampled surface-regularity measure `|∂u × ∂v|` over
    /// all pieces.
    pub min_regularity: f64,
    /// The constrained bi-9 resampling of the cap (PK-599), or why it
    /// declined — the bi-6 `sectors` then serve.
    pub bi9:
        Result<super::kp20_bi9_resample::Kp20Bi9Resample, super::kp20_bi9_resample::Kp20Bi9Decline>,
}

/// Errors from the paper cap assembly.
#[derive(Debug)]
pub enum Kp20PaperCapError {
    /// Guide fit failed.
    Guide(super::kp20_bi6_guide_lsq::Kp20GuideLsqError),
    /// τ solve failed.
    Tau(super::kp20_bi6_tau_solve_assemble::Kp20PlanarTauPaperSolveError),
    /// The step (iv) strip system was singular.
    StripSolve(usize),
}

impl std::fmt::Display for Kp20PaperCapError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Guide(err) => write!(f, "paper cap guide fit failed: {err}"),
            Self::Tau(err) => write!(f, "paper cap tau solve failed: {err}"),
            Self::StripSolve(curve) => {
                write!(f, "paper cap strip solve singular at curve {curve}")
            }
        }
    }
}

impl std::error::Error for Kp20PaperCapError {}

/// Production adapter: assemble the paper cap from a smooth
/// single-cubic-per-sector tensor-border (the kernel's production
/// input shape) and return it in the macro-cap output
/// contract consumed by the fillet dispatch (`macro_patches` +
/// flattened `patches_bi6`; no per-strip diagnostics — the paper
/// path's own gates live in [`Kp20PaperCapReport`]).
pub fn assemble_kp20_paper_macro_cap(
    star_point: Point3,
    tensor_border: &super::kp_polar_boundary_input::TensorBorderBi3,
) -> Result<
    (
        super::kp20_bi6_assemble::Kp20Bi6MacroCapOutput,
        Kp20PaperCapReport,
    ),
    Kp20PaperCapError,
> {
    let border = Kp20PaperGuideBorder::from_tensor_border(tensor_border);
    let report = assemble_kp20_paper_cap(&border, star_point)?;
    let patches_bi6 = report
        .sectors
        .iter()
        .flat_map(|sector| {
            Kp20PieceId::all()
                .into_iter()
                .map(|piece_id| sector.piece(piece_id).cps)
        })
        .collect();
    let output = super::kp20_bi6_assemble::Kp20Bi6MacroCapOutput {
        n: report.n,
        macro_patches: report.sectors.clone(),
        patches_bi6,
        strip_projections: Vec::new(),
    };
    Ok((output, report))
}

/// Assemble the KP20 §4.4 paper cap for a two-piece-per-sector input
/// tensor-border.
pub fn assemble_kp20_paper_cap(
    border: &Kp20PaperGuideBorder,
    star_point: Point3,
) -> Result<Kp20PaperCapReport, Kp20PaperCapError> {
    assemble_kp20_paper_cap_with_enforcement(border, star_point, Kp20StripEnforcement::HardG2)
}

/// Assemble the KP20 §4.4 paper cap for a two-piece-per-sector input
/// tensor-border with an explicit strip-enforcement mode.
///
/// `HardG2` is the paper's construction and requires input in the
/// paper's smoothness class (C² ring data). `SoftG` is recorded
/// containment for out-of-class input (e.g., G¹ rolling-ball fillet
/// boundaries): the sector-separating strips keep C⁰/(11)/(12) as
/// hard constraints and minimize the eq (8)/(9) residuals in least
/// squares instead of enforcing an infeasible system, trading exact
/// G² for boundedness. The report's residuals tell the truth either
/// way.
pub fn assemble_kp20_paper_cap_with_enforcement(
    border: &Kp20PaperGuideBorder,
    star_point: Point3,
    enforcement: Kp20StripEnforcement,
) -> Result<Kp20PaperCapReport, Kp20PaperCapError> {
    assemble_kp20_paper_cap_conditioned(border, star_point, enforcement, 0.0)
}

/// Full-control assembly: strip enforcement plus the recorded guide
/// fairing weight (see `solve_kp20_paper_guide_with_fairing`).
pub fn assemble_kp20_paper_cap_conditioned(
    border: &Kp20PaperGuideBorder,
    star_point: Point3,
    enforcement: Kp20StripEnforcement,
    guide_fairing: f64,
) -> Result<Kp20PaperCapReport, Kp20PaperCapError> {
    let n = border.sectors.len();
    let guide = super::kp20_paper_guide::solve_kp20_paper_guide_with_fairing(
        border,
        star_point,
        guide_fairing,
    )
    .map_err(Kp20PaperCapError::Guide)?;
    let tau = solve_kp20_paper_tau(n).map_err(Kp20PaperCapError::Tau)?;

    // Step (iii): exact Hermite sampling per sector, with τ rotated
    // into each sector.
    let alpha = core::f64::consts::TAU / n as f64;
    let mut sectors: Vec<Kp20MacroPatch> = (0..n)
        .map(|s| {
            let rot = [(s as f64 * alpha).pcos(), (s as f64 * alpha).psin()];
            let mut tau_sector = tau.tau;
            for piece_id in Kp20PieceId::all() {
                let piece = tau_sector.piece_mut(piece_id);
                for i in 0..7 {
                    for j in 0..7 {
                        let cp = piece.cps[i][j];
                        piece.cps[i][j] = Point3::new(
                            cp.x * rot[0] - cp.y * rot[1],
                            cp.x * rot[1] + cp.y * rot[0],
                            0.0,
                        );
                    }
                }
            }
            kp20_sample_paper_cap_sector(n, s, &guide.sectors[s], &tau_sector)
        })
        .collect();
    let sampled = sectors.clone();

    // Step (ii): overwrite the gray rows from the reparameterized
    // split input border.
    let c = alpha.pcos();
    let mut placements_per_sector = Vec::with_capacity(n);
    let mut input_corner_mismatch = 0.0_f64;
    for sector_border in border.sectors.iter() {
        input_corner_mismatch =
            input_corner_mismatch.max(super::kp20_paper_tau::kp20_paper_border_corner_mismatch(
                n,
                &sector_border[0].rows,
                &sector_border[1].rows,
                2.0 * c / 3.0,
            ));
        let placements = kp20_paper_border_band_placements(
            n,
            &sector_border[0].rows,
            &sector_border[1].rows,
            2.0 * c / 3.0,
        );
        placements_per_sector.push(placements);
    }

    // Cross-sector seam reconciliation BEFORE pinning: the two
    // sectors flanking each separating curve both band-place the
    // shared seam-row slots (prev P2 (i,0) ↔ curr P4 (0,i) near m,
    // prev P1 (i,0) ↔ curr P1 (0,i) near o), each from ITS adjacent
    // input border piece. On real (non-C²-at-m) data the two claims
    // disagree by the m cross-row mismatch (measured 1.6e-3..4e-3 on
    // the pyramid apex), and step (iv) cannot repair it — both slots
    // are band-frozen, so the C⁰ row degenerates to a constant and is
    // dropped. Averaging the claims makes the seam C⁰-exact by
    // construction; the strip projection then re-solves the G rows
    // around the reconciled pins.
    {
        use super::kp20_bi6_cap::Kp20PieceId;
        const P2: Kp20PieceId = Kp20PieceId::NearOriginRight;
        const P4: Kp20PieceId = Kp20PieceId::NearMidpointLeft;
        for curve in 0..n {
            let prev = (curve + n - 1) % n;
            let curr = curve;
            for i in 4..=6_usize {
                {
                    let (prev_piece, curr_piece) = (P2, P4);
                    let a = placements_per_sector[prev]
                        .iter()
                        .position(|&(p, pi, pj, _)| p == prev_piece && pi == i && pj == 0);
                    let b = placements_per_sector[curr]
                        .iter()
                        .position(|&(p, pi, pj, _)| p == curr_piece && pi == 0 && pj == i);
                    if let (Some(a), Some(b)) = (a, b) {
                        let va = placements_per_sector[prev][a].3;
                        let vb = placements_per_sector[curr][b].3;
                        let mean = Point3::new(
                            0.5 * (va.x + vb.x),
                            0.5 * (va.y + vb.y),
                            0.5 * (va.z + vb.z),
                        );
                        placements_per_sector[prev][a].3 = mean;
                        placements_per_sector[curr][b].3 = mean;
                    }
                }
            }
        }
    }
    for (s, placements) in placements_per_sector.iter().enumerate() {
        for &(piece, i, j, value) in placements {
            sectors[s].piece_mut(piece).cps[i][j] = value;
        }
    }

    // Step (iv): strip projection per sector-separating curve.
    let mut max_strip_movement = 0.0_f64;
    for curve in 0..n {
        let movement = project_strip(&mut sectors, n, curve, enforcement)
            .ok_or(Kp20PaperCapError::StripSolve(curve))?;
        max_strip_movement = max_strip_movement.max(movement);
    }

    let mut report = super::kp20_paper_cap_measure::measure_paper_cap(
        n,
        &sectors,
        &placements_per_sector,
        guide.max_fit_residual,
        max_strip_movement,
    );
    report.input_corner_mismatch = input_corner_mismatch;
    report.bi9 = match enforcement {
        Kp20StripEnforcement::HardG2 => super::kp20_bi9_resample::resample_kp20_cap_bi9(
            n,
            &guide.sectors,
            &tau.tau,
            &sampled,
            &report.sectors,
        ),
        Kp20StripEnforcement::SoftG => {
            Err(super::kp20_bi9_resample::Kp20Bi9Decline::SoftEnforcement)
        }
    };
    Ok(report)
}

#[cfg(test)]
#[path = "kp20_paper_cap_tests.rs"]
mod tests;
