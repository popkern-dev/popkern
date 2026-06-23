//! KP20 side-by-side macro-cap assembly wiring.
//!
//! This entry point wires the already-implemented KP20 pieces in the
//! paper order available today: guide-sampled macro patches, input
//! tensor-border split, internal C2 propagation, and S0/G2 strip
//! projection. It deliberately remains separate from the KP15
//! assembler and takes guide-sampled patches as input until the full
//! KP18 c-net guide builder is wired.

use super::kp_polar_boundary_input::TensorBorderBi3;
use super::kp20_bi6_cap::{
    Kp20Bi6Patch, Kp20MacroPatch, Kp20PieceId, initialise_kp20_macro_patches_from_tensor_border,
    kp20_propagate_c2_upper_from_lower,
};
use super::kp20_bi6_cnet::{Kp20RegularCNet, kp20_tensor_border_from_regular_cnet};
use super::kp20_bi6_guide_build::{Kp20GuideBuildInput, build_kp20_guide_ring};
use super::kp20_bi6_math::kp20_b_coeffs_near_o;
use super::kp20_bi6_piece_maps::sample_kp20_macro_patch_from_default_maps;
use super::kp20_bi6_s0_solve::{Kp20S0Projection, Kp20S0SolveError, project_kp20_g2_strip_pair};
use crate::geometry::Point3;

/// Input for end-to-end KP20 bi-6 cap assembly from tensor-border data.
#[derive(Debug, Clone)]
pub struct Kp20Bi6CapInput {
    /// Extraordinary point used as the guide-ring apex sample.
    pub star_point: Point3,
    /// Bi-3 tensor-border around the hole.
    pub tensor_border: TensorBorderBi3,
}

/// Input for KP20 bi-6 cap assembly from a regular cubic c-net ring.
#[derive(Debug, Clone)]
pub struct Kp20Bi6CNetCapInput {
    /// Extraordinary point used as the guide-ring apex sample.
    pub star_point: Point3,
    /// Regular cubic c-net stencils around the cap.
    pub cnet: Kp20RegularCNet,
}

/// Input for KP20 macro-cap assembly from guide-sampled patches.
#[derive(Debug, Clone)]
pub struct Kp20Bi6MacroCapInput {
    /// Bi-3 tensor-border around the hole.
    pub tensor_border: TensorBorderBi3,
    /// One guide-sampled 2x2 macro-patch per sector.
    pub sampled_macro_patches: Vec<Kp20MacroPatch>,
}

/// Output from KP20 macro-cap assembly.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20Bi6MacroCapOutput {
    /// Valence of the extraordinary vertex.
    pub n: usize,
    /// Assembled sector macro-patches after boundary, C2, and G2 passes.
    pub macro_patches: Vec<Kp20MacroPatch>,
    /// Flattened `4n` bi-6 patch grids in sector-major, [`Kp20PieceId::all`] order.
    pub patches_bi6: Vec<[[Point3; 7]; 7]>,
    /// Per-sector S0/G2 strip projection diagnostics.
    pub strip_projections: Vec<Kp20S0Projection>,
}

/// Error from KP20 macro-cap assembly.
#[derive(Debug, Clone, PartialEq)]
pub enum Kp20Bi6AssembleError {
    /// Tensor-border valence and sampled macro-patch count differ.
    ValenceMismatch {
        /// Number of tensor-border sectors.
        tensor_border: usize,
        /// Number of guide-sampled macro-patches supplied.
        sampled_macro_patches: usize,
    },
    /// Guide-ring construction failed.
    Guide(String),
    /// Tensor-border split or validation failed.
    TensorBorder(String),
    /// Regular c-net conversion failed.
    CNet(String),
    /// One sector's S0 strip solve failed.
    S0 {
        /// Sector whose strip projection failed.
        sector: usize,
        /// Underlying projection error.
        source: Kp20S0SolveError,
    },
}

impl std::fmt::Display for Kp20Bi6AssembleError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ValenceMismatch {
                tensor_border,
                sampled_macro_patches,
            } => write!(
                f,
                "KP20 macro-cap valence mismatch: tensor-border has {tensor_border}, sampled patches have {sampled_macro_patches}"
            ),
            Self::Guide(e) => write!(f, "KP20 guide-ring: {e}"),
            Self::TensorBorder(e) => write!(f, "KP20 tensor-border: {e}"),
            Self::CNet(e) => write!(f, "KP20 c-net: {e}"),
            Self::S0 { sector, source } => write!(f, "KP20 sector {sector} S0 solve: {source}"),
        }
    }
}

impl std::error::Error for Kp20Bi6AssembleError {}

/// Assemble a KP20 2x2 bi-6 macro cap from tensor-border data.
pub fn assemble_kp20_bi6_cap(
    input: &Kp20Bi6CapInput,
) -> Result<Kp20Bi6MacroCapOutput, Kp20Bi6AssembleError> {
    let guide = build_kp20_guide_ring(&Kp20GuideBuildInput {
        star_point: input.star_point,
        tensor_border: input.tensor_border.clone(),
    })
    .map_err(Kp20Bi6AssembleError::Guide)?;
    let sampled_macro_patches = guide
        .sectors
        .iter()
        .map(sample_kp20_macro_patch_from_default_maps)
        .collect();
    assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: input.tensor_border.clone(),
        sampled_macro_patches,
    })
}

/// Assemble a KP20 2x2 bi-6 macro cap from regular cubic c-net data.
pub fn assemble_kp20_bi6_cap_from_cnet(
    input: &Kp20Bi6CNetCapInput,
) -> Result<Kp20Bi6MacroCapOutput, Kp20Bi6AssembleError> {
    let tensor_border =
        kp20_tensor_border_from_regular_cnet(&input.cnet).map_err(Kp20Bi6AssembleError::CNet)?;
    assemble_kp20_bi6_cap(&Kp20Bi6CapInput {
        star_point: input.star_point,
        tensor_border,
    })
}

/// Assemble a KP20 2x2 bi-6 macro cap from guide-sampled patches.
pub fn assemble_kp20_bi6_macro_cap(
    input: &Kp20Bi6MacroCapInput,
) -> Result<Kp20Bi6MacroCapOutput, Kp20Bi6AssembleError> {
    let n = input.tensor_border.sectors.len();
    if input.sampled_macro_patches.len() != n {
        return Err(Kp20Bi6AssembleError::ValenceMismatch {
            tensor_border: n,
            sampled_macro_patches: input.sampled_macro_patches.len(),
        });
    }
    let mut macro_patches = input.sampled_macro_patches.clone();
    let boundary = initialise_kp20_macro_patches_from_tensor_border(&input.tensor_border)
        .map_err(Kp20Bi6AssembleError::TensorBorder)?;

    // Save the originally-sampled apex BB CPs of NearOriginLeft and
    // NearOriginRight (the ◦-marked region in KP20 Fig. 6c, indices
    // i + j ≤ 3 of each piece's bi-6 BB grid). These come from the
    // Hermite-sample of the LSQ guide which carries Lemma-4-coherent
    // structure (verified by probe
    // `lsq_guide_full_apex_jet_transpose_match_across_sectors`). The
    // C² propagation below overwrites these with NearMidpoint-derived
    // values that destroy that coherence; we restore them after step
    // (iii) C² runs, realising paper §4.4 step (i).
    let lemma4_overlay_left: Vec<[[Point3; 4]; 4]> = macro_patches
        .iter()
        .map(|mp| {
            let cps = &mp.piece(Kp20PieceId::NearOriginLeft).cps;
            let mut block = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
            for i in 0..4 {
                for j in 0..(4 - i) {
                    block[i][j] = cps[i][j];
                }
            }
            block
        })
        .collect();
    let lemma4_overlay_right: Vec<[[Point3; 4]; 4]> = macro_patches
        .iter()
        .map(|mp| {
            let cps = &mp.piece(Kp20PieceId::NearOriginRight).cps;
            let mut block = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
            for i in 0..4 {
                for j in 0..(4 - i) {
                    block[i][j] = cps[i][j];
                }
            }
            block
        })
        .collect();

    // Save FULL NearOrigin pieces too so we can restore the
    // C²-overwritten rows after step (i) overlay. Paper §4.4 does
    // not explicitly list C² propagation as a step — it's a kernel
    // optimisation that creates inconsistency with the overlay's
    // Lemma-4-coherent apex.
    let near_origin_left_full: Vec<[[Point3; 7]; 7]> = macro_patches
        .iter()
        .map(|mp| mp.piece(Kp20PieceId::NearOriginLeft).cps)
        .collect();
    let near_origin_right_full: Vec<[[Point3; 7]; 7]> = macro_patches
        .iter()
        .map(|mp| mp.piece(Kp20PieceId::NearOriginRight).cps)
        .collect();

    apply_tensor_border_rows(&mut macro_patches, &boundary);
    apply_internal_c2_to_all(&mut macro_patches);

    // Paper §4.4 step (i): restore the ◦-marked CPs from the
    // pre-C² Hermite sample. The LSQ guide's Lemma-4-coherent apex
    // jet (verified machine-zero) propagates through to the bi-6
    // macro-patch BB form via the identical kp20_corner_block_from_jet
    // conversion used at sampling time.
    apply_lemma4_cubic_jet_overlay(
        &mut macro_patches,
        &lemma4_overlay_left,
        &lemma4_overlay_right,
    );
    // Restore non-apex cps that the C² propagation overwrote.
    // Specifically NearOrigin.cps[5][*] (row 5) which is outside
    // the apex overlay region (i + j > 3 always at i = 5). Without
    // this restoration the v=0 row of NearOrigin is a Frankenstein
    // of overlay (i ≤ 3) and C² (i = 5), which the eq(8) polynomial
    // residual amplifies. We also restore cps[0][i] and cps[1][i]
    // for i > 3 (outside overlay reach) to keep those rows
    // polynomially consistent with the overlay region.
    for (k, macro_patch) in macro_patches.iter_mut().enumerate() {
        let l = macro_patch.piece_mut(Kp20PieceId::NearOriginLeft);
        for j in 4..7 {
            l.cps[0][j] = near_origin_left_full[k][0][j];
            l.cps[1][j] = near_origin_left_full[k][1][j];
        }
        for j in 0..7 {
            l.cps[5][j] = near_origin_left_full[k][5][j];
        }
        let r = macro_patch.piece_mut(Kp20PieceId::NearOriginRight);
        for j in 4..7 {
            r.cps[0][j] = near_origin_right_full[k][0][j];
            r.cps[1][j] = near_origin_right_full[k][1][j];
        }
        for j in 0..7 {
            r.cps[5][j] = near_origin_right_full[k][5][j];
        }
    }

    let strip_projections = apply_g2_strip_projection_to_all(&mut macro_patches)?;
    let patches_bi6 = flatten_macro_patches(&macro_patches);

    Ok(Kp20Bi6MacroCapOutput {
        n,
        macro_patches,
        patches_bi6,
        strip_projections,
    })
}

fn apply_tensor_border_rows(macro_patches: &mut [Kp20MacroPatch], boundary: &[Kp20MacroPatch]) {
    for (dst, src) in macro_patches.iter_mut().zip(boundary) {
        dst.piece_mut(Kp20PieceId::NearMidpointLeft).cps[0..3]
            .copy_from_slice(&src.piece(Kp20PieceId::NearMidpointLeft).cps[0..3]);
        dst.piece_mut(Kp20PieceId::NearMidpointRight).cps[0..3]
            .copy_from_slice(&src.piece(Kp20PieceId::NearMidpointRight).cps[0..3]);
    }
}

fn apply_internal_c2_to_all(macro_patches: &mut [Kp20MacroPatch]) {
    for macro_patch in macro_patches {
        let lower_left = *macro_patch.piece(Kp20PieceId::NearMidpointLeft);
        let upper_left = macro_patch.piece_mut(Kp20PieceId::NearOriginLeft);
        kp20_propagate_c2_upper_from_lower(&lower_left, upper_left);

        let lower_right = *macro_patch.piece(Kp20PieceId::NearMidpointRight);
        let upper_right = macro_patch.piece_mut(Kp20PieceId::NearOriginRight);
        kp20_propagate_c2_upper_from_lower(&lower_right, upper_right);
    }
}

fn apply_g2_strip_projection_to_all(
    macro_patches: &mut [Kp20MacroPatch],
) -> Result<Vec<Kp20S0Projection>, Kp20Bi6AssembleError> {
    let n = macro_patches.len();
    let c = (core::f64::consts::TAU / n as f64).cos();
    let b0 = kp20_b_coeffs_near_o(c)[0];
    let mut out = Vec::with_capacity(n);
    for (sector, macro_patch) in macro_patches.iter_mut().enumerate() {
        let (left, right) = paired_origin_pieces_mut(macro_patch);
        let projection = project_kp20_g2_strip_pair(left, right, b0)
            .map_err(|source| Kp20Bi6AssembleError::S0 { sector, source })?;
        out.push(projection);
    }
    Ok(out)
}

fn paired_origin_pieces_mut(
    macro_patch: &mut Kp20MacroPatch,
) -> (&mut Kp20Bi6Patch, &mut Kp20Bi6Patch) {
    let [left, right, _, _] = &mut macro_patch.pieces;
    (left, right)
}

/// Paper §4.4 step (i): overlay the Lemma 4 ◦-marked CPs onto the
/// NearOrigin pieces' apex region (i + j ≤ 3) so the cubic
/// Hermite jet is preserved after C² propagation.
fn apply_lemma4_cubic_jet_overlay(
    macro_patches: &mut [Kp20MacroPatch],
    overlay_left: &[[[Point3; 4]; 4]],
    overlay_right: &[[[Point3; 4]; 4]],
) {
    for (k, macro_patch) in macro_patches.iter_mut().enumerate() {
        let block_l = &overlay_left[k];
        let block_r = &overlay_right[k];
        for i in 0..4 {
            for j in 0..(4 - i) {
                macro_patch.piece_mut(Kp20PieceId::NearOriginLeft).cps[i][j] = block_l[i][j];
                macro_patch.piece_mut(Kp20PieceId::NearOriginRight).cps[i][j] = block_r[i][j];
            }
        }
    }
}

fn flatten_macro_patches(macro_patches: &[Kp20MacroPatch]) -> Vec<[[Point3; 7]; 7]> {
    let mut out = Vec::with_capacity(4 * macro_patches.len());
    for macro_patch in macro_patches {
        for id in Kp20PieceId::all() {
            out.push(macro_patch.piece(id).cps);
        }
    }
    out
}

#[cfg(test)]
#[path = "kp20_bi6_assemble_tests.rs"]
mod tests;
