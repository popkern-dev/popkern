//! PR08 Catmull-Clark characteristic-ring coefficient substrate.
//!
//! # ⚠ SUPERSEDED — transcribes a book misprint (recovery plan RC5)
//!
//! This module's PR08 §6.1 subdivision blocks transcribe the printed
//! row for `q'₄` (PR08 p.112), which carries an UNPUBLISHED misprint:
//! the printed row is `(q₂δ, q₁+q₂w, q₁ | …)` where the true
//! Catmull-Clark edge-midpoint stencil requires `(q₁δ, q₁+q₂w, q₂ | …)`.
//! The misprint preserves row sums and the eigenvector here is solved
//! from the same wrong matrix, so the module's internal gates (affine
//! invariance, eigen-relation, λ-scaling, rotational covariance) all
//! pass while every produced ring value `q₄..q₁₃` is wrong at every
//! valence. All tensor-border source candidates ranked through this
//! module are therefore void. See
//! `docs/architecture/kp20-paper-faithful-recovery-plan.md` (RC5) and
//! the independent verification in `tools/reference/cc-characteristic/`.
//!
//! The reference-grade replacement is `cc_characteristic_ring`, which
//! derives the subdivision rows from mesh-connectivity stencils and is
//! pinned against the external oracle. This module is retained only so
//! the historical diagnostic probes still compile until the WP6
//! triage; do not add new consumers.
//!
//! KP20 needs the tensor-border of the Catmull-Clark characteristic
//! ring, not a direct row fill from the six KP15/KP18 cap-net labels.
//! This module starts from the PR08 §6.1 block-circulant
//! Catmull-Clark subdivision matrix and computes the Fourier-index-one
//! eigenvector for the 13 control coefficients of one segment.

use super::kp_polar_boundary_input::TensorBorderBi3;
use super::kp20_bi6_cnet::{
    Kp20RegularTensorBorderCNet, kp20_tensor_border_from_regular_cnet_rows,
};
use super::kp20_bi6_math::Kp20CharacteristicScalars;
use crate::geometry::{Point3, PortableMath};

#[derive(Debug, Clone, PartialEq)]
/// PR08 frequency-one Catmull-Clark characteristic ring.
pub struct Kp20CatmullClarkCharacteristicRing {
    /// One PR08 Fig. 6.3 segment per sector, with coefficients
    /// `q_{j,1}..q_{j,13}`.
    pub sectors: Vec<[Point3; 13]>,
    /// Catmull-Clark subdominant eigenvalue.
    pub lambda: f64,
}

/// Compute the PR08 frequency-one Catmull-Clark characteristic ring.
pub fn kp20_catmull_clark_characteristic_ring(
    n: usize,
) -> Result<Kp20CatmullClarkCharacteristicRing, String> {
    let scalars = Kp20CharacteristicScalars::for_valence(n)?;
    let v_hat = pr08_catmull_clark_fourier_one_eigenvector(n)?;
    let mut sectors = Vec::with_capacity(n);
    for j in 0..n {
        let theta = core::f64::consts::TAU * j as f64 / n as f64;
        let wj = Complex::from_angle(theta);
        let mut sector = [Point3::new(0.0, 0.0, 0.0); 13];
        for (k, cp) in sector.iter_mut().enumerate() {
            let value = wj * v_hat[k];
            *cp = Point3::new(value.re, value.im, 0.0);
        }
        sectors.push(sector);
    }
    Ok(Kp20CatmullClarkCharacteristicRing {
        sectors,
        lambda: scalars.lambda,
    })
}

/// Build the accepted PR08/Stam Appendix-A Q1 characteristic tensor-border.
///
/// This is the production source-construction surface for the KP20
/// characteristic tensor-border campaign. It evaluates PR08's
/// frequency-one Catmull-Clark characteristic ring, extracts Stam
/// Appendix-B `q^1` regular subpatch interior rows from the Appendix-A
/// subdivided characteristic labels, then converts those regular cubic
/// B-spline rows to the bi-3 tensor-border consumed by KP20 equation (14).
pub fn kp20_pr08_stam_q1_characteristic_tensor_border(n: usize) -> Result<TensorBorderBi3, String> {
    let ring = kp20_catmull_clark_characteristic_ring(n)?;
    let mut sectors = Vec::with_capacity(n);
    for sector in 0..n {
        sectors.push(stam_appendix_a_regular_subpatch_interior_rows(
            ring.lambda,
            n,
            sector,
            0,
        )?);
    }
    let rows = Kp20RegularTensorBorderCNet { sectors };
    kp20_tensor_border_from_regular_cnet_rows(&rows)
}

/// Build PR08/Stam Appendix-A/B `q^1` regular tensor-border rows from the
/// PR08 frequency-one inner ring.
///
/// This follows the source correspondence now used by the KP20 Fig. 8
/// campaign: PR08 `q_{j,1}` and `q_{j,2}` seed Stam Appendix-A's top vector,
/// Stam Appendix-A subdivision produces the extended characteristic labels,
/// and Stam Appendix-B `q^1` selects the regular subpatch grid whose interior
/// three rows feed the KP20 tensor-border conversion.
pub fn kp20_pr08_stam_inner_q1_regular_tensor_border_rows(
    n: usize,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let ring = kp20_catmull_clark_characteristic_ring(n)?;
    let mut sectors = Vec::with_capacity(n);
    for sector in 0..n {
        sectors.push(stam_appendix_a_pr08_inner_ring_subpatch_rows(
            &ring, sector, 0,
        )?);
    }
    Ok(Kp20RegularTensorBorderCNet { sectors })
}

/// Build the PR08 Fig. 6.3 endpoint-side regular rows.
///
/// The rows are the side of the Fig. 6.3 lower-right regular B-spline
/// patch whose endpoints are the two neighboring `q13` controls:
/// `[q13, q11, q8, prev q13]`, `[q12, q10, q5, prev q4]`,
/// `[q7, q6, q2, prev q3]`. This is the paper-derived endpoint-edge
/// source for KP15 Fig. 3's `v_{k-1} -> m_k` tensor-border experiments;
/// it is intentionally separate from the inner `q1` source whose first
/// endpoint is the replicated extraordinary center.
#[cfg(test)]
pub fn kp20_pr08_fig63_endpoint_side_regular_tensor_border_rows(
    n: usize,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let ring = kp20_catmull_clark_characteristic_ring(n)?;
    let mut sectors = Vec::with_capacity(n);
    for sector_index in 0..n {
        let sector = &ring.sectors[sector_index];
        let previous = &ring.sectors[(sector_index + n - 1) % n];
        sectors.push([
            [sector[12], sector[10], sector[7], previous[12]],
            [sector[11], sector[9], sector[4], previous[3]],
            [sector[6], sector[5], sector[1], previous[2]],
        ]);
    }
    Ok(Kp20RegularTensorBorderCNet { sectors })
}

/// Build the endpoint-normalized PR08/Stam inner-ring `q^1` regular rows.
///
/// The normalization is a global complex similarity applied to all regular
/// rows so sector `0`'s nonzero first-row endpoint matches the unit sector
/// spoke order used by the KP20 τ equations. The opposite first-row endpoint
/// is the extraordinary origin for this source and remains the origin. This is
/// recorded as a separate API because it is a chart normalization of the
/// PR08/Stam source, not a different source construction.
pub fn kp20_pr08_stam_inner_q1_endpoint_normalized_regular_tensor_border_rows(
    n: usize,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let rows = kp20_pr08_stam_inner_q1_regular_tensor_border_rows(n)?;
    endpoint_normalized_regular_tensor_border_rows(&rows, n, 0, true, true)
}

/// Convert the endpoint-normalized PR08/Stam inner-ring `q^1` regular rows to
/// the bi-3 tensor-border form consumed by KP20 equation (14).
pub fn kp20_pr08_stam_inner_q1_characteristic_tensor_border(
    n: usize,
) -> Result<TensorBorderBi3, String> {
    let rows = kp20_pr08_stam_inner_q1_endpoint_normalized_regular_tensor_border_rows(n)?;
    kp20_tensor_border_from_regular_cnet_rows(&rows)
}

#[cfg(test)]
// Index form: `delta` is simultaneously the circulant block address and the
// DFT phase exponent of the PR08 block algebra; iterators would obscure that.
#[allow(clippy::needless_range_loop)]
pub fn pr08_catmull_clark_subdivision_matrix(n: usize) -> Result<Vec<Vec<f64>>, String> {
    if n < 3 {
        return Err(format!(
            "PR08 Catmull-Clark subdivision matrix needs n >= 3; got {n}"
        ));
    }
    let block_count = 13;
    let size = n * block_count;
    let mut blocks = vec![[[Complex::zero(); 13]; 13]; n];
    for mode in 0..n {
        let a_hat = pr08_catmull_clark_a_hat_block(n, mode);
        for delta in 0..n {
            let phase = Complex::from_angle(
                -core::f64::consts::TAU * mode as f64 * delta as f64 / n as f64,
            );
            for row in 0..block_count {
                for col in 0..block_count {
                    blocks[delta][row][col] = blocks[delta][row][col] + a_hat[row][col] * phase;
                }
            }
        }
    }

    let mut matrix = vec![vec![0.0_f64; size]; size];
    for dest_sector in 0..n {
        for delta in 0..n {
            let src_sector = (dest_sector + delta) % n;
            for row in 0..block_count {
                for col in 0..block_count {
                    matrix[dest_sector * block_count + row][src_sector * block_count + col] =
                        blocks[delta][row][col].re / n as f64;
                }
            }
        }
    }
    Ok(matrix)
}

#[cfg(test)]
pub fn pr08_catmull_clark_subdivision_stencil_matrix(n: usize) -> Result<Vec<Vec<f64>>, String> {
    if n < 3 {
        return Err(format!(
            "PR08 Catmull-Clark subdivision stencil matrix needs n >= 3; got {n}"
        ));
    }
    let size = n * 13;
    let mut matrix = vec![vec![0.0_f64; size]; size];

    fn add(
        matrix: &mut [Vec<f64>],
        n: usize,
        dest_sector: usize,
        dest_label: usize,
        source_delta: usize,
        source_label: usize,
        weight: f64,
    ) {
        if weight == 0.0 {
            return;
        }
        let source_sector = (dest_sector + source_delta) % n;
        matrix[dest_sector * 13 + dest_label][source_sector * 13 + source_label] += weight;
    }

    fn add_zero_mode(
        matrix: &mut [Vec<f64>],
        n: usize,
        dest_sector: usize,
        dest_label: usize,
        source_label: usize,
        weight: f64,
    ) {
        for source_delta in 0..n {
            add(
                matrix,
                n,
                dest_sector,
                dest_label,
                source_delta,
                source_label,
                weight / n as f64,
            );
        }
    }

    fn add_same(
        matrix: &mut [Vec<f64>],
        n: usize,
        dest_sector: usize,
        dest_label: usize,
        source_label: usize,
        weight: f64,
    ) {
        add(matrix, n, dest_sector, dest_label, 0, source_label, weight);
    }

    fn add_next(
        matrix: &mut [Vec<f64>],
        n: usize,
        dest_sector: usize,
        dest_label: usize,
        source_label: usize,
        weight: f64,
    ) {
        add(matrix, n, dest_sector, dest_label, 1, source_label, weight);
    }

    fn add_previous(
        matrix: &mut [Vec<f64>],
        n: usize,
        dest_sector: usize,
        dest_label: usize,
        source_label: usize,
        weight: f64,
    ) {
        add(
            matrix,
            n,
            dest_sector,
            dest_label,
            n - 1,
            source_label,
            weight,
        );
    }

    fn add_cosine(
        matrix: &mut [Vec<f64>],
        n: usize,
        dest_sector: usize,
        dest_label: usize,
        source_label: usize,
        weight: f64,
    ) {
        add_next(
            matrix,
            n,
            dest_sector,
            dest_label,
            source_label,
            0.5 * weight,
        );
        add_previous(
            matrix,
            n,
            dest_sector,
            dest_label,
            source_label,
            0.5 * weight,
        );
    }

    fn add_one_plus_next(
        matrix: &mut [Vec<f64>],
        n: usize,
        dest_sector: usize,
        dest_label: usize,
        source_label: usize,
        weight: f64,
    ) {
        add_same(matrix, n, dest_sector, dest_label, source_label, weight);
        add_next(matrix, n, dest_sector, dest_label, source_label, weight);
    }

    fn add_one_plus_previous(
        matrix: &mut [Vec<f64>],
        n: usize,
        dest_sector: usize,
        dest_label: usize,
        source_label: usize,
        weight: f64,
    ) {
        add_same(matrix, n, dest_sector, dest_label, source_label, weight);
        add_previous(matrix, n, dest_sector, dest_label, source_label, weight);
    }

    let alpha = 1.0 - 7.0 / (4.0 * n as f64);
    let beta = 3.0 / (2.0 * n as f64);
    let gamma = 1.0 / (4.0 * n as f64);
    let p1 = 1.0 / 64.0;
    let p2 = 3.0 / 32.0;
    let p3 = 9.0 / 16.0;
    let q1 = 1.0 / 16.0;
    let q2 = 3.0 / 8.0;
    let r = 1.0 / 4.0;

    for sector in 0..n {
        add_zero_mode(&mut matrix, n, sector, 0, 0, alpha);
        add_zero_mode(&mut matrix, n, sector, 0, 1, beta);
        add_zero_mode(&mut matrix, n, sector, 0, 2, gamma);

        add_zero_mode(&mut matrix, n, sector, 1, 0, q2);
        add_cosine(&mut matrix, n, sector, 1, 1, 2.0 * q1);
        add_same(&mut matrix, n, sector, 1, 1, q2);
        add_one_plus_next(&mut matrix, n, sector, 1, 2, q1);

        add_zero_mode(&mut matrix, n, sector, 2, 0, r);
        add_one_plus_previous(&mut matrix, n, sector, 2, 1, r);
        add_same(&mut matrix, n, sector, 2, 2, r);

        add_zero_mode(&mut matrix, n, sector, 3, 0, q2);
        add_same(&mut matrix, n, sector, 3, 1, q1);
        add_next(&mut matrix, n, sector, 3, 1, q2);
        add_same(&mut matrix, n, sector, 3, 2, q1);
        add_same(&mut matrix, n, sector, 3, 3, q1);
        add_next(&mut matrix, n, sector, 3, 4, q1);

        add_zero_mode(&mut matrix, n, sector, 4, 0, p2);
        add_cosine(&mut matrix, n, sector, 4, 1, 2.0 * p1);
        add_same(&mut matrix, n, sector, 4, 1, p3);
        add_one_plus_previous(&mut matrix, n, sector, 4, 2, p2);
        add_previous(&mut matrix, n, sector, 4, 3, p1);
        add_same(&mut matrix, n, sector, 4, 4, p2);
        add_same(&mut matrix, n, sector, 4, 5, p1);

        add_zero_mode(&mut matrix, n, sector, 5, 0, q1);
        add_next(&mut matrix, n, sector, 5, 1, q1);
        add_same(&mut matrix, n, sector, 5, 1, q2);
        add_same(&mut matrix, n, sector, 5, 2, q2);
        add_same(&mut matrix, n, sector, 5, 4, q1);
        add_same(&mut matrix, n, sector, 5, 5, q1);

        add_zero_mode(&mut matrix, n, sector, 6, 0, p1);
        add_one_plus_next(&mut matrix, n, sector, 6, 1, p2);
        add_same(&mut matrix, n, sector, 6, 2, p3);
        add_same(&mut matrix, n, sector, 6, 3, p2);
        add_one_plus_next(&mut matrix, n, sector, 6, 4, p1);
        add_same(&mut matrix, n, sector, 6, 5, p2);
        add_same(&mut matrix, n, sector, 6, 6, p1);

        add_same(&mut matrix, n, sector, 7, 1, q2);
        add_one_plus_previous(&mut matrix, n, sector, 7, 2, q1);
        add_previous(&mut matrix, n, sector, 7, 3, q1);
        add_same(&mut matrix, n, sector, 7, 4, q2);
        add_same(&mut matrix, n, sector, 7, 5, q1);

        add_same(&mut matrix, n, sector, 8, 1, r);
        add_same(&mut matrix, n, sector, 8, 2, r);
        add_same(&mut matrix, n, sector, 8, 4, r);
        add_same(&mut matrix, n, sector, 8, 5, r);

        add_same(&mut matrix, n, sector, 9, 1, q1);
        add_same(&mut matrix, n, sector, 9, 2, q2);
        add_same(&mut matrix, n, sector, 9, 3, q1);
        add_same(&mut matrix, n, sector, 9, 4, q1);
        add_same(&mut matrix, n, sector, 9, 5, q2);
        add_same(&mut matrix, n, sector, 9, 6, q1);

        add_same(&mut matrix, n, sector, 10, 2, r);
        add_same(&mut matrix, n, sector, 10, 3, r);
        add_same(&mut matrix, n, sector, 10, 5, r);
        add_same(&mut matrix, n, sector, 10, 6, r);

        add_next(&mut matrix, n, sector, 11, 1, q1);
        add_same(&mut matrix, n, sector, 11, 2, q2);
        add_same(&mut matrix, n, sector, 11, 3, q2);
        add_next(&mut matrix, n, sector, 11, 4, q1);
        add_same(&mut matrix, n, sector, 11, 5, q1);
        add_same(&mut matrix, n, sector, 11, 6, q1);

        add_next(&mut matrix, n, sector, 12, 1, r);
        add_same(&mut matrix, n, sector, 12, 2, r);
        add_same(&mut matrix, n, sector, 12, 3, r);
        add_next(&mut matrix, n, sector, 12, 4, r);
    }

    Ok(matrix)
}

#[cfg(test)]
// Index form: `(dest_sector, dest_label)` also form the flattened subdivision
// matrix row index `dest_sector * 13 + dest_label`; keep that arithmetic visible.
#[allow(clippy::needless_range_loop)]
pub fn pr08_subdivide_characteristic_ring(
    ring: &Kp20CatmullClarkCharacteristicRing,
) -> Result<Vec<[Point3; 13]>, String> {
    let n = ring.sectors.len();
    let matrix = pr08_catmull_clark_subdivision_matrix(n)?;
    let mut out = vec![[Point3::new(0.0, 0.0, 0.0); 13]; n];
    for dest_sector in 0..n {
        for dest_label in 0..13 {
            let row = dest_sector * 13 + dest_label;
            let mut value = Point3::new(0.0, 0.0, 0.0);
            for src_sector in 0..n {
                for src_label in 0..13 {
                    let weight = matrix[row][src_sector * 13 + src_label];
                    let source = ring.sectors[src_sector][src_label];
                    value.x += weight * source.x;
                    value.y += weight * source.y;
                    value.z += weight * source.z;
                }
            }
            out[dest_sector][dest_label] = value;
        }
    }
    Ok(out)
}

#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kp20Pr08TensorBorderCandidate {
    /// The visually direct PR08 Fig. 6.3 strip
    /// `(q4,q7,q12,q13)`, `(q3,q6,q10,q11)`, `(q2,q5,q9,q8)`.
    VisibleStrip,
    /// Same Fig. 6.3 strip with depth order reversed.
    VisibleStripReversed,
    /// PR08 normalization text's middle-patch 3x3 tensor grid, extended by
    /// repeating the last visible column to make a diagnostic 4x3 row set.
    MiddlePatchFormulaRepeated,
    /// PR08 middle-patch formula grid with the fourth column completed from
    /// the visible right column of Fig. 6.3.
    MiddlePatchFormulaVisibleRight,
    /// Same as [`Self::MiddlePatchFormulaVisibleRight`] with depth order
    /// reversed.
    MiddlePatchFormulaVisibleRightReversed,
    /// PR08 middle-patch formula grid completed by a next-sector column.
    MiddlePatchFormulaNextColumn,
    /// Same next-sector completion with depth order reversed.
    MiddlePatchFormulaNextColumnReversed,
    /// PR08 middle-patch formula grid completed by the visible left column
    /// `(q2,q5,q8)` from Fig. 6.3.
    MiddlePatchFormulaLeftColumn,
    /// Same left-column completion with depth order reversed.
    MiddlePatchFormulaLeftColumnReversed,
    /// PR08 middle-patch 3x3 grid transposed across the adjacent `j+1`
    /// row visible in Fig. 6.3.
    MiddlePatchFormulaNextRowTransposed,
    /// PR08 middle-patch 3x3 grid transposed across the adjacent `j+1`
    /// outer row ending at the shared `q13` node.
    MiddlePatchFormulaNextOuterRowTransposed,
    /// Same transposed adjacent `j+1` outer row, but following PR08
    /// Fig. 6.3's visible outer continuation through `q_{j,11}` instead
    /// of repeating the shared `q_{j,13}` endpoint.
    MiddlePatchFormulaNextOuterRowTransposedToQ11,
    /// Same adjacent `j+1` outer-row transposed completion with depth order
    /// reversed.
    MiddlePatchFormulaNextOuterRowTransposedReversed,
    /// PR08 middle-patch 3x3 grid transposed with the first column treated
    /// as Stam Fig. 6 / Appendix-B subdivided `\bar C_n` labels
    /// `(4, 2n+7, 2n+16)`, not direct PR08 `q` aliases.
    MiddlePatchFormulaStamSubdividedFirstColumn,
    /// Stam Appendix-B `q^1` regular subpatch rows, evaluated from
    /// Appendix-A subdivided `\bar C_n` labels.
    StamAppendixASubpatchQ1InteriorRows,
    /// Stam Appendix-B `q^2` regular subpatch rows, evaluated from
    /// Appendix-A subdivided `\bar C_n` labels.
    StamAppendixASubpatchQ2InteriorRows,
    /// Stam Appendix-B `q^3` regular subpatch rows, evaluated from
    /// Appendix-A subdivided `\bar C_n` labels.
    StamAppendixASubpatchQ3InteriorRows,
    /// Stam Appendix-B `q^1` rows after feeding PR08 `q_{j,1..3}` into
    /// Appendix-A's extended subdivision vector.
    StamAppendixAPr08InnerRingQ1Rows,
    /// Stam Appendix-B `q^2` rows after feeding PR08 `q_{j,1..3}` into
    /// Appendix-A's extended subdivision vector.
    StamAppendixAPr08InnerRingQ2Rows,
    /// Stam Appendix-B `q^3` rows after feeding PR08 `q_{j,1..3}` into
    /// Appendix-A's extended subdivision vector.
    StamAppendixAPr08InnerRingQ3Rows,
    /// Stam Appendix-A/B rows after deriving the frequency-one top vector
    /// from PR08 `q13` (odd ring) and `q2` (even ring), with the complex
    /// similarity scale required by Stam's printed characteristic top vector.
    StamAppendixAPr08ScaledQ13Q2TopRows,
    /// Same adjacent `j+1` transposed completion with depth order reversed.
    MiddlePatchFormulaNextRowTransposedReversed,
    /// PR08 middle-patch 3x3 grid transposed across the adjacent `j-1`
    /// row visible in Fig. 6.3.
    MiddlePatchFormulaPreviousRowTransposed,
    /// Same adjacent `j-1` transposed completion with depth order reversed.
    MiddlePatchFormulaPreviousRowTransposedReversed,
    /// Stam Fig. 6 lower-right regular patch rows mapped to PR08 Fig. 6.3:
    /// `(q4,q7,q12,q13)`, `(q3,q6,q10,q11)`, `(q1,q2,q5,q8)`.
    StamLowerRightPatchRows,
    /// Same Stam lower-right row set with depth order reversed.
    StamLowerRightPatchRowsReversed,
    /// Adjacent three-row band of the Stam lower-right patch using the
    /// PR08 `j-1` row visible below Fig. 6.3.
    StamLowerRightStarSideRows,
    /// Same star-side lower-right row band with depth order reversed.
    StamLowerRightStarSideRowsReversed,
    /// Left-side band of the same full Stam lower-right patch grid.
    StamLowerRightLeftSideRows,
    /// Same left-side band with depth order reversed.
    StamLowerRightLeftSideRowsReversed,
    /// Right-side band of the same full Stam lower-right patch grid.
    StamLowerRightRightSideRows,
    /// Same right-side band with depth order reversed.
    StamLowerRightRightSideRowsReversed,
}

#[cfg(test)]
impl Kp20Pr08TensorBorderCandidate {
    pub const fn all() -> [Self; 32] {
        [
            Self::VisibleStrip,
            Self::VisibleStripReversed,
            Self::MiddlePatchFormulaRepeated,
            Self::MiddlePatchFormulaVisibleRight,
            Self::MiddlePatchFormulaVisibleRightReversed,
            Self::MiddlePatchFormulaNextColumn,
            Self::MiddlePatchFormulaNextColumnReversed,
            Self::MiddlePatchFormulaLeftColumn,
            Self::MiddlePatchFormulaLeftColumnReversed,
            Self::MiddlePatchFormulaNextRowTransposed,
            Self::MiddlePatchFormulaNextOuterRowTransposed,
            Self::MiddlePatchFormulaNextOuterRowTransposedToQ11,
            Self::MiddlePatchFormulaNextOuterRowTransposedReversed,
            Self::MiddlePatchFormulaStamSubdividedFirstColumn,
            Self::StamAppendixASubpatchQ1InteriorRows,
            Self::StamAppendixASubpatchQ2InteriorRows,
            Self::StamAppendixASubpatchQ3InteriorRows,
            Self::StamAppendixAPr08InnerRingQ1Rows,
            Self::StamAppendixAPr08InnerRingQ2Rows,
            Self::StamAppendixAPr08InnerRingQ3Rows,
            Self::StamAppendixAPr08ScaledQ13Q2TopRows,
            Self::MiddlePatchFormulaNextRowTransposedReversed,
            Self::MiddlePatchFormulaPreviousRowTransposed,
            Self::MiddlePatchFormulaPreviousRowTransposedReversed,
            Self::StamLowerRightPatchRows,
            Self::StamLowerRightPatchRowsReversed,
            Self::StamLowerRightStarSideRows,
            Self::StamLowerRightStarSideRowsReversed,
            Self::StamLowerRightLeftSideRows,
            Self::StamLowerRightLeftSideRowsReversed,
            Self::StamLowerRightRightSideRows,
            Self::StamLowerRightRightSideRowsReversed,
        ]
    }

    pub const fn name(self) -> &'static str {
        match self {
            Self::VisibleStrip => "pr08-visible-strip",
            Self::VisibleStripReversed => "pr08-visible-strip-reversed",
            Self::MiddlePatchFormulaRepeated => "pr08-middle-patch-formula-repeated",
            Self::MiddlePatchFormulaVisibleRight => "pr08-middle-patch-formula-visible-right",
            Self::MiddlePatchFormulaVisibleRightReversed => {
                "pr08-middle-patch-formula-visible-right-reversed"
            }
            Self::MiddlePatchFormulaNextColumn => "pr08-middle-patch-formula-next-column",
            Self::MiddlePatchFormulaNextColumnReversed => {
                "pr08-middle-patch-formula-next-column-reversed"
            }
            Self::MiddlePatchFormulaLeftColumn => "pr08-middle-patch-formula-left-column",
            Self::MiddlePatchFormulaLeftColumnReversed => {
                "pr08-middle-patch-formula-left-column-reversed"
            }
            Self::MiddlePatchFormulaNextRowTransposed => {
                "pr08-middle-patch-formula-next-row-transposed"
            }
            Self::MiddlePatchFormulaNextOuterRowTransposed => {
                "pr08-middle-patch-formula-next-outer-row-transposed"
            }
            Self::MiddlePatchFormulaNextOuterRowTransposedToQ11 => {
                "pr08-middle-patch-formula-next-outer-row-transposed-to-q11"
            }
            Self::MiddlePatchFormulaNextOuterRowTransposedReversed => {
                "pr08-middle-patch-formula-next-outer-row-transposed-reversed"
            }
            Self::MiddlePatchFormulaStamSubdividedFirstColumn => {
                "pr08-middle-patch-formula-stam-subdivided-first-column"
            }
            Self::StamAppendixASubpatchQ1InteriorRows => {
                "stam-appendix-a-subpatch-q1-interior-rows"
            }
            Self::StamAppendixASubpatchQ2InteriorRows => {
                "stam-appendix-a-subpatch-q2-interior-rows"
            }
            Self::StamAppendixASubpatchQ3InteriorRows => {
                "stam-appendix-a-subpatch-q3-interior-rows"
            }
            Self::StamAppendixAPr08InnerRingQ1Rows => "stam-appendix-a-pr08-inner-ring-q1-rows",
            Self::StamAppendixAPr08InnerRingQ2Rows => "stam-appendix-a-pr08-inner-ring-q2-rows",
            Self::StamAppendixAPr08InnerRingQ3Rows => "stam-appendix-a-pr08-inner-ring-q3-rows",
            Self::StamAppendixAPr08ScaledQ13Q2TopRows => {
                "stam-appendix-a-pr08-scaled-q13-q2-top-rows"
            }
            Self::MiddlePatchFormulaNextRowTransposedReversed => {
                "pr08-middle-patch-formula-next-row-transposed-reversed"
            }
            Self::MiddlePatchFormulaPreviousRowTransposed => {
                "pr08-middle-patch-formula-previous-row-transposed"
            }
            Self::MiddlePatchFormulaPreviousRowTransposedReversed => {
                "pr08-middle-patch-formula-previous-row-transposed-reversed"
            }
            Self::StamLowerRightPatchRows => "pr08-stam-lower-right-patch-rows",
            Self::StamLowerRightPatchRowsReversed => "pr08-stam-lower-right-patch-rows-reversed",
            Self::StamLowerRightStarSideRows => "pr08-stam-lower-right-star-side-rows",
            Self::StamLowerRightStarSideRowsReversed => {
                "pr08-stam-lower-right-star-side-rows-reversed"
            }
            Self::StamLowerRightLeftSideRows => "pr08-stam-lower-right-left-side-rows",
            Self::StamLowerRightLeftSideRowsReversed => {
                "pr08-stam-lower-right-left-side-rows-reversed"
            }
            Self::StamLowerRightRightSideRows => "pr08-stam-lower-right-right-side-rows",
            Self::StamLowerRightRightSideRowsReversed => {
                "pr08-stam-lower-right-right-side-rows-reversed"
            }
        }
    }
}

#[cfg(test)]
pub fn kp20_pr08_characteristic_tensor_border_candidate_rows(
    ring: &Kp20CatmullClarkCharacteristicRing,
    candidate: Kp20Pr08TensorBorderCandidate,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    if ring.sectors.len() < 3 {
        return Err(format!(
            "PR08 characteristic tensor-border candidate needs at least 3 sectors, got {}",
            ring.sectors.len()
        ));
    }
    let mut sectors = Vec::with_capacity(ring.sectors.len());
    for (sector_index, sector) in ring.sectors.iter().enumerate() {
        validate_pr08_sector(sector_index, sector)?;
        let previous = &ring.sectors[(sector_index + ring.sectors.len() - 1) % ring.sectors.len()];
        let next = &ring.sectors[(sector_index + 1) % ring.sectors.len()];
        let visible = [
            [sector[3], sector[6], sector[11], sector[12]],
            [sector[2], sector[5], sector[9], sector[10]],
            [sector[1], sector[4], sector[8], sector[7]],
        ];
        let formula_repeated = [
            [sector[2], sector[3], sector[12], sector[12]],
            [sector[5], sector[6], sector[11], sector[11]],
            [sector[8], sector[9], sector[12], sector[12]],
        ];
        let formula_visible_right = [
            [sector[2], sector[3], sector[12], sector[12]],
            [sector[5], sector[6], sector[11], sector[10]],
            [sector[8], sector[9], sector[12], sector[7]],
        ];
        let formula_next_column = [
            [sector[2], sector[3], sector[12], next[7]],
            [sector[5], sector[6], sector[11], next[8]],
            [sector[8], sector[9], sector[12], next[12]],
        ];
        let formula_left_column = [
            [sector[1], sector[2], sector[3], sector[12]],
            [sector[4], sector[5], sector[6], sector[11]],
            [sector[7], sector[8], sector[9], sector[12]],
        ];
        let formula_next_row_transposed = [
            [next[1], sector[2], sector[5], sector[8]],
            [next[4], sector[3], sector[6], sector[9]],
            [next[7], sector[12], sector[11], sector[12]],
        ];
        let formula_next_outer_row_transposed = [
            [next[3], sector[2], sector[5], sector[8]],
            [next[6], sector[3], sector[6], sector[9]],
            [next[12], sector[12], sector[11], sector[12]],
        ];
        let formula_next_outer_row_transposed_to_q11 = [
            [next[3], sector[2], sector[5], sector[8]],
            [next[6], sector[3], sector[6], sector[9]],
            [next[12], sector[12], sector[11], sector[10]],
        ];
        let stam_subdivided_first_column = stam_appendix_a_middle_patch_first_column(
            ring.lambda,
            ring.sectors.len(),
            sector_index,
        )?;
        let formula_stam_subdivided_first_column = [
            [
                stam_subdivided_first_column[0],
                sector[2],
                sector[5],
                sector[8],
            ],
            [
                stam_subdivided_first_column[1],
                sector[3],
                sector[6],
                sector[9],
            ],
            [
                stam_subdivided_first_column[2],
                sector[12],
                sector[11],
                sector[12],
            ],
        ];
        let stam_appendix_a_q1 = stam_appendix_a_regular_subpatch_interior_rows(
            ring.lambda,
            ring.sectors.len(),
            sector_index,
            0,
        )?;
        let stam_appendix_a_q2 = stam_appendix_a_regular_subpatch_interior_rows(
            ring.lambda,
            ring.sectors.len(),
            sector_index,
            1,
        )?;
        let stam_appendix_a_q3 = stam_appendix_a_regular_subpatch_interior_rows(
            ring.lambda,
            ring.sectors.len(),
            sector_index,
            2,
        )?;
        let stam_pr08_inner_q1 =
            stam_appendix_a_pr08_inner_ring_subpatch_rows(ring, sector_index, 0)?;
        let stam_pr08_inner_q2 =
            stam_appendix_a_pr08_inner_ring_subpatch_rows(ring, sector_index, 1)?;
        let stam_pr08_inner_q3 =
            stam_appendix_a_pr08_inner_ring_subpatch_rows(ring, sector_index, 2)?;
        let stam_pr08_scaled_q13_q2_top =
            stam_appendix_a_pr08_scaled_q13_q2_top_rows(ring, sector_index, 0, 1, true)?;
        let formula_previous_row_transposed = [
            [sector[2], sector[5], sector[8], previous[12]],
            [sector[3], sector[6], sector[9], previous[3]],
            [sector[12], sector[11], sector[12], previous[0]],
        ];
        let stam_lower_right = [
            [sector[3], sector[6], sector[11], sector[12]],
            [sector[2], sector[5], sector[9], sector[10]],
            [sector[0], sector[1], sector[4], sector[7]],
        ];
        let stam_lower_right_star_side = [
            [sector[2], sector[5], sector[9], sector[10]],
            [sector[0], sector[1], sector[4], sector[7]],
            [previous[0], previous[2], previous[3], previous[12]],
        ];
        let stam_lower_right_grid = [
            [sector[3], sector[6], sector[11], sector[12]],
            [sector[2], sector[5], sector[9], sector[10]],
            [sector[0], sector[1], sector[4], sector[7]],
            [previous[0], previous[2], previous[3], previous[12]],
        ];
        let stam_lower_right_left_side = [
            [
                stam_lower_right_grid[0][0],
                stam_lower_right_grid[1][0],
                stam_lower_right_grid[2][0],
                stam_lower_right_grid[3][0],
            ],
            [
                stam_lower_right_grid[0][1],
                stam_lower_right_grid[1][1],
                stam_lower_right_grid[2][1],
                stam_lower_right_grid[3][1],
            ],
            [
                stam_lower_right_grid[0][2],
                stam_lower_right_grid[1][2],
                stam_lower_right_grid[2][2],
                stam_lower_right_grid[3][2],
            ],
        ];
        let stam_lower_right_right_side = [
            [
                stam_lower_right_grid[0][3],
                stam_lower_right_grid[1][3],
                stam_lower_right_grid[2][3],
                stam_lower_right_grid[3][3],
            ],
            [
                stam_lower_right_grid[0][2],
                stam_lower_right_grid[1][2],
                stam_lower_right_grid[2][2],
                stam_lower_right_grid[3][2],
            ],
            [
                stam_lower_right_grid[0][1],
                stam_lower_right_grid[1][1],
                stam_lower_right_grid[2][1],
                stam_lower_right_grid[3][1],
            ],
        ];
        sectors.push(match candidate {
            Kp20Pr08TensorBorderCandidate::VisibleStrip => visible,
            Kp20Pr08TensorBorderCandidate::VisibleStripReversed => {
                [visible[2], visible[1], visible[0]]
            }
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaRepeated => formula_repeated,
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaVisibleRight => formula_visible_right,
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaVisibleRightReversed => [
                formula_visible_right[2],
                formula_visible_right[1],
                formula_visible_right[0],
            ],
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextColumn => formula_next_column,
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextColumnReversed => [
                formula_next_column[2],
                formula_next_column[1],
                formula_next_column[0],
            ],
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaLeftColumn => formula_left_column,
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaLeftColumnReversed => [
                formula_left_column[2],
                formula_left_column[1],
                formula_left_column[0],
            ],
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextRowTransposed => {
                formula_next_row_transposed
            }
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposed => {
                formula_next_outer_row_transposed
            }
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedToQ11 => {
                formula_next_outer_row_transposed_to_q11
            }
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextOuterRowTransposedReversed => [
                formula_next_outer_row_transposed[2],
                formula_next_outer_row_transposed[1],
                formula_next_outer_row_transposed[0],
            ],
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaStamSubdividedFirstColumn => {
                formula_stam_subdivided_first_column
            }
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ1InteriorRows => {
                stam_appendix_a_q1
            }
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ2InteriorRows => {
                stam_appendix_a_q2
            }
            Kp20Pr08TensorBorderCandidate::StamAppendixASubpatchQ3InteriorRows => {
                stam_appendix_a_q3
            }
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ1Rows => stam_pr08_inner_q1,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ2Rows => stam_pr08_inner_q2,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08InnerRingQ3Rows => stam_pr08_inner_q3,
            Kp20Pr08TensorBorderCandidate::StamAppendixAPr08ScaledQ13Q2TopRows => {
                stam_pr08_scaled_q13_q2_top
            }
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaNextRowTransposedReversed => [
                formula_next_row_transposed[2],
                formula_next_row_transposed[1],
                formula_next_row_transposed[0],
            ],
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaPreviousRowTransposed => {
                formula_previous_row_transposed
            }
            Kp20Pr08TensorBorderCandidate::MiddlePatchFormulaPreviousRowTransposedReversed => [
                formula_previous_row_transposed[2],
                formula_previous_row_transposed[1],
                formula_previous_row_transposed[0],
            ],
            Kp20Pr08TensorBorderCandidate::StamLowerRightPatchRows => stam_lower_right,
            Kp20Pr08TensorBorderCandidate::StamLowerRightPatchRowsReversed => [
                stam_lower_right[2],
                stam_lower_right[1],
                stam_lower_right[0],
            ],
            Kp20Pr08TensorBorderCandidate::StamLowerRightStarSideRows => stam_lower_right_star_side,
            Kp20Pr08TensorBorderCandidate::StamLowerRightStarSideRowsReversed => [
                stam_lower_right_star_side[2],
                stam_lower_right_star_side[1],
                stam_lower_right_star_side[0],
            ],
            Kp20Pr08TensorBorderCandidate::StamLowerRightLeftSideRows => stam_lower_right_left_side,
            Kp20Pr08TensorBorderCandidate::StamLowerRightLeftSideRowsReversed => [
                stam_lower_right_left_side[2],
                stam_lower_right_left_side[1],
                stam_lower_right_left_side[0],
            ],
            Kp20Pr08TensorBorderCandidate::StamLowerRightRightSideRows => {
                stam_lower_right_right_side
            }
            Kp20Pr08TensorBorderCandidate::StamLowerRightRightSideRowsReversed => [
                stam_lower_right_right_side[2],
                stam_lower_right_right_side[1],
                stam_lower_right_right_side[0],
            ],
        });
    }
    Ok(Kp20RegularTensorBorderCNet { sectors })
}

#[cfg(test)]
pub fn stam_appendix_a_middle_patch_first_column(
    lambda: f64,
    valence: usize,
    sector_index: usize,
) -> Result<[Point3; 3], String> {
    let labels = stam_appendix_a_subdivided_extended_characteristic_labels(lambda, valence)?;
    let rotate = |p: Point3| {
        let theta = core::f64::consts::TAU * sector_index as f64 / valence as f64;
        let (sin_theta, cos_theta) = theta.sin_cos();
        rotate_xy(p, sin_theta, cos_theta)
    };
    Ok([
        rotate(labels[3]),
        rotate(labels[2 * valence + 6]),
        rotate(labels[2 * valence + 15]),
    ])
}

fn stam_appendix_a_pr08_inner_ring_subpatch_rows(
    ring: &Kp20CatmullClarkCharacteristicRing,
    sector_index: usize,
    patch_index: usize,
) -> Result<[[Point3; 4]; 3], String> {
    let grid = stam_appendix_a_pr08_inner_ring_subpatch_full_grid(ring, sector_index, patch_index)?;
    Ok([grid[1], grid[2], grid[3]])
}

fn stam_appendix_a_pr08_inner_ring_subpatch_full_grid(
    ring: &Kp20CatmullClarkCharacteristicRing,
    sector_index: usize,
    patch_index: usize,
) -> Result<[[Point3; 4]; 4], String> {
    let labels = stam_appendix_a_pr08_inner_ring_subdivided_labels(ring)?;
    stam_appendix_b_subpatch_full_grid_from_labels(
        &labels,
        ring.sectors.len(),
        sector_index,
        patch_index,
    )
}

#[cfg(test)]
pub fn kp20_pr08_stam_appendix_b_inner_ring_full_grid_for_test(
    ring: &Kp20CatmullClarkCharacteristicRing,
    sector_index: usize,
    patch_index: usize,
) -> Result<[[Point3; 4]; 4], String> {
    stam_appendix_a_pr08_inner_ring_subpatch_full_grid(ring, sector_index, patch_index)
}

#[cfg(test)]
pub fn kp20_pr08_stam_appendix_b_projected_rows_from_top_vectors_for_test(
    ring: &Kp20CatmullClarkCharacteristicRing,
    u0_x: &[f64],
    u0_y: &[f64],
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    let valence = ring.sectors.len();
    let labels = stam_appendix_a_subdivide_extended_labels_from_top(ring.lambda, u0_x, u0_y)?;
    let mut sectors = Vec::with_capacity(valence);
    for sector_index in 0..valence {
        let grid = stam_appendix_b_subpatch_full_grid_from_labels(
            &labels,
            valence,
            sector_index,
            patch_index,
        )?;
        let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for (out_row, row) in rows.iter_mut().enumerate() {
            let source_row = if reverse_rows { 2 - out_row } else { out_row };
            *row = grid[row_start + source_row];
        }
        sectors.push(rows);
    }
    Ok(Kp20RegularTensorBorderCNet { sectors })
}

#[cfg(test)]
fn stam_appendix_a_pr08_scaled_q13_q2_top_rows(
    ring: &Kp20CatmullClarkCharacteristicRing,
    sector_index: usize,
    patch_index: usize,
    row_start: usize,
    reverse_rows: bool,
) -> Result<[[Point3; 4]; 3], String> {
    let grid = stam_appendix_a_pr08_scaled_q13_q2_top_full_grid(ring, sector_index, patch_index)?;
    let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
    for (out_row, row) in rows.iter_mut().enumerate() {
        let source_row = if reverse_rows { 2 - out_row } else { out_row };
        *row = grid[row_start + source_row];
    }
    Ok(rows)
}

#[cfg(test)]
fn stam_appendix_a_pr08_scaled_q13_q2_top_full_grid(
    ring: &Kp20CatmullClarkCharacteristicRing,
    sector_index: usize,
    patch_index: usize,
) -> Result<[[Point3; 4]; 4], String> {
    let labels = stam_appendix_a_pr08_scaled_q13_q2_top_labels(ring)?;
    stam_appendix_b_subpatch_full_grid_from_labels(
        &labels,
        ring.sectors.len(),
        sector_index,
        patch_index,
    )
}

#[cfg(test)]
fn stam_appendix_a_pr08_scaled_q13_q2_top_labels(
    ring: &Kp20CatmullClarkCharacteristicRing,
) -> Result<Vec<Point3>, String> {
    let valence = ring.sectors.len();
    if valence < 4 {
        return Err(format!(
            "Stam Appendix-A PR08 q13/q2 source needs valence >= 4, got {valence}"
        ));
    }
    let top_count = 2 * valence + 1;
    let stam_top = stam_appendix_a_frequency_one_top_vector(ring);
    let q13 = pr08_label_ring_as_complex(ring, 12)?;
    let q2 = pr08_label_ring_as_complex(ring, 1)?;
    let odd_scale = complex_similarity_ring_scale(&q13, &stam_top.odd)
        .ok_or_else(|| "degenerate PR08 q13 ring for Stam top-vector scaling".to_string())?;
    let even_scale = complex_similarity_ring_scale(&q2, &stam_top.even)
        .ok_or_else(|| "degenerate PR08 q2 ring for Stam top-vector scaling".to_string())?;

    let mut u0_x = vec![0.0_f64; top_count];
    let mut u0_y = vec![0.0_f64; top_count];
    u0_x[0] = ring.sectors[0][0].x;
    u0_y[0] = ring.sectors[0][0].y;
    for sector in 0..valence {
        let odd = Complex::new(ring.sectors[sector][12].x, ring.sectors[sector][12].y) * odd_scale;
        let even = Complex::new(ring.sectors[sector][1].x, ring.sectors[sector][1].y) * even_scale;
        u0_x[1 + 2 * sector] = odd.re;
        u0_y[1 + 2 * sector] = odd.im;
        u0_x[2 + 2 * sector] = even.re;
        u0_y[2 + 2 * sector] = even.im;
    }
    stam_appendix_a_subdivide_extended_labels_from_top(ring.lambda, &u0_x, &u0_y)
}

fn stam_appendix_a_regular_subpatch_interior_rows(
    lambda: f64,
    valence: usize,
    sector_index: usize,
    patch_index: usize,
) -> Result<[[Point3; 4]; 3], String> {
    let grid =
        stam_appendix_a_regular_subpatch_full_grid(lambda, valence, sector_index, patch_index)?;
    Ok([grid[1], grid[2], grid[3]])
}

fn stam_appendix_a_regular_subpatch_full_grid(
    lambda: f64,
    valence: usize,
    sector_index: usize,
    patch_index: usize,
) -> Result<[[Point3; 4]; 4], String> {
    let labels = stam_appendix_a_subdivided_extended_characteristic_labels(lambda, valence)?;
    stam_appendix_b_subpatch_full_grid_from_labels(&labels, valence, sector_index, patch_index)
}

fn stam_appendix_b_subpatch_full_grid_from_labels(
    labels: &[Point3],
    valence: usize,
    sector_index: usize,
    patch_index: usize,
) -> Result<[[Point3; 4]; 4], String> {
    let label = |one_based: usize| -> Result<Point3, String> {
        labels.get(one_based - 1).copied().ok_or_else(|| {
            format!(
                "Stam Appendix-B label {one_based} is outside subdivided label count {}",
                labels.len()
            )
        })
    };
    let n2 = 2 * valence;
    let grid_labels = match patch_index {
        0 => [
            [8, 7, n2 + 5, n2 + 13],
            [1, 6, n2 + 4, n2 + 12],
            [4, 5, n2 + 3, n2 + 11],
            [n2 + 7, n2 + 6, n2 + 2, n2 + 10],
        ],
        1 => [
            [1, 6, n2 + 4, n2 + 12],
            [4, 5, n2 + 3, n2 + 11],
            [n2 + 7, n2 + 6, n2 + 2, n2 + 10],
            [n2 + 16, n2 + 15, n2 + 14, n2 + 9],
        ],
        2 => [
            [2, 1, 6, n2 + 4],
            [3, 4, 5, n2 + 3],
            [n2 + 8, n2 + 7, n2 + 6, n2 + 2],
            [n2 + 17, n2 + 16, n2 + 15, n2 + 14],
        ],
        _ => {
            return Err(format!(
                "Stam Appendix-B subpatch index must be 0, 1, or 2, got {patch_index}"
            ));
        }
    };
    let theta = core::f64::consts::TAU * sector_index as f64 / valence as f64;
    let (sin_theta, cos_theta) = theta.psin_cos();
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for row in 0..4 {
        for col in 0..4 {
            grid[row][col] = rotate_xy(label(grid_labels[row][col])?, sin_theta, cos_theta);
        }
    }
    Ok(grid)
}

fn rotate_xy(p: Point3, sin_theta: f64, cos_theta: f64) -> Point3 {
    Point3::new(
        cos_theta * p.x - sin_theta * p.y,
        sin_theta * p.x + cos_theta * p.y,
        p.z,
    )
}

fn endpoint_normalized_regular_tensor_border_rows(
    rows: &Kp20RegularTensorBorderCNet,
    valence: usize,
    sector: usize,
    swap_targets: bool,
    conjugate_source: bool,
) -> Result<Kp20RegularTensorBorderCNet, String> {
    if valence < 3 {
        return Err(format!(
            "PR08/Stam endpoint normalization needs valence >= 3; got {valence}"
        ));
    }
    let sector_rows = rows
        .sectors
        .get(sector)
        .ok_or_else(|| format!("endpoint normalization sector {sector} is missing"))?;
    let alpha = core::f64::consts::TAU / valence as f64;
    let spoke0 = Complex::from_angle(alpha * sector as f64);
    let spoke1 = Complex::from_angle(alpha * (sector + 1) as f64);
    let target0 = if swap_targets { spoke1 } else { spoke0 };
    let target1 = if swap_targets { spoke0 } else { spoke1 };
    let source0 = maybe_conjugate_complex(point_as_complex(sector_rows[0][0]), conjugate_source);
    let source1 = maybe_conjugate_complex(point_as_complex(sector_rows[0][3]), conjugate_source);
    let scale = complex_similarity_scale(&[(source0, target0), (source1, target1)])
        .ok_or_else(|| "degenerate PR08/Stam endpoint normalization source".to_string())?;
    Ok(complex_scaled_regular_tensor_border_rows(
        rows,
        scale,
        conjugate_source,
    ))
}

fn point_as_complex(point: Point3) -> Complex {
    Complex::new(point.x, point.y)
}

fn maybe_conjugate_complex(value: Complex, conjugate: bool) -> Complex {
    if conjugate { value.conj() } else { value }
}

fn complex_similarity_scale(pairs: &[(Complex, Complex)]) -> Option<Complex> {
    let mut numerator = Complex::zero();
    let mut denominator = 0.0_f64;
    for &(source, target) in pairs {
        numerator = numerator + target * source.conj();
        denominator += source.norm2();
    }
    if denominator <= 1.0e-20 {
        None
    } else {
        Some(numerator / denominator)
    }
}

fn complex_scaled_regular_tensor_border_rows(
    rows: &Kp20RegularTensorBorderCNet,
    scale: Complex,
    conjugate_source: bool,
) -> Kp20RegularTensorBorderCNet {
    let sectors = rows
        .sectors
        .iter()
        .map(|sector_rows| {
            let mut out = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
            for row in 0..3 {
                for col in 0..4 {
                    let value = maybe_conjugate_complex(
                        point_as_complex(sector_rows[row][col]),
                        conjugate_source,
                    ) * scale;
                    out[row][col] = Point3::new(value.re, value.im, sector_rows[row][col].z);
                }
            }
            out
        })
        .collect();
    Kp20RegularTensorBorderCNet { sectors }
}

fn stam_appendix_a_subdivided_extended_characteristic_labels(
    lambda: f64,
    valence: usize,
) -> Result<Vec<Point3>, String> {
    if valence < 4 {
        return Err(format!(
            "Stam Appendix-A diagnostic currently needs valence >= 4, got {valence}"
        ));
    }
    let top_count = 2 * valence + 1;
    let mut u0_x = vec![0.0_f64; top_count];
    let mut u0_y = vec![0.0_f64; top_count];
    let gain = 4.0 * lambda - 1.0;
    u0_x[0] = 0.0;
    u0_y[0] = 0.0;
    for sector in 0..valence {
        let theta = core::f64::consts::TAU * sector as f64 / valence as f64;
        let next_theta = core::f64::consts::TAU * (sector + 1) as f64 / valence as f64;
        let (sin_theta, cos_theta) = theta.psin_cos();
        let (sin_next, cos_next) = next_theta.psin_cos();
        u0_x[1 + 2 * sector] = gain * cos_theta;
        u0_y[1 + 2 * sector] = gain * sin_theta;
        u0_x[2 + 2 * sector] = cos_theta + cos_next;
        u0_y[2 + 2 * sector] = sin_theta + sin_next;
    }

    let u1_x = stam_appendix_a_outer_controls(lambda, &u0_x)?;
    let u1_y = stam_appendix_a_outer_controls(lambda, &u0_y)?;
    let mut labels = Vec::with_capacity(2 * valence + 17);
    for idx in 0..top_count {
        labels.push(Point3::new(lambda * u0_x[idx], lambda * u0_y[idx], 0.0));
    }
    for idx in 0..7 {
        labels.push(Point3::new(lambda * u1_x[idx], lambda * u1_y[idx], 0.0));
    }
    let extra_x = stam_appendix_a_extra_controls(&u0_x, &u1_x);
    let extra_y = stam_appendix_a_extra_controls(&u0_y, &u1_y);
    for idx in 0..9 {
        labels.push(Point3::new(extra_x[idx], extra_y[idx], 0.0));
    }
    Ok(labels)
}

fn stam_appendix_a_pr08_inner_ring_subdivided_labels(
    ring: &Kp20CatmullClarkCharacteristicRing,
) -> Result<Vec<Point3>, String> {
    let valence = ring.sectors.len();
    if valence < 4 {
        return Err(format!(
            "Stam Appendix-A PR08 inner-ring diagnostic currently needs valence >= 4, got {valence}"
        ));
    }
    let top_count = 2 * valence + 1;
    let mut u0_x = vec![0.0_f64; top_count];
    let mut u0_y = vec![0.0_f64; top_count];
    u0_x[0] = ring.sectors[0][0].x;
    u0_y[0] = ring.sectors[0][0].y;
    for sector in 0..valence {
        u0_x[1 + 2 * sector] = ring.sectors[sector][1].x;
        u0_y[1 + 2 * sector] = ring.sectors[sector][1].y;
        u0_x[2 + 2 * sector] = ring.sectors[sector][2].x;
        u0_y[2 + 2 * sector] = ring.sectors[sector][2].y;
    }
    stam_appendix_a_subdivide_extended_labels_from_top(ring.lambda, &u0_x, &u0_y)
}

#[cfg(test)]
struct StamAppendixATopVectorRings {
    odd: Vec<Complex>,
    even: Vec<Complex>,
}

#[cfg(test)]
fn stam_appendix_a_frequency_one_top_vector(
    ring: &Kp20CatmullClarkCharacteristicRing,
) -> StamAppendixATopVectorRings {
    let valence = ring.sectors.len();
    let gain = 4.0 * ring.lambda - 1.0;
    let mut odd = Vec::with_capacity(valence);
    let mut even = Vec::with_capacity(valence);
    for sector in 0..valence {
        let theta = core::f64::consts::TAU * sector as f64 / valence as f64;
        let next_theta = core::f64::consts::TAU * (sector + 1) as f64 / valence as f64;
        odd.push(Complex::new(gain * theta.cos(), gain * theta.sin()));
        even.push(Complex::new(
            theta.cos() + next_theta.cos(),
            theta.sin() + next_theta.sin(),
        ));
    }
    StamAppendixATopVectorRings { odd, even }
}

#[cfg(test)]
fn pr08_label_ring_as_complex(
    ring: &Kp20CatmullClarkCharacteristicRing,
    label_index: usize,
) -> Result<Vec<Complex>, String> {
    if label_index >= 13 {
        return Err(format!(
            "PR08 label index must be in 0..13, got {label_index}"
        ));
    }
    Ok(ring
        .sectors
        .iter()
        .map(|sector| Complex::new(sector[label_index].x, sector[label_index].y))
        .collect())
}

#[cfg(test)]
fn complex_similarity_ring_scale(source: &[Complex], target: &[Complex]) -> Option<Complex> {
    if source.len() != target.len() {
        return None;
    }
    let mut numerator = Complex::zero();
    let mut denominator = 0.0_f64;
    for (s, t) in source.iter().zip(target) {
        numerator = numerator + *t * s.conj();
        denominator += s.norm2();
    }
    if denominator < 1.0e-20 {
        return None;
    }
    Some(numerator / denominator)
}

fn stam_appendix_a_subdivide_extended_labels_from_top(
    lambda: f64,
    u0_x: &[f64],
    u0_y: &[f64],
) -> Result<Vec<Point3>, String> {
    let top_count = u0_x.len();
    if u0_y.len() != top_count {
        return Err(format!(
            "Stam Appendix-A x/y top rows must have the same length, got {} and {}",
            top_count,
            u0_y.len()
        ));
    }
    let u1_x = stam_appendix_a_outer_controls(lambda, u0_x)?;
    let u1_y = stam_appendix_a_outer_controls(lambda, u0_y)?;
    let mut labels = Vec::with_capacity(top_count + 16);
    for idx in 0..top_count {
        labels.push(Point3::new(lambda * u0_x[idx], lambda * u0_y[idx], 0.0));
    }
    for idx in 0..7 {
        labels.push(Point3::new(lambda * u1_x[idx], lambda * u1_y[idx], 0.0));
    }
    let extra_x = stam_appendix_a_extra_controls(u0_x, &u1_x);
    let extra_y = stam_appendix_a_extra_controls(u0_y, &u1_y);
    for idx in 0..9 {
        labels.push(Point3::new(extra_x[idx], extra_y[idx], 0.0));
    }
    Ok(labels)
}

fn stam_appendix_a_outer_controls(lambda: f64, u0: &[f64]) -> Result<Vec<f64>, String> {
    let s11 = stam_appendix_a_s11(u0.len());
    let s12 = stam_appendix_a_s12();
    let mut matrix = s12;
    for (idx, row) in matrix.iter_mut().enumerate() {
        row[idx] -= lambda;
    }
    let mut rhs = vec![0.0_f64; 7];
    for row in 0..7 {
        rhs[row] = -s11[row]
            .iter()
            .zip(u0)
            .map(|(weight, value)| weight * value)
            .sum::<f64>();
    }
    solve_real_dense(matrix.iter().map(|row| row.to_vec()).collect(), rhs)
        .ok_or_else(|| "singular Stam Appendix-A outer-control solve".to_string())
}

fn stam_appendix_a_extra_controls(u0: &[f64], u1: &[f64]) -> [f64; 9] {
    let s21 = stam_appendix_a_s21(u0.len());
    let s22 = stam_appendix_a_s22();
    let mut out = [0.0_f64; 9];
    for row in 0..9 {
        out[row] += s21[row]
            .iter()
            .zip(u0)
            .map(|(weight, value)| weight * value)
            .sum::<f64>();
        out[row] += s22[row]
            .iter()
            .zip(u1)
            .map(|(weight, value)| weight * value)
            .sum::<f64>();
    }
    out
}

fn stam_appendix_a_s11(top_count: usize) -> Vec<Vec<f64>> {
    let a = 9.0 / 16.0;
    let b = 3.0 / 32.0;
    let c = 1.0 / 64.0;
    let d = 3.0 / 8.0;
    let e = 1.0 / 16.0;
    let mut rows = vec![vec![0.0_f64; top_count]; 7];
    let mut set = |row: usize, col: usize, value: f64| {
        if col < top_count {
            rows[row][col] = value;
        }
    };
    set(0, 0, c);
    set(0, 3, b);
    set(0, 4, a);
    set(0, 5, b);
    set(1, 0, e);
    set(1, 3, e);
    set(1, 4, d);
    set(1, 5, d);
    set(2, 0, b);
    set(2, 3, c);
    set(2, 4, b);
    set(2, 5, a);
    set(2, 6, b);
    set(2, 7, c);
    set(3, 0, e);
    set(3, 5, d);
    set(3, 6, d);
    set(3, 7, e);
    set(4, 0, e);
    set(4, 3, d);
    set(4, 4, d);
    set(4, 5, e);
    set(5, 0, b);
    set(5, 1, c);
    set(5, 2, b);
    set(5, 3, a);
    set(5, 4, b);
    set(5, 5, c);
    set(6, 0, e);
    set(6, 1, e);
    set(6, 2, d);
    set(6, 3, d);
    rows
}

fn stam_appendix_a_s12() -> [[f64; 7]; 7] {
    let b = 3.0 / 32.0;
    let c = 1.0 / 64.0;
    let e = 1.0 / 16.0;
    [
        [c, b, c, 0.0, b, c, 0.0],
        [0.0, e, e, 0.0, 0.0, 0.0, 0.0],
        [0.0, c, b, c, 0.0, 0.0, 0.0],
        [0.0, 0.0, e, e, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, e, e, 0.0],
        [0.0, 0.0, 0.0, 0.0, c, b, c],
        [0.0, 0.0, 0.0, 0.0, 0.0, e, e],
    ]
}

fn stam_appendix_a_s21(top_count: usize) -> Vec<Vec<f64>> {
    let d = 3.0 / 8.0;
    let e = 1.0 / 16.0;
    let f = 1.0 / 4.0;
    let mut rows = vec![vec![0.0_f64; top_count]; 9];
    let mut set = |row: usize, col: usize, value: f64| {
        if col < top_count {
            rows[row][col] = value;
        }
    };
    set(0, 4, f);
    set(1, 4, d);
    set(1, 5, e);
    set(2, 4, f);
    set(2, 5, f);
    set(3, 4, e);
    set(3, 5, d);
    set(3, 6, e);
    set(4, 5, f);
    set(4, 6, f);
    set(5, 3, e);
    set(5, 4, d);
    set(6, 3, f);
    set(6, 4, f);
    set(7, 2, e);
    set(7, 3, d);
    set(7, 4, e);
    set(8, 2, f);
    set(8, 3, f);
    rows
}

fn stam_appendix_a_s22() -> [[f64; 7]; 9] {
    let d = 3.0 / 8.0;
    let e = 1.0 / 16.0;
    let f = 1.0 / 4.0;
    [
        [f, f, 0.0, 0.0, f, 0.0, 0.0],
        [e, d, e, 0.0, e, 0.0, 0.0],
        [0.0, f, f, 0.0, 0.0, 0.0, 0.0],
        [0.0, e, d, e, 0.0, 0.0, 0.0],
        [0.0, 0.0, f, f, 0.0, 0.0, 0.0],
        [e, e, 0.0, 0.0, d, e, 0.0],
        [0.0, 0.0, 0.0, 0.0, f, f, 0.0],
        [0.0, 0.0, 0.0, 0.0, e, d, e],
        [0.0, 0.0, 0.0, 0.0, 0.0, f, f],
    ]
}

fn solve_real_dense(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Option<Vec<f64>> {
    let n = rhs.len();
    for pivot in 0..n {
        let mut pivot_row = pivot;
        let mut pivot_abs = matrix[pivot][pivot].abs();
        for (row, values) in matrix.iter().enumerate().skip(pivot + 1) {
            let abs = values[pivot].abs();
            if abs > pivot_abs {
                pivot_abs = abs;
                pivot_row = row;
            }
        }
        if pivot_abs < 1.0e-14 {
            return None;
        }
        matrix.swap(pivot, pivot_row);
        rhs.swap(pivot, pivot_row);
        let diagonal = matrix[pivot][pivot];
        for value in &mut matrix[pivot][pivot..n] {
            *value /= diagonal;
        }
        rhs[pivot] /= diagonal;
        let pivot_values = matrix[pivot].clone();
        let pivot_rhs = rhs[pivot];
        for row in 0..n {
            if row == pivot {
                continue;
            }
            let factor = matrix[row][pivot];
            if factor == 0.0 {
                continue;
            }
            for (col, pivot_value) in pivot_values.iter().enumerate().skip(pivot) {
                matrix[row][col] -= factor * pivot_value;
            }
            rhs[row] -= factor * pivot_rhs;
        }
    }
    Some(rhs)
}

#[cfg(test)]
fn validate_pr08_sector(sector_index: usize, sector: &[Point3; 13]) -> Result<(), String> {
    for (coeff, p) in sector.iter().enumerate() {
        if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
            return Err(format!(
                "PR08 characteristic ring sector {sector_index} coeff q{} is non-finite",
                coeff + 1
            ));
        }
    }
    Ok(())
}

fn pr08_catmull_clark_fourier_one_eigenvector(n: usize) -> Result<[Complex; 13], String> {
    if n < 3 {
        return Err(format!(
            "PR08 Catmull-Clark characteristic ring needs n >= 3; got {n}"
        ));
    }
    let scalars = Kp20CharacteristicScalars::for_valence(n)?;
    let c = scalars.c;
    let lambda = scalars.lambda;
    let w = Complex::from_angle(core::f64::consts::TAU / n as f64);
    let zero = Complex::zero();
    let one = Complex::new(1.0, 0.0);

    let v0 = [
        zero,
        one + w,
        Complex::new(16.0 * lambda - 2.0 * c - 6.0, 0.0),
    ];
    let (a10, a11, a20, a21) = pr08_catmull_clark_lower_blocks(w);
    let mut rhs = [Complex::zero(); 4];
    for row in 0..4 {
        for col in 0..3 {
            rhs[row] = rhs[row] - a10[row][col] * v0[col];
        }
    }
    let mut system = a11;
    for (i, row) in system.iter_mut().enumerate() {
        row[i] = row[i] - Complex::new(lambda, 0.0);
    }
    let v1 = solve_complex_4(system, rhs)?;

    let mut v2 = [Complex::zero(); 6];
    for row in 0..6 {
        let mut value = Complex::zero();
        for col in 0..3 {
            value = value + a20[row][col] * v0[col];
        }
        for col in 0..4 {
            value = value + a21[row][col] * v1[col];
        }
        v2[row] = value / lambda;
    }

    let mut out = [Complex::zero(); 13];
    out[0..3].copy_from_slice(&v0);
    out[3..7].copy_from_slice(&v1);
    out[7..13].copy_from_slice(&v2);
    Ok(out)
}

#[cfg(test)]
type Matrix3 = [[Complex; 3]; 3];
type Matrix4 = [[Complex; 4]; 4];
type Matrix4x3 = [[Complex; 3]; 4];
type Matrix6x3 = [[Complex; 3]; 6];
type Matrix6x4 = [[Complex; 4]; 6];

#[cfg(test)]
fn pr08_catmull_clark_a00(n: usize) -> Matrix3 {
    pr08_catmull_clark_a00_for_mode(n, 1)
}

#[cfg(test)]
fn pr08_catmull_clark_a00_for_mode(n: usize, mode: usize) -> Matrix3 {
    let w = Complex::from_angle(core::f64::consts::TAU * mode as f64 / n as f64);
    let c = (core::f64::consts::TAU * mode as f64 / n as f64).cos();
    let delta = if mode == 0 { 1.0 } else { 0.0 };
    let alpha = 1.0 - 7.0 / (4.0 * n as f64);
    let beta = 3.0 / (2.0 * n as f64);
    let gamma = 1.0 / (4.0 * n as f64);
    [
        [
            Complex::new(alpha * delta, 0.0),
            Complex::new(beta * delta, 0.0),
            Complex::new(gamma * delta, 0.0),
        ],
        [
            Complex::new(3.0 / 8.0 * delta, 0.0),
            Complex::new(1.0 / 8.0 * c + 3.0 / 8.0, 0.0),
            (Complex::new(1.0, 0.0) + w) * (1.0 / 16.0),
        ],
        [
            Complex::new(1.0 / 4.0 * delta, 0.0),
            (Complex::new(1.0, 0.0) + w.conj()) * (1.0 / 4.0),
            Complex::new(1.0 / 4.0, 0.0),
        ],
    ]
}

#[cfg(test)]
fn pr08_catmull_clark_a_hat_block(n: usize, mode: usize) -> [[Complex; 13]; 13] {
    let w = Complex::from_angle(core::f64::consts::TAU * mode as f64 / n as f64);
    let a00 = pr08_catmull_clark_a00_for_mode(n, mode);
    let (a10, a11, a20, a21) = pr08_catmull_clark_lower_blocks_for_mode(w, mode);
    let mut block = [[Complex::zero(); 13]; 13];
    for row in 0..3 {
        for col in 0..3 {
            block[row][col] = a00[row][col];
        }
    }
    for row in 0..4 {
        for col in 0..3 {
            block[3 + row][col] = a10[row][col];
        }
        for col in 0..4 {
            block[3 + row][3 + col] = a11[row][col];
        }
    }
    for row in 0..6 {
        for col in 0..3 {
            block[7 + row][col] = a20[row][col];
        }
        for col in 0..4 {
            block[7 + row][3 + col] = a21[row][col];
        }
    }
    block
}

#[cfg(test)]
fn pr08_catmull_clark_lower_blocks_for_mode(
    w: Complex,
    mode: usize,
) -> (Matrix4x3, Matrix4, Matrix6x3, Matrix6x4) {
    let (mut a10, a11, a20, a21) = pr08_catmull_clark_lower_blocks(w);
    if mode == 0 {
        a10[0][0] = (3.0 / 8.0).into();
        a10[1][0] = (3.0 / 32.0).into();
        a10[2][0] = (1.0 / 16.0).into();
        a10[3][0] = (1.0 / 64.0).into();
    }
    (a10, a11, a20, a21)
}

fn pr08_catmull_clark_lower_blocks(w: Complex) -> (Matrix4x3, Matrix4, Matrix6x3, Matrix6x4) {
    let z = Complex::zero();
    let wc = w.conj();
    let c = w.re;
    let p1 = 1.0 / 64.0;
    let p2 = 3.0 / 32.0;
    let p3 = 9.0 / 16.0;
    let q1 = 1.0 / 16.0;
    let q2 = 3.0 / 8.0;
    let r = 1.0 / 4.0;
    let one = Complex::new(1.0, 0.0);

    let a10 = [
        [z, Complex::new(q1, 0.0) + w * q2, q1.into()],
        [z, (2.0 * p1 * c + p3).into(), (one + wc) * p2],
        [z, w * q1 + q2.into(), q2.into()],
        [z, (one + w) * p2, p3.into()],
    ];
    let a11 = [
        [q1.into(), w * q1, z, z],
        [wc * p1, p2.into(), p1.into(), z],
        [z, q1.into(), q1.into(), z],
        [p2.into(), (one + w) * p1, p2.into(), p1.into()],
    ];
    let a20 = [
        [z, q2.into(), (one + wc) * q1],
        [z, r.into(), r.into()],
        [z, q1.into(), q2.into()],
        [z, z, r.into()],
        [z, w * q1, q2.into()],
        [z, w * r, r.into()],
    ];
    let a21 = [
        [wc * q1, q2.into(), q1.into(), z],
        [z, r.into(), r.into(), z],
        [q1.into(), q1.into(), q2.into(), q1.into()],
        [r.into(), z, r.into(), r.into()],
        [q2.into(), w * q1, q1.into(), q1.into()],
        [r.into(), w * r, z, z],
    ];
    (a10, a11, a20, a21)
}

fn solve_complex_4(
    mut matrix: [[Complex; 4]; 4],
    mut rhs: [Complex; 4],
) -> Result<[Complex; 4], String> {
    for pivot in 0..4 {
        let mut pivot_row = pivot;
        let mut pivot_norm = matrix[pivot][pivot].norm2();
        for (row, values) in matrix.iter().enumerate().skip(pivot + 1) {
            let norm = values[pivot].norm2();
            if norm > pivot_norm {
                pivot_row = row;
                pivot_norm = norm;
            }
        }
        if pivot_norm < 1e-28 {
            return Err("singular PR08 Catmull-Clark eigenvector solve".to_string());
        }
        if pivot_row != pivot {
            matrix.swap(pivot, pivot_row);
            rhs.swap(pivot, pivot_row);
        }
        let diagonal = matrix[pivot][pivot];
        for value in &mut matrix[pivot][pivot..4] {
            *value = *value / diagonal;
        }
        rhs[pivot] = rhs[pivot] / diagonal;
        let pivot_row_values = matrix[pivot];
        let pivot_rhs = rhs[pivot];
        for row in 0..4 {
            if row == pivot {
                continue;
            }
            let factor = matrix[row][pivot];
            if factor.norm2() == 0.0 {
                continue;
            }
            for (col, pivot_value) in pivot_row_values.iter().enumerate().skip(pivot) {
                matrix[row][col] = matrix[row][col] - factor * *pivot_value;
            }
            rhs[row] = rhs[row] - factor * pivot_rhs;
        }
    }
    Ok(rhs)
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    const fn zero() -> Self {
        Self::new(0.0, 0.0)
    }

    fn from_angle(theta: f64) -> Self {
        let (im, re) = theta.psin_cos();
        Self::new(re, im)
    }

    const fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    fn norm2(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
}

impl From<f64> for Complex {
    fn from(value: f64) -> Self {
        Self::new(value, 0.0)
    }
}

impl core::ops::Add for Complex {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl core::ops::Sub for Complex {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl core::ops::Mul for Complex {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl core::ops::Mul<f64> for Complex {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self::Output {
        Self::new(self.re * rhs, self.im * rhs)
    }
}

impl core::ops::Div<f64> for Complex {
    type Output = Self;
    fn div(self, rhs: f64) -> Self::Output {
        Self::new(self.re / rhs, self.im / rhs)
    }
}

impl core::ops::Div for Complex {
    type Output = Self;
    fn div(self, rhs: Self) -> Self::Output {
        let denom = rhs.norm2();
        Self::new(
            (self.re * rhs.re + self.im * rhs.im) / denom,
            (self.im * rhs.re - self.re * rhs.im) / denom,
        )
    }
}

#[cfg(test)]
#[path = "kp20_bi6_cc_ring_tests.rs"]
mod tests;
