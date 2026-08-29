//! The chart-orientation repair's polish problem (PK-249) — split
//! from `three_quad_chart_orientation` for the 500-line structural
//! cap. Owns the TFI anchor, the localized pull-weight field, the
//! polish L-M evaluator, and the DoF packing.

use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};

use super::super::super::kp_binonic_optimizer::{
    analytic_ring_jacobian, outer_fillet_jacobian_sampled, ring_idx,
};
use super::super::super::kp_lm_solver::LmEval;
use super::three_quad_chart_orientation_metrics::{
    any_self_intersects, g2_residual_at, outer_fillet_sin,
};
use super::{binonic_reversal_band_count, interior_opposed_count};

const TETHER_WEIGHT: f64 = 3e-2;

/// Weight of the folded curvature-agreement rows in the polish —
/// the G2 gate is 5e-2 like G1, so the two fields carry equal
/// weight per sample.
const G2_WEIGHT: f64 = 1.0;

/// Per-CP pull weights toward the TFI anchor: each strict-interior
/// CP is weighted by its Greville point's proximity to the patch's
/// opposed samples (falloff radius 0.35 in UV), so the unfold only
/// moves geometry near the measured violations. Clean patches get
/// an all-zero field.
pub(super) fn opposed_proximity_weights(
    ring: &[RectangularGregoryBiNonic],
    references: &[(Vec3, f64)],
) -> Vec<[[f64; 10]; 10]> {
    const RADIUS: f64 = 0.35;
    let mut fields = Vec::with_capacity(ring.len());
    for (binonic, (reference, _)) in ring.iter().zip(references.iter()) {
        let normal = |u: f64, v: f64| binonic.partial_u(u, v).cross(&binonic.partial_v(u, v));
        let mut opposed: Vec<(f64, f64)> = Vec::new();
        for k in 0..=24 {
            let t = 0.02 + 0.96 * (k as f64) / 24.0;
            for (u, v) in [(0.02, t), (0.98, t), (t, 0.02), (t, 0.98)] {
                if normal(u, v).dot(reference) < 0.0 {
                    opposed.push((u, v));
                }
            }
        }
        for i in 0..9 {
            let u = 0.1 + 0.8 * (i as f64) / 8.0;
            for j in 0..9 {
                let v = 0.1 + 0.8 * (j as f64) / 8.0;
                if normal(u, v).dot(reference) < 0.0 {
                    opposed.push((u, v));
                }
            }
        }
        let mut field = [[0.0_f64; 10]; 10];
        if !opposed.is_empty() {
            for (i, row) in field.iter_mut().enumerate() {
                let gu = (i as f64) / 9.0;
                for (j, w) in row.iter_mut().enumerate() {
                    let gv = (j as f64) / 9.0;
                    let mut best = 0.0_f64;
                    for &(u, v) in &opposed {
                        let d = ((gu - u).powi(2) + (gv - v).powi(2)).sqrt();
                        best = best.max(1.0 - d / RADIUS);
                    }
                    *w = best.clamp(0.0, 1.0);
                }
            }
        }
        fields.push(field);
    }
    fields
}

/// Per-patch transfinite Coons (TFI) fill of the four boundary
/// rows: `P[i][j] = (1−a)·row0 + a·row9 + (1−b)·col0 + b·col9 −
/// bilinear(corners)` — the orientation-clean anchor whenever the
/// boundary itself is healthy.
pub(super) fn boundary_tfi_ring(
    ring: &[RectangularGregoryBiNonic],
) -> Vec<RectangularGregoryBiNonic> {
    let mut out = ring.to_vec();
    for patch in out.iter_mut() {
        let grid = patch.grid;
        for i in 1..=8 {
            let a = (i as f64) / 9.0;
            for j in 1..=8 {
                let b = (j as f64) / 9.0;
                let mix = |c: fn(&Point3) -> f64| -> f64 {
                    (1.0 - a) * c(&grid[0][j])
                        + a * c(&grid[9][j])
                        + (1.0 - b) * c(&grid[i][0])
                        + b * c(&grid[i][9])
                        - (1.0 - a) * (1.0 - b) * c(&grid[0][0])
                        - a * (1.0 - b) * c(&grid[9][0])
                        - (1.0 - a) * b * c(&grid[0][9])
                        - a * b * c(&grid[9][9])
                };
                patch.grid[i][j] = Point3::new(mix(|p| p.x), mix(|p| p.y), mix(|p| p.z));
            }
        }
    }
    out
}

pub(super) struct SeamPolishEval {
    pub base: Vec<RectangularGregoryBiNonic>,
    /// Per patch: fixed initial mid-face unit normal + its magnitude.
    pub references: Vec<(Vec3, f64)>,
    pub x0: Vec<f64>,
    pub outer_fillets: Option<Vec<SurfaceGeometry>>,
    pub char_len: f64,
    /// Lazily built 0/1 selector `∂ring_CP/∂x` for the seam-row
    /// parameterization (shared spoke blocks contribute TWO entries
    /// per column). Cached because it is constant across iterates
    /// (PK-275 mixed analytic Jacobian).
    pub selector_m: Option<Vec<Vec<f64>>>,
}

impl SeamPolishEval {
    fn build(&self, x: &[f64]) -> Vec<RectangularGregoryBiNonic> {
        let mut ring = self.base.clone();
        unpack_seam_rows(x, &mut ring);
        ring
    }

    /// The folded-G2 residual block, exactly as `residual` emits it
    /// (weighted, zero-padded on unreadable samples).
    fn g2_block(&self, ring: &[RectangularGregoryBiNonic]) -> Vec<f64> {
        let n = ring.len();
        let mut out = Vec::with_capacity(n * 25);
        for k in 0..n {
            let a = &ring[k];
            let b = &ring[(k + 1) % n];
            for sample in 0..25 {
                let t = 0.05 + 0.90 * (sample as f64) / 24.0;
                out.push(G2_WEIGHT * g2_residual_at(a, b, t).unwrap_or(0.0));
            }
        }
        out
    }
}

impl LmEval for SeamPolishEval {
    fn residual(&mut self, x: &[f64]) -> Option<Vec<f64>> {
        let ring = self.build(x);
        let mut r = coplanarity_residuals(&ring);
        // Folded curvature agreement across the spokes — the polish
        // must optimize what the G2 gate measures, or it sacrifices
        // curvature for coplanarity (unreadable samples contribute
        // zero so the vector length stays fixed).
        r.extend(self.g2_block(&ring));
        if let Some(fillets) = &self.outer_fillets {
            for (k, binonic) in ring.iter().enumerate() {
                let Some(fillet) = fillets.get(k) else {
                    r.extend(std::iter::repeat_n(0.0, 25));
                    continue;
                };
                for s in 0..25 {
                    let t = 0.05 + 0.90 * (s as f64) / 24.0;
                    r.push(outer_fillet_sin(binonic, fillet, t).unwrap_or(0.0));
                }
            }
        }
        for (xi, x0i) in x.iter().zip(self.x0.iter()) {
            r.push(TETHER_WEIGHT * (xi - x0i) / self.char_len.max(1e-12));
        }
        Some(r)
    }

    /// Feasibility = the orientation stays clean and the ring stays
    /// fold-free with margin — the polish can improve continuity
    /// only INSIDE the repaired region.
    fn fold_safe(&mut self, x: &[f64]) -> bool {
        let ring = self.build(x);
        if binonic_reversal_band_count(&ring) + interior_opposed_count(&ring) > 0 {
            return false;
        }
        // Cheap orientation re-check against the fixed references
        // guards the mid-face flip the band counter cannot see
        // (its reference is recomputed per candidate).
        for (binonic, (reference, _)) in ring.iter().zip(self.references.iter()) {
            let mid = binonic
                .partial_u(0.5, 0.5)
                .cross(&binonic.partial_v(0.5, 0.5));
            if mid.dot(reference) <= 0.0 {
                return false;
            }
        }
        !any_self_intersects(&ring, 1.0e-7)
    }

    /// Mixed analytic Jacobian (PK-275). Full-DoF FD cost one
    /// residual evaluation — including 25N NURBS fillet projections
    /// — per column per iteration (measured D = 504, 40 iterations,
    /// twice per repair on the extreme ladder). Per block:
    ///
    /// - coplanarity rows: analytic, shared ring Jacobian chained
    ///   through the constant seam-row selector;
    /// - folded-G2 rows: per-column FD of ONLY this block — the
    ///   rebuild is a CP overwrite and the samples are projection-
    ///   free arithmetic, so this keeps the branch-gated curvature
    ///   formula's exact FD behavior at negligible cost;
    /// - outer-fillet rows: analytic (footpoint chain), 6
    ///   projections per SAMPLE instead of one full projection
    ///   sweep per COLUMN;
    /// - tether rows: exact diagonal.
    fn jacobian(&mut self, x: &[f64]) -> Option<Vec<Vec<f64>>> {
        let n = self.base.len();
        let d = x.len();
        if self.selector_m.is_none() {
            self.selector_m = Some(seam_row_selector_m(n));
        }
        let m = self.selector_m.as_ref()?;
        let ring = self.build(x);

        let mut jac = analytic_ring_jacobian(&ring, m, d);

        // Folded-G2 block: stack-linear forward FD at the same eps
        // the full-FD driver used (LmConfig.fd_eps mirrors this
        // formula). The rebuild-based per-column FD re-evaluated
        // every seam sample's bidegree-9 derivatives D times per
        // iteration and was the measured binding term once the
        // Gregory SI gate was cut (PK-282); the Bézier stack is
        // linear in the CPs, so the same forward differences come
        // from closed-form stack deltas plus the pure gated formula.
        let eps = 1e-6 * self.char_len.max(1e-6);
        jac.extend(
            super::three_quad_chart_orientation_g2_jacobian::g2_jacobian_stack_linear(
                &ring, d, eps, G2_WEIGHT,
            ),
        );

        if let Some(fillets) = &self.outer_fillets {
            jac.extend(outer_fillet_jacobian_sampled(&ring, fillets, m, d, 25));
        }

        // Tether block: exact diagonal.
        let tether_scale = TETHER_WEIGHT / self.char_len.max(1e-12);
        for i in 0..d {
            let mut row = vec![0.0; d];
            row[i] = tether_scale;
            jac.push(row);
        }
        Some(jac)
    }
}

/// 0/1 selector matrix mapping the packed seam-row parameter vector
/// onto flat ring-CP coordinates, mirroring [`pack_seam_rows`]'s
/// order exactly: per-patch seam cells first, then the shared
/// per-spoke blocks — each spoke column carries TWO unit entries
/// because `unpack_seam_rows` writes `grid[9][j]` of patch `k` and
/// `grid[9−j][0]` of patch `k+1` from the same parameter.
#[cfg(test)]
#[path = "three_quad_chart_orientation_polish_fd_tests.rs"]
mod three_quad_chart_orientation_polish_fd_tests;

fn seam_row_selector_m(n: usize) -> Vec<Vec<f64>> {
    let cells = seam_row_cells();
    let d = n * cells.len() * 3 + n * 8 * 3;
    let mut m = vec![vec![0.0_f64; d]; n * 300];
    let mut k = 0;
    for p in 0..n {
        for &(i, j) in &cells {
            for axis in 0..3 {
                m[ring_idx(p, i, j, axis)][k] = 1.0;
                k += 1;
            }
        }
    }
    for p in 0..n {
        for j in 1..=8 {
            for axis in 0..3 {
                m[ring_idx(p, 9, j, axis)][k] = 1.0;
                m[ring_idx((p + 1) % n, 9 - j, 0, axis)][k] = 1.0;
                k += 1;
            }
        }
    }
    m
}

/// Per-sample spoke coplanarity, matching the promotion metric's
/// sampling so the no-regression pressure measures the same field.
fn coplanarity_residuals(ring: &[RectangularGregoryBiNonic]) -> Vec<f64> {
    let n = ring.len();
    let mut out = Vec::with_capacity(n * 50);
    for k in 0..n {
        let a = &ring[k];
        let b = &ring[(k + 1) % n];
        for sample in 0..50 {
            let t = 0.05 + 0.90 * (sample as f64) / 49.0;
            let t_b = 1.0 - t;
            let xtan_a = a.partial_u(1.0, t);
            let xtan_b = b.partial_v(t_b, 0.0);
            let t_spoke = a.partial_v(1.0, t);
            let numer = xtan_a.cross(&xtan_b).dot(&t_spoke).abs();
            let denom = xtan_a.length() * xtan_b.length() * t_spoke.length();
            out.push(if denom > 1e-300 { numer / denom } else { 0.0 });
        }
    }
    out
}

// ── Parameterization ─────────────────────────────────────────

/// The polish's per-patch degrees of freedom: the seam-adjacent
/// interior rows — `row1`/`row8` (all interior j) plus the
/// `col1`/`col8` interior cells not already covered — the rows that
/// own spoke coplanarity and outer tangency. Outer boundary rows
/// never enter the parameter vector; the spoke rows enter only via
/// the SHARED per-spoke blocks in [`pack_seam_rows`].
pub(super) fn seam_row_cells() -> Vec<(usize, usize)> {
    let mut cells = Vec::with_capacity(48);
    for j in 1..=8 {
        cells.push((1, j));
        cells.push((8, j));
    }
    for i in 2..=7 {
        cells.push((i, 1));
        cells.push((i, 8));
    }
    // The G2 rows: curvature agreement across the spokes and the
    // outer boundary is shaped by the second ring of CPs (measured
    // 1:50 s0.85: the unfolded candidate cleared G0/G1 and declined
    // at g2 5.36e-2 with these rows frozen).
    for j in 2..=7 {
        cells.push((2, j));
        cells.push((7, j));
    }
    for i in 3..=6 {
        cells.push((i, 2));
        cells.push((i, 7));
    }
    cells
}

/// Pack the polish DoF: the per-patch seam-adjacent interior rows
/// PLUS one shared block per spoke — the spoke curve's interior CPs
/// (`grid[9][j]` of patch k ≡ `grid[9−j][0]` of patch k+1). The
/// spoke interior is decomposition-authored free geometry (its
/// endpoints at C and the hex vertex stay fixed), and it is exactly
/// where a flat overfold must OPEN: with the spoke frozen, unfolding
/// trades the fold for a crease along it. Writing the shared block
/// to both patches keeps ring G0 exact by construction.
pub(super) fn pack_seam_rows(ring: &[RectangularGregoryBiNonic]) -> Vec<f64> {
    let cells = seam_row_cells();
    let mut x = Vec::with_capacity((ring.len() * cells.len() + ring.len() * 8) * 3);
    for patch in ring {
        for &(i, j) in &cells {
            x.push(patch.grid[i][j].x);
            x.push(patch.grid[i][j].y);
            x.push(patch.grid[i][j].z);
        }
    }
    for patch in ring {
        for j in 1..=8 {
            x.push(patch.grid[9][j].x);
            x.push(patch.grid[9][j].y);
            x.push(patch.grid[9][j].z);
        }
    }
    x
}

pub(super) fn unpack_seam_rows(x: &[f64], ring: &mut [RectangularGregoryBiNonic]) {
    let cells = seam_row_cells();
    let n = ring.len();
    let mut idx = 0;
    for patch in ring.iter_mut() {
        for &(i, j) in &cells {
            patch.grid[i][j] = Point3::new(x[idx], x[idx + 1], x[idx + 2]);
            idx += 3;
        }
    }
    for k in 0..n {
        for j in 1..=8 {
            let p = Point3::new(x[idx], x[idx + 1], x[idx + 2]);
            idx += 3;
            ring[k].grid[9][j] = p;
            ring[(k + 1) % n].grid[9 - j][0] = p;
        }
    }
}
