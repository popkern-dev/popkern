//! The spoke rows' null space in the spokes' structure (PK-620).
//!
//! A free coefficient is
//!
//! - **interior** when no spoke row touches it: its own null-space
//!   direction, the unit vector;
//! - on spoke p's **strip** when only spoke p's rows touch it (the rows
//!   beside the spoke, away from the centre);
//! - a **separator** when two or more spokes' rows touch it (the corner
//!   blocks at the centre, which every spoke reaches through its two
//!   sectors).
//!
//! Per spoke, with `C_S` its rows over its strip and `C_Q` over the
//! separators, the strip's own null space `N_p = null(C_S)` gives spoke
//! p's **local modes**: they move only that spoke's strip. A separator
//! motion `x_Q` is admissible when every spoke's rows can absorb it,
//! `Ũ_pᵀ C_Q x_Q = 0` with `Ũ_p` spanning `C_S`'s left null space; its
//! strips then follow by the pseudo-inverse, `x_S = −C_S⁺ C_Q x_Q`. The
//! admissible separator motions with their strip responses are the
//! **global** block, orthonormalised.
//!
//! Interior and local-mode columns are unit or orthonormal on disjoint
//! supports, and the responses lie in each `C_S`'s row space, orthogonal
//! to `N_p`: the basis is orthonormal, so every quantity the fit takes
//! from `Z` (the energy's normalisation, the warm-start projection) is
//! the dense basis's. Measured at n = 16, degree 9 with corner cells,
//! per coordinate: 640 interior, 240 local-mode and 38 global columns,
//! against a dense 278-column block over every coupled coefficient.

use std::ops::Range;

use faer::Mat;

use super::linalg::svd;
use super::nullspace::RANK_TOL;
use super::space::{Functional, SectorSpace};

/// A free coefficient's basis row: `(first column, length, offset into
/// the values)` segments.
pub(super) type BasisRow = Vec<(usize, usize, usize)>;

/// The null-space columns' blocks, per coordinate.
#[derive(Debug, Clone)]
pub(super) struct ColumnLayout {
    /// Each sector's interior columns.
    pub interior: Vec<Range<usize>>,
    /// Each spoke's local modes (spoke p lies between sectors p and
    /// p + 1).
    pub modes: Vec<Range<usize>>,
    /// The global block.
    pub global: Range<usize>,
    /// Null-space dimension.
    pub k: usize,
}

/// The basis rows and the layout, with the particular solve.
pub(super) struct NestedBasis {
    pub layout: ColumnLayout,
    pub rows: Vec<BasisRow>,
    pub values: Vec<f64>,
    pub particular: ParticularSolve,
}

/// A solution of the spoke rows `C x = b` in the same structure: the
/// separators from the stacked compatibility rows by their
/// pseudo-inverse, then each strip from its own (PK-620; the dense SVD
/// of every coupled coefficient this replaces cost 3.3 s at n = 16,
/// degree 9).
#[derive(Debug)]
pub(super) struct ParticularSolve {
    spokes: Vec<SpokeParticular>,
    separators: Vec<usize>,
    /// `K⁺ = k_v k_uᵀ` (`k_v` scaled by `1/σ`).
    k_u: Mat<f64>,
    k_v: Mat<f64>,
}

/// One spoke's part of the particular solve.
#[derive(Debug)]
struct SpokeParticular {
    /// The spoke's rows in the spoke row list.
    rows: Range<usize>,
    strip: Vec<usize>,
    /// `C_S⁺ = pinv_v pinv_uᵀ` (`pinv_v` scaled by `1/σ`).
    pinv_u: Mat<f64>,
    pinv_v: Mat<f64>,
    /// `C_S`'s left null space (rows × nullity).
    left_null: Mat<f64>,
    /// The rows over the separators.
    c_q: Mat<f64>,
}

impl ParticularSolve {
    /// A solution of the spoke rows with right-hand side `b` (one entry
    /// per spoke row) over the free rows (`nfree` of them); interior
    /// coefficients 0. When `b` is outside the rows' range this is the
    /// separators' least-squares choice.
    pub fn solve(&self, b: &[f64], nfree: usize) -> Vec<f64> {
        let mut e = Vec::new();
        for spoke in &self.spokes {
            let bp = &b[spoke.rows.clone()];
            for j in 0..spoke.left_null.ncols() {
                e.push(
                    (0..bp.len())
                        .map(|t| spoke.left_null.read(t, j) * bp[t])
                        .sum::<f64>(),
                );
            }
        }
        let rank = self.k_u.ncols();
        let coef: Vec<f64> = (0..rank)
            .map(|j| (0..e.len()).map(|i| self.k_u.read(i, j) * e[i]).sum())
            .collect();
        let xq: Vec<f64> = (0..self.separators.len())
            .map(|q| (0..rank).map(|j| self.k_v.read(q, j) * coef[j]).sum())
            .collect();
        let mut x = vec![0.0; nfree];
        for (q, &r) in self.separators.iter().enumerate() {
            x[r] = xq[q];
        }
        for spoke in &self.spokes {
            let bp = &b[spoke.rows.clone()];
            let y: Vec<f64> = (0..bp.len())
                .map(|t| {
                    bp[t]
                        - (0..xq.len())
                            .map(|q| spoke.c_q.read(t, q) * xq[q])
                            .sum::<f64>()
                })
                .collect();
            let rank = spoke.pinv_u.ncols();
            let coef: Vec<f64> = (0..rank)
                .map(|j| (0..y.len()).map(|t| spoke.pinv_u.read(t, j) * y[t]).sum())
                .collect();
            for (a, &r) in spoke.strip.iter().enumerate() {
                x[r] = (0..rank).map(|j| spoke.pinv_v.read(a, j) * coef[j]).sum();
            }
        }
        x
    }
}

/// Where a free coefficient sits.
#[derive(Clone, Copy, PartialEq)]
enum Role {
    Interior,
    Strip(usize),
    Separator,
}

/// One spoke's pieces: its strip (free rows), local modes, and the
/// strip response to a separator motion.
struct SpokeSolve {
    strip: Vec<usize>,
    modes: Mat<f64>,
    /// `−C_S⁺ C_Q` (strip × separators).
    response: Mat<f64>,
    particular: SpokeParticular,
}

/// The nested basis; `None` when the spoke rows do not group by spoke
/// or an SVD workspace cannot be sized.
pub(super) fn nested_basis(
    space: &SectorSpace,
    free: &[usize],
    slots: &[(bool, usize)],
    spoke_rows: &[Functional],
) -> Option<NestedBasis> {
    let n = space.n;
    if n == 0 || !spoke_rows.len().is_multiple_of(n) {
        return None;
    }
    let per_spoke = spoke_rows.len() / n;
    let mut touching: Vec<Vec<usize>> = vec![Vec::new(); free.len()];
    for (index, row) in spoke_rows.iter().enumerate() {
        let spoke = index / per_spoke;
        for &(i, w) in row {
            if w != 0.0
                && let (false, r) = slots[i]
                && !touching[r].contains(&spoke)
            {
                touching[r].push(spoke);
            }
        }
    }
    let role: Vec<Role> = touching
        .iter()
        .map(|t| match t.as_slice() {
            [] => Role::Interior,
            [p] => Role::Strip(*p),
            _ => Role::Separator,
        })
        .collect();
    let separators: Vec<usize> = (0..free.len())
        .filter(|&r| role[r] == Role::Separator)
        .collect();
    let mut separator_col = vec![usize::MAX; free.len()];
    for (q, &r) in separators.iter().enumerate() {
        separator_col[r] = q;
    }
    let nq_all = separators.len();
    let mut spokes = Vec::with_capacity(n);
    let mut compatibility: Vec<Vec<f64>> = Vec::new();
    for p in 0..n {
        let strip: Vec<usize> = (0..free.len())
            .filter(|&r| role[r] == Role::Strip(p))
            .collect();
        let mut strip_col = std::collections::HashMap::new();
        for (a, &r) in strip.iter().enumerate() {
            strip_col.insert(r, a);
        }
        let rows = &spoke_rows[p * per_spoke..(p + 1) * per_spoke];
        let mut c_s = Mat::<f64>::zeros(rows.len(), strip.len());
        let mut c_q = Mat::<f64>::zeros(rows.len(), nq_all);
        for (t, row) in rows.iter().enumerate() {
            for &(i, w) in row {
                let (false, r) = slots[i] else { continue };
                if let Some(&a) = strip_col.get(&r) {
                    c_s.write(t, a, c_s.read(t, a) + w);
                } else if separator_col[r] != usize::MAX {
                    let q = separator_col[r];
                    c_q.write(t, q, c_q.read(t, q) + w);
                }
            }
        }
        let row_range = p * per_spoke..(p + 1) * per_spoke;
        if strip.is_empty() {
            // Every row constrains the separators directly.
            for t in 0..rows.len() {
                compatibility.push((0..nq_all).map(|q| c_q.read(t, q)).collect());
            }
            let identity = Mat::from_fn(rows.len(), rows.len(), |i, j| f64::from(u8::from(i == j)));
            spokes.push(SpokeSolve {
                strip: Vec::new(),
                modes: Mat::zeros(0, 0),
                response: Mat::zeros(0, nq_all),
                particular: SpokeParticular {
                    rows: row_range,
                    strip,
                    pinv_u: Mat::zeros(rows.len(), 0),
                    pinv_v: Mat::zeros(0, 0),
                    left_null: identity,
                    c_q,
                },
            });
            continue;
        }
        let right = svd(&c_s)?;
        let s_max = right.s.first().copied().unwrap_or(0.0);
        let rank = right.s.iter().filter(|&&s| s > RANK_TOL * s_max).count();
        let modes = Mat::from_fn(strip.len(), strip.len() - rank, |a, j| {
            right.v.read(a, rank + j)
        });
        // −C_S⁺ C_Q = −V_ρ Σ⁻¹ U_ρᵀ C_Q.
        let mut ut_cq = Mat::<f64>::zeros(rank, nq_all);
        for j in 0..rank {
            for q in 0..nq_all {
                let v: f64 = (0..rows.len())
                    .map(|t| right.u.read(t, j) * c_q.read(t, q))
                    .sum();
                ut_cq.write(j, q, v / right.s[j]);
            }
        }
        let response = Mat::from_fn(strip.len(), nq_all, |a, q| {
            -(0..rank)
                .map(|j| right.v.read(a, j) * ut_cq.read(j, q))
                .sum::<f64>()
        });
        // C_S's left null space: the full left basis from the
        // transpose's right vectors.
        let left = svd(&c_s.transpose().to_owned())?;
        let left_null = Mat::from_fn(rows.len(), rows.len() - rank, |t, j| {
            left.v.read(t, rank + j)
        });
        let pinv_u = Mat::from_fn(rows.len(), rank, |t, j| right.u.read(t, j));
        let pinv_v = Mat::from_fn(strip.len(), rank, |a, j| right.v.read(a, j) / right.s[j]);
        for j in rank..rows.len() {
            compatibility.push(
                (0..nq_all)
                    .map(|q| {
                        (0..rows.len())
                            .map(|t| left.v.read(t, j) * c_q.read(t, q))
                            .sum()
                    })
                    .collect(),
            );
        }
        spokes.push(SpokeSolve {
            strip: strip.clone(),
            modes,
            response,
            particular: SpokeParticular {
                rows: row_range,
                strip,
                pinv_u,
                pinv_v,
                left_null,
                c_q,
            },
        });
    }
    let (admissible, k_u, k_v) = admissible_separators(&compatibility, nq_all)?;
    let global = orthonormal_global(&admissible, &spokes, &separators, free.len())?;
    let nq = global.ncols();
    // Layout and rows.
    let per = space.per();
    let mut rows: Vec<BasisRow> = vec![Vec::new(); free.len()];
    let mut values = vec![1.0];
    let mut column = 0;
    let mut interior = Vec::with_capacity(n);
    for s in 0..n {
        let first = column;
        for (r, &i) in free.iter().enumerate() {
            if i / per == s && role[r] == Role::Interior {
                rows[r].push((column, 1, 0));
                column += 1;
            }
        }
        interior.push(first..column);
    }
    let mut modes = Vec::with_capacity(n);
    for spoke in &spokes {
        let count = spoke.modes.ncols();
        for (a, &r) in spoke.strip.iter().enumerate() {
            if count > 0 {
                rows[r].push((column, count, values.len()));
                values.extend((0..count).map(|j| spoke.modes.read(a, j)));
            }
        }
        modes.push(column..column + count);
        column += count;
    }
    let global_range = column..column + nq;
    if nq > 0 {
        for r in 0..free.len() {
            if role[r] != Role::Interior {
                rows[r].push((column, nq, values.len()));
                values.extend((0..nq).map(|g| global.read(r, g)));
            }
        }
    }
    column += nq;
    let particular = ParticularSolve {
        spokes: spokes.into_iter().map(|s| s.particular).collect(),
        separators,
        k_u,
        k_v,
    };
    Some(NestedBasis {
        particular,
        layout: ColumnLayout {
            interior,
            modes,
            global: global_range,
            k: column,
        },
        rows,
        values,
    })
}

/// The admissible separator motions, the null space of the stacked
/// compatibility rows `K` (separators × nq), with `K⁺`'s factors.
fn admissible_separators(
    compatibility: &[Vec<f64>],
    nq_all: usize,
) -> Option<(Mat<f64>, Mat<f64>, Mat<f64>)> {
    if compatibility.is_empty() {
        let identity = Mat::from_fn(nq_all, nq_all, |i, j| f64::from(u8::from(i == j)));
        return Some((identity, Mat::zeros(0, 0), Mat::zeros(nq_all, 0)));
    }
    let k = Mat::from_fn(compatibility.len(), nq_all, |i, j| compatibility[i][j]);
    let d = svd(&k)?;
    let s_max = d.s.first().copied().unwrap_or(0.0);
    let rank = d.s.iter().filter(|&&s| s > RANK_TOL * s_max).count();
    let null = Mat::from_fn(nq_all, nq_all - rank, |i, j| d.v.read(i, rank + j));
    let k_u = Mat::from_fn(compatibility.len(), rank, |i, j| d.u.read(i, j));
    let k_v = Mat::from_fn(nq_all, rank, |i, j| d.v.read(i, j) / d.s[j]);
    Some((null, k_u, k_v))
}

/// The global columns over the free rows: each admissible separator
/// motion with every strip's response, orthonormalised (modified
/// Gram–Schmidt, twice).
fn orthonormal_global(
    admissible: &Mat<f64>,
    spokes: &[SpokeSolve],
    separators: &[usize],
    nfree: usize,
) -> Option<Mat<f64>> {
    let nq = admissible.ncols();
    let mut z = Mat::<f64>::zeros(nfree, nq);
    for g in 0..nq {
        for (q, &r) in separators.iter().enumerate() {
            z.write(r, g, admissible.read(q, g));
        }
        for spoke in spokes {
            for (a, &r) in spoke.strip.iter().enumerate() {
                let v: f64 = (0..separators.len())
                    .map(|q| spoke.response.read(a, q) * admissible.read(q, g))
                    .sum();
                z.write(r, g, v);
            }
        }
    }
    for _ in 0..2 {
        for g in 0..nq {
            for h in 0..g {
                let dot: f64 = (0..nfree).map(|r| z.read(r, g) * z.read(r, h)).sum();
                for r in 0..nfree {
                    z.write(r, g, z.read(r, g) - dot * z.read(r, h));
                }
            }
            let norm = (0..nfree).map(|r| z.read(r, g).powi(2)).sum::<f64>().sqrt();
            if norm.is_nan() || norm <= 1e-12 {
                return None;
            }
            for r in 0..nfree {
                z.write(r, g, z.read(r, g) / norm);
            }
        }
    }
    Some(z)
}
