//! Geometric multigrid preconditioned conjugate gradients (map §1.5).
//!
//! Levels live on the 2× node lattice: a level-l dof sits at a finest
//! node whose coordinates are multiples of 2^l. Prolongation is
//! trilinear (weights 1 or ½ per axis). A coarse node gets one dof per
//! connected component of the finer dofs it interpolates to, using the
//! finer operator's graph, so separate solids and gaps wider than a
//! finer cell stay decoupled on every level (de Prenter 2019). Coarse
//! operators are Galerkin products PᵀAP; the coarsest is factored
//! densely. Each dof carries a label (the displacement component for
//! vector problems): coarse dofs group by (node, label), so components
//! never merge.

pub mod pcg;
pub mod smoother;

use crate::linalg::CholeskyError;
use crate::linalg::sparse::{Csr, SparseError};

use smoother::Chebyshev;

pub mod coarse;

/// Most levels in a hierarchy (map §1.5).
pub const MAX_LEVELS: usize = 6;

/// Largest connected block of the coarsest level factored densely:
/// 2,800² × 8 bytes = 62.7 MB. A hierarchy whose coarsening stalls with
/// a larger block declines instead of allocating gigabytes.
pub const COARSEST_HARD_CAP: usize = 2800;

// The cap keeps the dense coarsest factor under 64 MB.
const _: () = assert!(COARSEST_HARD_CAP * COARSEST_HARD_CAP * 8 < 64 * 1024 * 1024);

/// Stop coarsening when a level keeps more than this fraction of dofs.
pub const MIN_COARSENING: f64 = 0.7;

/// One multigrid level.
#[derive(Debug, Clone)]
pub struct Level {
    /// The operator.
    pub a: Csr,
    /// Finest-lattice node coordinates of each dof.
    pub positions: Vec<[usize; 3]>,
    /// Label of each dof (component); coarse dofs keep their label.
    pub labels: Vec<u8>,
    /// Prolongation from the next coarser level (absent on the
    /// coarsest).
    pub p: Option<Csr>,
    /// Restriction Pᵀ to the next coarser level.
    pub r: Option<Csr>,
    /// The smoother (absent on the coarsest).
    pub smoother: Option<Chebyshev>,
}

/// A multigrid hierarchy, finest level first.
#[derive(Debug, Clone)]
pub struct Hierarchy {
    /// Levels, finest first.
    pub levels: Vec<Level>,
    /// Dense factors of the coarsest operator's connected blocks.
    pub coarsest: coarse::BlockCholesky,
    /// Levels whose smoother fell back to damped Jacobi.
    pub jacobi_fallbacks: usize,
}

/// Why a hierarchy could not be built.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum MultigridError {
    /// A sparse product failed.
    #[error(transparent)]
    Sparse(#[from] SparseError),
    /// The coarsest operator is not positive definite.
    #[error("coarsest level: {0}")]
    Coarsest(#[from] CholeskyError),
    /// Coarsening stalled with the coarsest level above
    /// [`COARSEST_HARD_CAP`].
    #[error("coarsest level has a connected block of {dofs} dofs, above the cap of {cap}")]
    CoarsestTooLarge {
        /// Dofs of the coarsest level's largest connected block.
        dofs: usize,
        /// The cap.
        cap: usize,
    },
    /// Operator and positions disagree in size.
    #[error("{dofs} dofs but {positions} positions")]
    Positions {
        /// Operator size.
        dofs: usize,
        /// Positions given.
        positions: usize,
    },
}

/// Parameters of the hierarchy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HierarchyOptions {
    /// Coarsen until at most this many dofs remain.
    pub coarsest_dofs: usize,
    /// Split coarse dofs by connected component (always true in
    /// production; false only to show what splitting prevents).
    pub split_components: bool,
}

/// A prolongation with the coarser dofs' positions and labels.
type Prolongation = (Csr, Vec<[usize; 3]>, Vec<u8>);

/// Union-find with path halving and smallest-root union, so labels are
/// deterministic.
fn find(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

/// The prolongation from a coarser level to the level with operator `a`
/// and dof `positions` on lattice spacing `spacing`, and the coarser
/// dofs' positions.
fn prolongation(
    a: &Csr,
    positions: &[[usize; 3]],
    labels: &[u8],
    spacing: usize,
    split: bool,
) -> Result<Prolongation, SparseError> {
    let coarse = 2 * spacing;
    // ((coarse node, label), fine dof, weight), one per contribution.
    let mut contributions: Vec<(([usize; 3], u8), usize, f64)> = Vec::new();
    for (i, p) in positions.iter().enumerate() {
        let mut corners: Vec<([usize; 3], f64)> = vec![([0; 3], 1.0)];
        for axis in 0..3 {
            let base = p[axis] / coarse * coarse;
            let mut next = Vec::with_capacity(corners.len() * 2);
            for (c, w) in corners {
                if p[axis] == base {
                    let mut q = c;
                    q[axis] = base;
                    next.push((q, w));
                } else {
                    for node in [base, base + coarse] {
                        let mut q = c;
                        q[axis] = node;
                        next.push((q, 0.5 * w));
                    }
                }
            }
            corners = next;
        }
        for (q, w) in corners {
            contributions.push(((q, labels[i]), i, w));
        }
    }
    // Group by coarse node in lexicographic order (x fastest).
    contributions.sort_by_key(|&((q, l), i, _)| (q[2], q[1], q[0], l, i));
    let n = positions.len();
    let mut member = vec![usize::MAX; n];
    let mut parent: Vec<usize> = (0..n).collect();
    let mut triplets = Vec::new();
    let mut coarse_positions = Vec::new();
    let mut coarse_labels = Vec::new();
    let mut start = 0;
    while start < contributions.len() {
        let key = contributions[start].0;
        let end = start
            + contributions[start..]
                .iter()
                .take_while(|c| c.0 == key)
                .count();
        let group = &contributions[start..end];
        for (k, &(_, i, _)) in group.iter().enumerate() {
            member[i] = start + k;
            parent[i] = i;
        }
        if split {
            for &(_, i, _) in group {
                for (j, v) in a.row(i) {
                    if v != 0.0 && member[j] != usize::MAX && (start..end).contains(&member[j]) {
                        let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                        let (lo, hi) = (ri.min(rj), ri.max(rj));
                        parent[hi] = lo;
                    }
                }
            }
        }
        // Components labelled by their smallest fine dof.
        let mut roots: Vec<usize> = group
            .iter()
            .map(|&(_, i, _)| {
                if split {
                    find(&mut parent, i)
                } else {
                    group[0].1
                }
            })
            .collect();
        let mut labels_here = roots.clone();
        labels_here.sort_unstable();
        labels_here.dedup();
        let first = coarse_positions.len();
        coarse_positions.extend(std::iter::repeat_n(key.0, labels_here.len()));
        coarse_labels.extend(std::iter::repeat_n(key.1, labels_here.len()));
        for (&(_, i, w), root) in group.iter().zip(roots.iter_mut()) {
            let column = first + labels_here.partition_point(|&l| l < *root);
            triplets.push((i, column, w));
        }
        for &(_, i, _) in group {
            member[i] = usize::MAX;
        }
        start = end;
    }
    let p = Csr::from_triplets(n, coarse_positions.len(), &triplets)?;
    Ok((p, coarse_positions, coarse_labels))
}

impl Hierarchy {
    /// Builds the hierarchy for the SPD operator `a` whose dofs sit at
    /// finest-lattice `positions` with component `labels`.
    pub fn build(
        a: Csr,
        positions: Vec<[usize; 3]>,
        labels: Vec<u8>,
        options: &HierarchyOptions,
    ) -> Result<Hierarchy, MultigridError> {
        Hierarchy::rebuild(a, positions, labels, options, None)
    }

    /// [`Hierarchy::build`], taking each level's prolongation from
    /// `previous` where that level's operator has the same graph and
    /// zeros, positions and labels (the prolongation reads nothing
    /// else, so it is the same): a design update rescales the entries
    /// and keeps the graph.
    pub fn rebuild(
        a: Csr,
        positions: Vec<[usize; 3]>,
        labels: Vec<u8>,
        options: &HierarchyOptions,
        previous: Option<&Hierarchy>,
    ) -> Result<Hierarchy, MultigridError> {
        if a.rows() != positions.len() || labels.len() != positions.len() {
            return Err(MultigridError::Positions {
                dofs: a.rows(),
                positions: positions.len(),
            });
        }
        let mut levels = vec![Level {
            a,
            positions,
            labels,
            p: None,
            r: None,
            smoother: None,
        }];
        let mut spacing = 1;
        while levels.len() < MAX_LEVELS {
            let fine = levels.len() - 1;
            let n = levels[fine].a.rows();
            if n <= options.coarsest_dofs {
                break;
            }
            let reused = previous.and_then(|h| {
                let (old, next) = (h.levels.get(fine)?, h.levels.get(fine + 1)?);
                let same = old.a.same_graph(&levels[fine].a)
                    && old.positions == levels[fine].positions
                    && old.labels == levels[fine].labels;
                if !same {
                    return None;
                }
                Some((
                    (old.p.clone()?, next.positions.clone(), next.labels.clone()),
                    old.r.clone()?,
                ))
            });
            let ((p, coarse_positions, coarse_labels), r) = match reused {
                Some(found) => found,
                None => {
                    let prolonged = prolongation(
                        &levels[fine].a,
                        &levels[fine].positions,
                        &levels[fine].labels,
                        spacing,
                        options.split_components,
                    )?;
                    let r = prolonged.0.transpose();
                    (prolonged, r)
                }
            };
            if p.cols() as f64 > MIN_COARSENING * n as f64 {
                break;
            }
            let ac = r.mul(&levels[fine].a.mul(&p)?)?;
            levels[fine].p = Some(p);
            levels[fine].r = Some(r);
            levels.push(Level {
                a: ac,
                positions: coarse_positions,
                labels: coarse_labels,
                p: None,
                r: None,
                smoother: None,
            });
            spacing *= 2;
        }
        let mut jacobi_fallbacks = 0;
        let last = levels.len() - 1;
        for level in &mut levels[..last] {
            let s = Chebyshev::new(&level.a);
            if s.is_jacobi() {
                jacobi_fallbacks += 1;
            }
            level.smoother = Some(s);
        }
        let coarsest = match coarse::BlockCholesky::factor(&levels[last].a, COARSEST_HARD_CAP) {
            Ok(factored) => factored?,
            Err(dofs) => {
                return Err(MultigridError::CoarsestTooLarge {
                    dofs,
                    cap: COARSEST_HARD_CAP,
                });
            }
        };
        Ok(Hierarchy {
            levels,
            coarsest,
            jacobi_fallbacks,
        })
    }

    /// One V-cycle approximating A⁻¹ b on level `l` (zero initial guess).
    pub fn v_cycle(&self, l: usize, b: &[f64]) -> Result<Vec<f64>, MultigridError> {
        let level = &self.levels[l];
        let (Some(p), Some(r), Some(smoother)) = (&level.p, &level.r, &level.smoother) else {
            return Ok(self.coarsest.solve(b)?);
        };
        let mut x = vec![0.0; b.len()];
        smoother.smooth_from_zero(&level.a, b, &mut x)?;
        let mut ax = vec![0.0; b.len()];
        level.a.mul_vec(&x, &mut ax)?;
        let residual: Vec<f64> = b.iter().zip(&ax).map(|(bi, ai)| bi - ai).collect();
        let mut rc = vec![0.0; r.rows()];
        r.mul_vec(&residual, &mut rc)?;
        let ec = self.v_cycle(l + 1, &rc)?;
        let mut e = vec![0.0; b.len()];
        p.mul_vec(&ec, &mut e)?;
        for (xi, ei) in x.iter_mut().zip(&e) {
            *xi += ei;
        }
        smoother.smooth(&level.a, b, &mut x)?;
        Ok(x)
    }
}

#[cfg(test)]
#[path = "solve_tests.rs"]
mod solve_tests;
