//! Sliding contacts between solids (PK-1073): a penalty on the normal
//! jump of the displacement across each shared patch, β ∫ (n·Δu)(n·Δv),
//! so the surfaces stay together along the normal and slide freely
//! (frictionless, no separation; `crate::contact`).
//!
//! β is [`SLIDING_PENALTY`] times the Nitsche scale (λ + 2μ)/h: stiff
//! enough that the contact's own compliance, a normal jump of σ/β under
//! a normal stress σ, stays a few thousandths of a cell's strain, and
//! soft enough for the multigrid solve. The jump is part of the answer,
//! not an error: under uniform compression each body's field stays
//! linear and the contact adds exactly σ/β.
//!
//! A group of coupled solids none of whose faces is held declines as
//! free; then each solid is checked for free motion on its block of the
//! coupled operator, its partners held: the contact restrains only its
//! normal, so a solid held through a contact alone reads free.

use crate::aggregate::Aggregation;
use crate::contact::{ContactPatch, Side, coupling_triplets};
use crate::element::ScalarTables;
use crate::error::AnalysisError;
use crate::grid::Grid;
use crate::linalg::sparse::Csr;
use crate::moments::SolidMoments;

use super::rigid::{RigidMotion, motions_below, rigid_spectrum};

/// The sliding penalty relative to the Nitsche scale (λ + 2μ)/h.
pub const SLIDING_PENALTY: f64 = 1e3;

fn internal(e: impl std::fmt::Display) -> AnalysisError {
    AnalysisError::Internal {
        stage: "sliding contact",
        detail: e.to_string(),
    }
}

/// `blocks` joined into one operator with every contact's coupling,
/// penalty `beta` (N/mm³).
pub(super) fn coupled(
    grid: &Grid,
    solids: &[(&SolidMoments, &Aggregation)],
    blocks: &[Csr],
    offsets: &[usize],
    contacts: &[ContactPatch],
    beta: f64,
) -> Result<Csr, AnalysisError> {
    let n: usize = blocks.iter().map(Csr::rows).sum();
    let mut triplets = Vec::new();
    for (block, &offset) in blocks.iter().zip(offsets) {
        for r in 0..block.rows() {
            for (c, v) in block.row(r) {
                triplets.push((offset + r, offset + c, v));
            }
        }
    }
    let reduce: Vec<Csr> = solids
        .iter()
        .map(|(_, g)| super::assemble::vector_extension(&g.extension).map_err(internal))
        .collect::<Result<_, _>>()?;
    let tables = ScalarTables::new();
    for patch in contacts {
        let side = |s: crate::model::SolidId| {
            let i = s.0 as usize;
            let ((moments, aggregation), reduce) =
                solids
                    .get(i)
                    .zip(reduce.get(i))
                    .ok_or(AnalysisError::Unsupported {
                        what: "a contact names a solid the model does not have",
                    })?;
            Ok::<_, AnalysisError>(Side {
                moments,
                extension: &aggregation.extension,
                reduce,
                offset: offsets[i],
            })
        };
        let (a, b) = (side(patch.solids.0)?, side(patch.solids.1)?);
        triplets.extend(coupling_triplets(grid, patch, (&a, &b), 3, |s, n| {
            let m = tables.facet_mass(s);
            let mut k = vec![vec![0.0; 24]; 24];
            for p in 0..8 {
                for q in 0..8 {
                    for i in 0..3 {
                        for j in 0..3 {
                            k[3 * p + i][3 * q + j] = beta * m[p][q] * n[i] * n[j];
                        }
                    }
                }
            }
            k
        })?);
    }
    Csr::from_triplets(n, n, &triplets).map_err(internal)
}

/// The solids of `count` joined by `contacts`, as groups of two or
/// more, each in ascending order.
pub(super) fn groups(count: usize, contacts: &[ContactPatch]) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..count).collect();
    fn root(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for c in contacts {
        let (a, b) = (c.solids.0.0 as usize, c.solids.1.0 as usize);
        if a < count && b < count {
            let (ra, rb) = (root(&mut parent, a), root(&mut parent, b));
            parent[ra.max(rb)] = ra.min(rb);
        }
    }
    let mut out: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
    for s in 0..count {
        let r = root(&mut parent, s);
        out.entry(r).or_default().push(s);
    }
    out.into_values().filter(|g| g.len() > 1).collect()
}

/// The free motions of solids `ranges` of the coupled `operator`, the
/// rest held: its rigid spectrum restricted to those dofs, at the free
/// nodes' `positions`, about `centre`.
pub(super) fn group_motions(
    operator: &Csr,
    ranges: &[(usize, usize)],
    positions: &[[f64; 3]],
    centre: [f64; 3],
) -> Result<Vec<RigidMotion>, AnalysisError> {
    let mut index = vec![usize::MAX; operator.rows()];
    let mut next = 0;
    for &(start, end) in ranges {
        for slot in &mut index[start..end] {
            *slot = next;
            next += 1;
        }
    }
    let mut triplets = Vec::new();
    for &(start, end) in ranges {
        for r in start..end {
            for (c, v) in operator.row(r) {
                if index[c] != usize::MAX {
                    triplets.push((index[r], index[c], v));
                }
            }
        }
    }
    let sub = Csr::from_triplets(next, next, &triplets).map_err(internal)?;
    let spectrum = rigid_spectrum(&sub, positions, centre).map_err(internal)?;
    Ok(motions_below(&spectrum, super::RIGID_MODE_TOL))
}

/// All six rigid motions about `centre`: what a group no face of which
/// is held is free to make (named in its decline, as the spectrum names
/// a single solid's).
pub(super) fn every_motion(centre: [f64; 3]) -> Vec<RigidMotion> {
    let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    axes.iter()
        .map(|&direction| RigidMotion::Translation { direction })
        .chain(axes.iter().map(|&direction| RigidMotion::Rotation {
            point: centre,
            direction,
        }))
        .collect()
}
