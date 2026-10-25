//! A design's system re-assembled at each update from what does not
//! change: each entry's unscaled block (shared with the sensitivities)
//! and loads, where its block lands in the node matrix, and the
//! extension. An update only scales and adds the blocks, in the order
//! `assemble_solid` adds them, and reduces, so the system is the same
//! to the bit (PH-120: the blocks were 0.09–0.12 s of a 0.19 s
//! re-assembly at 22k dof).

use crate::error::AnalysisError;
use crate::linalg::sparse::Csr;

use super::super::assemble::{
    SolidInput, cell_operators, cell_positions, unscaled_terms, vector_extension,
};
use super::super::solve::Prepared;
use super::super::stability::entry_betas;
use super::super::{ElasticProblem, internal};
use super::blocks::{Block, Blocks};

/// One solid's fixed part of the system.
struct ScaledSolid {
    /// The node matrix's pattern, zero.
    pattern: Csr,
    /// Per entry, per (node a, node b, row i): the value index of
    /// columns 3b..3b + 3 of row 3a + i.
    slots: Vec<[u32; 192]>,
    /// Per entry, its nodes.
    nodes: Vec<[usize; 8]>,
    /// Per entry, its scaled and applied loads.
    loads: Vec<([f64; 24], [f64; 24])>,
    /// The extension E and Eᵀ.
    extension: (Csr, Csr),
}

/// Every solid's fixed part.
pub(super) struct ScaledSystem {
    solids: Vec<ScaledSolid>,
}

impl ScaledSystem {
    /// The fixed parts of `prepared`'s systems.
    pub(super) fn new(
        prepared: &Prepared,
        problem: &ElasticProblem,
    ) -> Result<ScaledSystem, AnalysisError> {
        let grid = &prepared.grid;
        let operators = cell_operators(&prepared.vt, grid.cell_size(), &prepared.elastic)
            .map_err(|e| internal("generative", e))?;
        let mut solids = Vec::with_capacity(prepared.solids.len());
        for (s, (_, clipped, moments, aggregation)) in prepared.solids.iter().enumerate() {
            let input = SolidInput {
                grid,
                clipped,
                moments,
                aggregation,
                faces: &problem.faces,
                prescribed: None,
                thermal: prepared.strains[s].as_ref(),
                cores: prepared.cores[s].as_deref(),
                scales: None,
            };
            let e = &aggregation.extension;
            let mut nodes = Vec::with_capacity(moments.cells.len());
            for entry in 0..moments.cells.len() {
                nodes.push(
                    cell_positions(grid, moments, e, entry)
                        .map_err(|e| internal("generative", e))?,
                );
            }
            let pattern = Csr::block_pattern(e.nodes.len(), 3, &nodes)
                .map_err(|e| internal("generative", e))?;
            let slot = |n: &[usize; 8]| -> Result<[u32; 192], AnalysisError> {
                let mut out = [0u32; 192];
                for (a, &na) in n.iter().enumerate() {
                    for (b, &nb) in n.iter().enumerate() {
                        for i in 0..3 {
                            let at = pattern
                                .position(3 * na + i, 3 * nb)
                                .filter(|&k| {
                                    pattern.position(3 * na + i, 3 * nb + 2) == Some(k + 2)
                                })
                                .ok_or_else(|| {
                                    internal("generative", "block outside the pattern")
                                })?;
                            out[24 * a + 3 * b + i] = u32::try_from(at)
                                .map_err(|_| internal("generative", "pattern too large"))?;
                        }
                    }
                }
                Ok(out)
            };
            let slots = nodes.iter().map(slot).collect::<Result<Vec<_>, _>>()?;
            let betas = entry_betas(&operators, &prepared.bt, &input, &prepared.elastic)
                .map_err(|e| internal("generative", e))?;
            let load = |entry: usize| {
                unscaled_terms(
                    &operators,
                    &prepared.bt,
                    &input,
                    &prepared.elastic,
                    &betas,
                    entry,
                )
                .map(|(_, fs, fl)| (fs, fl))
                .map_err(|e| internal("generative", e))
            };
            let mut loads = Vec::with_capacity(moments.cells.len());
            for l in crate::parallel::map_range(0..moments.cells.len(), 256, load) {
                loads.push(l?);
            }
            let ev = vector_extension(e).map_err(|e| internal("generative", e))?;
            let et = ev.transpose();
            solids.push(ScaledSolid {
                pattern,
                slots,
                nodes,
                loads,
                extension: (ev, et),
            });
        }
        Ok(ScaledSystem { solids })
    }

    /// Each solid's reduced matrix and the stacked load, with each
    /// entry's block (from `blocks`) and scaled load times its scale
    /// in `scales` (per solid, per entry).
    pub(super) fn assemble(
        &self,
        blocks: &Blocks,
        scales: &[Vec<f64>],
    ) -> Result<(Vec<Csr>, Vec<f64>), AnalysisError> {
        let mut matrices = Vec::with_capacity(self.solids.len());
        let mut b = Vec::new();
        for (s, solid) in self.solids.iter().enumerate() {
            let mut a = solid.pattern.clone();
            let mut load = vec![0.0; solid.pattern.rows()];
            let values = a.values_mut();
            for (entry, (nodes, slots)) in solid.nodes.iter().zip(&solid.slots).enumerate() {
                let scale = scales
                    .get(s)
                    .and_then(|v| v.get(entry))
                    .copied()
                    .unwrap_or(1.0);
                let (fs, fl) = &solid.loads[entry];
                for r in 0..24 {
                    load[3 * nodes[r / 3] + r % 3] += scale * fs[r] + fl[r];
                }
                let k: &Block = blocks.own[s]
                    .get(entry)
                    .and_then(|b| b.as_deref())
                    .unwrap_or(&blocks.full);
                for a in 0..8 {
                    for bn in 0..8 {
                        for i in 0..3 {
                            let at = slots[24 * a + 3 * bn + i] as usize;
                            for j in 0..3 {
                                values[at + j] += scale * k[3 * a + i][3 * bn + j];
                            }
                        }
                    }
                }
            }
            let (ev, et) = &solid.extension;
            let reduced = et
                .mul(&a.mul(ev).map_err(|e| internal("assemble", e))?)
                .map_err(|e| internal("assemble", e))?;
            let mut bs = vec![0.0; ev.cols()];
            et.mul_vec(&load, &mut bs)
                .map_err(|e| internal("assemble", e))?;
            matrices.push(reduced);
            b.extend(bs);
        }
        Ok((matrices, b))
    }
}

#[cfg(test)]
#[path = "scaled_tests.rs"]
mod tests;
