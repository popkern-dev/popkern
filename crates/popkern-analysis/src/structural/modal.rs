//! Natural frequencies by shift-invert Lanczos (map Stage 2, item 4).
//!
//! K φ = λ M φ with K the condensed stiffness and M the consistent mass
//! ρ∫N_aN_b per component, both on the same moments and reduced by the
//! same AgFEM extension (E is linear, so a rigid motion maps exactly).
//! Lanczos (`lanczos`) runs on (K − σM)⁻¹M, grown until the modes
//! converge and restarted, deflated, until no repeated mode is missing
//! (`modal_restarts`). A supported part takes σ = 0; a free part (no
//! support: modal analysis of a loose body) takes σ < 0, and its six
//! rigid modes per solid are locked exactly before the runs (one start
//! vector would see one direction of them). λ is in N/(mm·kg):
//! ω² = 10³λ.

use crate::element::ScalarTables;
use crate::error::AnalysisError;
use crate::linalg::sparse::Csr;
use crate::model::AnalysisModel;
use crate::options::AnalysisOptions;
use crate::report::SolveReport;
use crate::solve::pcg::pcg;
use crate::solve::{Hierarchy, HierarchyOptions};
use crate::thermal::sample::point_value;

use super::assemble::{FaceMechanics, vector_extension};
use super::modal_restarts::{Locked, MAX_RESTARTS, lowest_modes};
use super::solve::prepare;
use super::{ElasticProblem, internal};

/// Lanczos residual (relative, `lanczos`) a mode must reach.
pub const MODE_TOLERANCE: f64 = 1.0e-6;

/// One mode.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct Mode {
    /// Natural frequency, Hz (0 for a rigid mode; λ below zero rounds to
    /// 0).
    pub frequency_hz: f64,
    /// ω², 1/s².
    pub omega_squared: f64,
    /// Lanczos's relative residual (bounds the relative error of ω² by
    /// its square over the gap to the next mode).
    pub residual: f64,
}

/// The modes of a part.
#[derive(Debug, Clone, PartialEq)]
pub struct ModalResult {
    /// The lowest modes, ascending.
    pub modes: Vec<Mode>,
    /// Each mode's displacement at the mesh vertices, per solid
    /// (mass-normalised, mm per √kg).
    pub shapes: Vec<Vec<Vec<[f64; 3]>>>,
    /// Total mass, kg.
    pub mass_kg: f64,
    /// What the solve did.
    pub report: SolveReport,
}

/// The `count` lowest modes of `problem` (its loads ignored) on `model`
/// with density `density` (kg/m³).
pub fn solve_modal(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    density: f64,
    count: usize,
    options: &AnalysisOptions,
) -> Result<ModalResult, AnalysisError> {
    solve_modal_with(model, problem, density, count, None, options)
}

/// [`solve_modal`] with an infill: the core's stiffness as the static
/// solve has it, and its mass the relative density's share.
pub fn solve_modal_with(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    density: f64,
    count: usize,
    infill: Option<&super::infill::Infill>,
    options: &AnalysisOptions,
) -> Result<ModalResult, AnalysisError> {
    solve_modal_capped(
        model,
        problem,
        density,
        (count, ModalCaps::default()),
        infill,
        options,
    )
}

/// Limits of a modal solve; the defaults are production's, the others
/// exist for the negative controls.
#[derive(Debug, Clone, Copy)]
pub(crate) struct ModalCaps {
    /// Lanczos steps per run.
    pub max_steps: usize,
    /// Deflated restarts after the first run.
    pub restarts: usize,
    /// Lock a free part's rigid modes before the runs.
    pub lock_rigid: bool,
}

impl Default for ModalCaps {
    fn default() -> Self {
        Self {
            max_steps: MAX_STEPS,
            restarts: MAX_RESTARTS,
            lock_rigid: true,
        }
    }
}

/// [`solve_modal_with`] under `caps`.
pub(crate) fn solve_modal_capped(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    density: f64,
    (count, caps): (usize, ModalCaps),
    infill: Option<&super::infill::Infill>,
    options: &AnalysisOptions,
) -> Result<ModalResult, AnalysisError> {
    if !(density.is_finite() && density > 0.0) || count == 0 {
        return Err(AnalysisError::Unsupported {
            what: "modal analysis needs a positive density and at least one mode",
        });
    }
    let supported = problem.faces.iter().any(|f| {
        matches!(
            f,
            FaceMechanics::Support(_)
                | FaceMechanics::Prescribed
                | FaceMechanics::PrescribedTraction
        )
    });
    let mut prepared = prepare(
        model,
        problem,
        options,
        (crate::moments::Splitting::Exact, !supported),
        &super::solve::Extras {
            infill,
            ..super::solve::Extras::default()
        },
    )?;
    // Mass per solid: ρ (kg/mm³) ∫ N_a N_b, per component, reduced.
    let rho = 1e-9 * density;
    let tables = ScalarTables::new();
    let mut mass_blocks = Vec::new();
    let solid_cores = prepared.cores.clone();
    for (s, (_, _, moments, aggregation)) in prepared.solids.iter().enumerate() {
        let e = &aggregation.extension;
        let cores = solid_cores.get(s).and_then(|c| c.as_deref());
        let mut triplets = Vec::new();
        for (entry, cell) in moments.cells.iter().enumerate() {
            let nodes = moments
                .node_keys(&prepared.grid, entry)
                .unwrap_or_default()
                .map(|n| e.position(n).unwrap_or(usize::MAX));
            if nodes.contains(&usize::MAX) {
                return Err(internal("mass", "cell node missing from the extension"));
            }
            // The core weighs its relative density's share.
            let share = infill.map_or(1.0, |f| f.density);
            let mass_moments: [f64; 27] = match cores.and_then(|c| c.get(entry)) {
                Some(crate::element::operators::Core::Full) => cell.m.map(|v| share * v),
                Some(crate::element::operators::Core::Partial(mc)) => {
                    std::array::from_fn(|i| cell.m[i] - (1.0 - share) * mc[i])
                }
                _ => cell.m,
            };
            let m = tables.mass(&mass_moments);
            for a in 0..8 {
                for b in 0..8 {
                    for i in 0..3 {
                        triplets.push((3 * nodes[a] + i, 3 * nodes[b] + i, rho * m[a][b]));
                    }
                }
            }
        }
        let n = 3 * e.nodes.len();
        let m_nodes = Csr::from_triplets(n, n, &triplets).map_err(|e| internal("mass", e))?;
        let ev = vector_extension(e).map_err(|e| internal("mass", e))?;
        mass_blocks.push(
            ev.transpose()
                .mul(&m_nodes.mul(&ev).map_err(|e| internal("mass", e))?)
                .map_err(|e| internal("mass", e))?,
        );
    }
    let k = Csr::block_diagonal(std::mem::take(&mut prepared.blocks))
        .map_err(|e| internal("assemble", e))?;
    let m = Csr::block_diagonal(mass_blocks).map_err(|e| internal("mass", e))?;
    let n = k.rows();
    prepared.report.dof = n;
    let ones_x: Vec<f64> = (0..n).map(|i| if i % 3 == 0 { 1.0 } else { 0.0 }).collect();
    let mut m_ones = vec![0.0; n];
    m.mul_vec(&ones_x, &mut m_ones)
        .map_err(|e| internal("mass", e))?;
    let mass_kg: f64 = ones_x.iter().zip(&m_ones).map(|(a, b)| a * b).sum();
    // σ: 0 for a supported part; below zero for a free one, small
    // against the mean stiffness-to-mass ratio of the diagonals.
    let (kd, md) = (k.diagonal(), m.diagonal());
    let ratio = kd.iter().sum::<f64>() / md.iter().sum::<f64>().max(f64::MIN_POSITIVE);
    let sigma = if supported { 0.0 } else { -1.0e-4 * ratio };
    let mut triplets = Vec::new();
    for r in 0..n {
        triplets.extend(k.row(r).map(|(c, v)| (r, c, v)));
        triplets.extend(m.row(r).map(|(c, v)| (r, c, -sigma * v)));
    }
    let shifted = Csr::from_triplets(n, n, &triplets).map_err(|e| internal("modal", e))?;
    // A free part's rigid motions, per solid, before the positions move
    // into the hierarchy.
    let ModalCaps {
        max_steps,
        restarts,
        lock_rigid,
    } = caps;
    let rigid = if supported || !lock_rigid {
        Vec::new()
    } else {
        rigid_vectors(&prepared.positions, &prepared.offsets, n)
    };
    let hierarchy = Hierarchy::build(
        shifted,
        std::mem::take(&mut prepared.positions),
        std::mem::take(&mut prepared.labels),
        &HierarchyOptions {
            coarsest_dofs: options.coarsest_dofs,
            split_components: true,
        },
    )
    .map_err(|e| prepared.report.singular(e.to_string()))?;
    prepared.report.levels = hierarchy.levels.len();
    // The rigid modes locked, then deflated Lanczos runs until a
    // restart finds nothing new (`modal_restarts`).
    let locked = rigid_modes(&m, &hierarchy, (sigma, options.cg_max_iters), rigid)?;
    let (found, spent) = lowest_modes(
        &m,
        &hierarchy,
        (sigma, options.cg_max_iters),
        (count, max_steps, restarts),
        locked,
    )?;
    if found.len() < count.min(n) {
        return Err(AnalysisError::SolverNotConverged {
            history: found.iter().map(|l| l.1).collect(),
        });
    }
    let modes: Vec<Mode> = found.iter().map(|l| mode(l.0, l.1)).collect();
    let vectors: Vec<Vec<f64>> = found.into_iter().map(|l| l.2).collect();
    prepared.report.lanczos_steps = spent.steps;
    prepared.report.lanczos_runs = spent.runs;
    prepared.report.cg_iterations = spent.iterations;
    let mut shapes = Vec::new();
    for x in &vectors {
        shapes.push(vertex_shapes(model, &prepared, x)?);
    }
    Ok(ModalResult {
        modes,
        shapes,
        mass_kg,
        report: prepared.report,
    })
}

/// Most Lanczos steps (the basis is kept whole for reorthogonalisation).
pub const MAX_STEPS: usize = 480;

/// The mode of eigenvalue `lambda` (N/(mm·kg)) with Lanczos residual
/// `residual`.
fn mode(lambda: f64, residual: f64) -> Mode {
    let omega_squared = 1e3 * lambda;
    Mode {
        frequency_hz: omega_squared.max(0.0).sqrt() / (2.0 * std::f64::consts::PI),
        omega_squared,
        residual,
    }
}

/// Each solid's six rigid motions on the free dofs (three per node at
/// its lattice position, solid `s` from `offsets[s]`): translations and
/// rotations about the solid's mean node. The extension is linear, so
/// they are the exact rigid motions of the part.
fn rigid_vectors(positions: &[[usize; 3]], offsets: &[usize], n: usize) -> Vec<Vec<f64>> {
    let mut out = Vec::new();
    for (s, &start) in offsets.iter().enumerate() {
        let end = offsets
            .get(s + 1)
            .copied()
            .unwrap_or(n)
            .min(positions.len());
        let nodes = (end.saturating_sub(start)) / 3;
        if nodes == 0 {
            continue;
        }
        let mut mean = [0.0; 3];
        for p in positions[start..end].iter().step_by(3) {
            for a in 0..3 {
                mean[a] += p[a] as f64 / nodes as f64;
            }
        }
        for motion in 0..6 {
            let mut v = vec![0.0; n];
            for dof in start..end {
                let c = (dof - start) % 3;
                let d = [0, 1, 2].map(|a| positions[dof][a] as f64 - mean[a]);
                v[dof] = match motion {
                    0..3 => f64::from(u8::from(c == motion)),
                    // e_k × d, component c.
                    _ => {
                        let k = motion - 3;
                        let (k1, k2) = ((k + 1) % 3, (k + 2) % 3);
                        if c == k1 {
                            -d[k2]
                        } else if c == k2 {
                            d[k1]
                        } else {
                            0.0
                        }
                    }
                };
            }
            out.push(v);
        }
    }
    out
}

/// The rigid motions as locked modes: M-orthonormalised (Gram–Schmidt,
/// twice), each with λ and the Lanczos-metric residual from one solve
/// with the shifted operator, θ = xᵀM y and ‖y − θx‖_M / θ for
/// y = (K − σM)⁻¹ M x.
fn rigid_modes(
    m: &Csr,
    hierarchy: &Hierarchy,
    (sigma, max_iters): (f64, usize),
    vectors: Vec<Vec<f64>>,
) -> Result<Vec<Locked>, AnalysisError> {
    let n = m.rows();
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let mul_m = |v: &[f64]| -> Result<Vec<f64>, AnalysisError> {
        let mut out = vec![0.0; n];
        m.mul_vec(v, &mut out).map_err(|e| internal("modal", e))?;
        Ok(out)
    };
    let mut basis: Vec<(Vec<f64>, Vec<f64>)> = Vec::new();
    for mut v in vectors {
        for _ in 0..2 {
            for (x, mx) in &basis {
                let c = dot(mx, &v);
                v.iter_mut().zip(x).for_each(|(a, b)| *a -= c * b);
            }
        }
        let mv = mul_m(&v)?;
        let norm = dot(&v, &mv).max(0.0).sqrt();
        if norm > 0.0 {
            basis.push((
                v.iter().map(|a| a / norm).collect(),
                mv.iter().map(|a| a / norm).collect(),
            ));
        }
    }
    let mut out = Vec::with_capacity(basis.len());
    for (x, mx) in basis {
        let solution = pcg(hierarchy, &mx, 1e-12, max_iters).map_err(|e| internal("pcg", e))?;
        if !solution.converged {
            return Err(AnalysisError::SolverNotConverged {
                history: solution.history,
            });
        }
        let y = solution.x;
        let theta = dot(&mx, &y);
        let r: Vec<f64> = y.iter().zip(&x).map(|(a, b)| a - theta * b).collect();
        let residual = dot(&r, &mul_m(&r)?).max(0.0).sqrt() / theta.abs().max(f64::MIN_POSITIVE);
        out.push((sigma + 1.0 / theta, residual, x, mx));
    }
    Ok(out)
}

/// Mode `x` (free dofs) at the mesh vertices, per solid.
fn vertex_shapes(
    model: &AnalysisModel,
    prepared: &super::solve::Prepared,
    x: &[f64],
) -> Result<Vec<Vec<[f64; 3]>>, AnalysisError> {
    let mut out = Vec::new();
    for (s, (_, _, moments, aggregation)) in prepared.solids.iter().enumerate() {
        let e = &aggregation.extension;
        let ev = vector_extension(e).map_err(|e| internal("extension", e))?;
        let mut u = vec![0.0; ev.rows()];
        let start = prepared.offsets.get(s).copied().unwrap_or(0);
        ev.mul_vec(&x[start..start + ev.cols()], &mut u)
            .map_err(|e| internal("extension", e))?;
        let components: [Vec<f64>; 3] =
            std::array::from_fn(|c| u.iter().skip(c).step_by(3).copied().collect());
        let mut shape = Vec::with_capacity(model.solids[s].vertices.len());
        for v in &model.solids[s].vertices {
            let mut d = [0.0; 3];
            for (slot, values) in d.iter_mut().zip(&components) {
                *slot = point_value(&prepared.grid, moments, e, values, *v)?.unwrap_or(0.0);
            }
            shape.push(d);
        }
        out.push(shape);
    }
    Ok(out)
}

#[cfg(test)]
#[path = "modal_tests.rs"]
mod tests;
