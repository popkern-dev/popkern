//! Generative design: minimum-compliance density topology optimisation
//! on the analysis grid (docs/architecture/analysis-generative-design.md;
//! Andreassen 2011 top88, Liu & Tovar 2014 top3d).
//!
//! Every active entry is a design variable x ∈ [0, 1]. The density
//! filter gives the physical densities x̃, and modified SIMP scales each
//! entry's condensed stiffness, and its Nitsche support terms with it,
//! by ε + x̃ᵖ(1 − ε): a support holds only the material that is there,
//! as a nodal support does in top88. Applied loads are not scaled.
//! Entries near faces the caller names are kept solid. The optimality
//! criteria update moves the design to the volume target; the loop
//! stops when no variable moves by more than the tolerance.

mod blocks;
pub mod body;
mod cells;
pub mod filter;
mod keep;
pub mod mma;
mod optimise;
mod options;
pub mod print;
pub mod projection;
mod scaled;
pub mod surface;
pub mod thickness;
pub mod update;

use popkern_math::PortableMath;

use crate::error::AnalysisError;
use crate::model::AnalysisModel;
use crate::options::AnalysisOptions;
use crate::report::{SolveReport, Warning};

use super::assemble::{FaceMechanics, vector_extension};
use super::solve::{Extras, Prepared, prepare, solve_system_reusing};
use super::{ElasticProblem, internal};

use blocks::{Block, Blocks, sensitivity_blocks};
pub use filter::{DensityFilter, DesignEntry, design_entries};
use keep::keep_floors;
pub use optimise::{Design, Iteration, Progress, optimise, optimise_observed};
pub use options::{GenerativeOptions, UpdateMethod};

/// The share of a part's volume thinner than the filter's smallest
/// member above which the design is warned of (`Warning::ThinDesignSpace`):
/// a block's corners and edges are 1–5% (12 × 4 × 2 at 0.5 and 0.25 mm
/// cells, r = 0.75 mm), a plate thinner than 2r all of it.
const THIN_SHARE: f64 = 0.2;

/// Cells per filter radius below which the design is warned of as not
/// converged in the grid (`Warning::FilterUnderResolved`).
pub const FILTER_CELLS: f64 = 3.0;

/// Compliance and its sensitivities for one design.
#[derive(Debug, Clone, PartialEq)]
pub struct Evaluation {
    /// Compliance fᵀu, N·mm (the eroded design's under the robust
    /// formulation: the largest of the three, as stiffness grows with
    /// density and eroded ≤ intermediate ≤ dilated entry by entry).
    pub compliance: f64,
    /// Physical densities x̄ (the intermediate design's under the
    /// robust formulation).
    pub physical: Vec<f64>,
    /// ∂c/∂x (through the filter; 0 on kept entries).
    pub gradient: Vec<f64>,
    /// ∂V/∂x, mm³ (likewise).
    pub volume_gradient: Vec<f64>,
    /// PCG iterations of the solve.
    pub cg_iterations: usize,
}

/// A design's solve, kept across evaluations: the geometry, the filter
/// and the kept entries are built once.
pub struct Evaluator {
    options: AnalysisOptions,
    /// CG's relative residual for the next evaluations (at least the
    /// analysis tolerance).
    tolerance: f64,
    design: GenerativeOptions,
    prepared: Prepared,
    positions: Vec<[usize; 3]>,
    labels: Vec<u8>,
    entries: Vec<DesignEntry>,
    filter: DensityFilter,
    floors: Vec<f64>,
    blocks: Blocks,
    scaled: scaled::ScaledSystem,
    /// The last solve's hierarchy, whose prolongations the next reuses.
    hierarchy: Option<crate::solve::Hierarchy>,
    projection: Option<projection::Projection>,
    robust: Option<f64>,
    print: Option<print::PrintFilter>,
}

/// A design through one projection and the print filter.
struct Stage {
    projected: Vec<f64>,
    printed: Option<print::Printed>,
    physical: Vec<f64>,
    held: Vec<bool>,
    slope: Vec<f64>,
}

impl Evaluator {
    /// The evaluator of `problem` on `model`, holding solid the material
    /// near the faces in `keep`. Supports must hold at zero
    /// displacement: a prescribed displacement's load would move with
    /// the densities.
    pub fn new(
        model: &AnalysisModel,
        problem: &ElasticProblem,
        design: GenerativeOptions,
        keep: &[crate::model::FaceId],
        options: &AnalysisOptions,
    ) -> Result<Evaluator, AnalysisError> {
        design.validate()?;
        let moving = problem.faces.iter().any(|f| match f {
            FaceMechanics::Support(crate::element::vector_boundary::PieceSupport::Fixed {
                u0,
                gradient,
            }) => u0
                .iter()
                .chain(gradient.iter().flatten())
                .any(|v| *v != 0.0),
            FaceMechanics::Prescribed | FaceMechanics::PrescribedTraction => true,
            _ => false,
        });
        if moving {
            return Err(AnalysisError::Unsupported {
                what: "generative design with a support at a non-zero displacement",
            });
        }
        let mut prepared = prepare(
            model,
            problem,
            options,
            (crate::moments::Splitting::Exact, false),
            &Extras::default(),
        )?;
        let positions = std::mem::take(&mut prepared.positions);
        let labels = std::mem::take(&mut prepared.labels);
        let moments: Vec<&crate::moments::SolidMoments> =
            prepared.solids.iter().map(|(_, _, m, _)| m).collect();
        let entries = design_entries(&prepared.grid, &moments);
        let filter = DensityFilter::new(&prepared.grid, &moments, &entries, design.filter_radius);
        let floors = keep_floors(&prepared, keep, &entries, design.keep_thickness);
        let cell = prepared.grid.cell_size();
        if design.filter_radius < FILTER_CELLS * cell {
            prepared.report.warnings.push(Warning::FilterUnderResolved {
                radius: design.filter_radius,
                cell,
            });
        }
        let ones = vec![1.0; entries.len()];
        let member = 2.0 * design.filter_radius;
        if let Some(thin) = cells::CellField::new(&entries, &ones, prepared.grid.cell_size())
            .and_then(|c| thickness::thin_region(&c, member))
            .filter(|t| t.fraction > THIN_SHARE)
        {
            prepared.report.warnings.push(Warning::ThinDesignSpace {
                fraction: thin.fraction,
                thickness: thin.thickness,
                member,
            });
        }
        let blocks = sensitivity_blocks(&prepared, problem)?;
        let scaled = scaled::ScaledSystem::new(&prepared, problem)?;
        let print = design
            .build_direction
            .map(|d| print::PrintFilter::new(&entries, prepared.grid.cell_size(), d));
        Ok(Evaluator {
            options: *options,
            tolerance: options.cg_rel_tol,
            design,
            prepared,
            positions,
            labels,
            entries,
            filter,
            floors,
            blocks,
            scaled,
            hierarchy: None,
            projection: None,
            robust: None,
            print,
        })
    }

    /// The design's entries.
    pub fn entries(&self) -> &[DesignEntry] {
        &self.entries
    }

    /// Each entry's density floor (the kept material).
    pub fn floors(&self) -> &[f64] {
        &self.floors
    }

    /// The density filter.
    pub fn filter(&self) -> &DensityFilter {
        &self.filter
    }

    /// The print filter, when on.
    pub fn print(&self) -> Option<&print::PrintFilter> {
        self.print.as_ref()
    }

    /// The solve's report (cell size, budget, warnings).
    pub fn report(&self) -> &SolveReport {
        &self.prepared.report
    }

    /// The projection the next evaluations use (`None`: off).
    pub fn set_projection(&mut self, projection: Option<projection::Projection>) {
        self.projection = projection;
    }

    /// CG's relative residual for the next evaluations (the analysis
    /// tolerance when looser than it).
    pub fn set_tolerance(&mut self, tolerance: f64) {
        self.tolerance = tolerance;
    }

    /// The robust formulation's offset threshold η the next evaluations
    /// use (`None`: off; with the projection on, the compliance is the
    /// eroded design's).
    pub fn set_robust(&mut self, eta: Option<f64>) {
        self.robust = eta;
    }

    /// The projection at threshold `threshold(η)` under the robust
    /// formulation, else the design's own.
    fn shifted(&self, threshold: impl Fn(f64) -> f64) -> Option<projection::Projection> {
        match (self.projection, self.robust) {
            (Some(pr), Some(eta)) => Some(projection::Projection {
                beta: pr.beta,
                kind: projection::ProjectionKind::Threshold(threshold(eta)),
            }),
            (pr, _) => pr,
        }
    }

    /// Design `filtered` through `projection` and the print filter, and
    /// its floors.
    fn stage(&self, filtered: &[f64], projection: Option<projection::Projection>) -> Stage {
        let projected: Vec<f64> = match projection {
            Some(pr) => filtered.iter().map(|&t| pr.apply(t)).collect(),
            None => filtered.to_vec(),
        };
        // The printed design when the print filter is on; where the floor
        // holds the density, the design does not move it.
        let printed = self
            .print
            .as_ref()
            .map(|f| f.apply(&projected, &self.floors));
        let (physical, held): (Vec<f64>, Vec<bool>) = match &printed {
            Some(p) => (p.densities.clone(), p.held().to_vec()),
            None => projected
                .iter()
                .zip(&self.floors)
                .map(|(r, f)| (r.max(*f), r < f))
                .unzip(),
        };
        // dx̄/dx̃ for the chain rule.
        let slope: Vec<f64> = match projection {
            Some(pr) => filtered.iter().map(|&t| pr.slope(t)).collect(),
            None => vec![1.0; filtered.len()],
        };
        Stage {
            projected,
            printed,
            physical,
            held,
            slope,
        }
    }

    /// ∂/∂x of a function of `stage`'s physical densities, from its
    /// derivative `d` per entry (0 where the floor holds).
    fn back(&self, stage: &Stage, d: &[f64]) -> Vec<f64> {
        let d: Vec<f64> = d
            .iter()
            .zip(&stage.held)
            .map(|(g, &k)| if k { 0.0 } else { *g })
            .collect();
        let d = match (&self.print, &stage.printed) {
            (Some(f), Some(p)) => f.back(&stage.projected, p, &d),
            _ => d,
        };
        let d: Vec<f64> = d.iter().zip(&stage.slope).map(|(g, s)| g * s).collect();
        self.filter.back(&d)
    }

    /// Compliance and sensitivities of design `x` (one per entry).
    pub fn evaluate(&mut self, x: &[f64]) -> Result<Evaluation, AnalysisError> {
        let (p, eps) = (self.design.penalty, self.design.void_stiffness);
        let filtered = self.filter.apply(x);
        // The design itself, and the one solved (the eroded under the
        // robust formulation).
        let design = self.stage(&filtered, self.projection);
        let solved = self
            .robust
            .map(|_| self.stage(&filtered, self.shifted(|eta| 1.0 - eta)));
        let solved = solved.as_ref().unwrap_or(&design);
        let mut scales: Vec<Vec<f64>> = self
            .prepared
            .solids
            .iter()
            .map(|(_, _, m, _)| vec![1.0; m.cells.len()])
            .collect();
        for (d, rho) in self.entries.iter().zip(&solved.physical) {
            scales[d.solid][d.entry] = eps + rho.ppowf(p) * (1.0 - eps);
        }
        let (blocks, b) = self.scaled.assemble(&self.blocks, &scales)?;
        let (solution, hierarchy) = solve_system_reusing(
            blocks,
            (self.positions.clone(), self.labels.clone()),
            &b,
            &AnalysisOptions {
                cg_rel_tol: self.tolerance.max(self.options.cg_rel_tol),
                ..self.options
            },
            &mut self.prepared.report,
            self.hierarchy.as_ref(),
        )?;
        self.hierarchy = Some(hierarchy);
        let compliance: f64 = b.iter().zip(&solution.x).map(|(f, u)| f * u).sum();
        let energies = self.energies(&solution.x)?;
        let dc: Vec<f64> = solved
            .physical
            .iter()
            .zip(&energies)
            .map(|(rho, energy)| -p * rho.ppowf(p - 1.0) * (1.0 - eps) * energy)
            .collect();
        let volumes: Vec<f64> = self.entries.iter().map(|d| d.volume).collect();
        let gradient = self.back(solved, &dc);
        let volume_gradient = self.back(&design, &volumes);
        // An infinite sensitivity (the print filter's P-norm once
        // overflowed) would reach the update as NaN.
        if gradient
            .iter()
            .chain(&volume_gradient)
            .any(|g| !g.is_finite())
        {
            return Err(internal("generative", "non-finite sensitivity"));
        }
        Ok(Evaluation {
            compliance,
            physical: design.physical,
            gradient,
            volume_gradient,
            cg_iterations: solution.iterations,
        })
    }

    /// u_eᵀ k_e u_e per entry for the unscaled blocks (material and
    /// supports), from the free-dof solution `x`.
    fn energies(&self, x: &[f64]) -> Result<Vec<f64>, AnalysisError> {
        let prepared = &self.prepared;
        let grid = &prepared.grid;
        let mut out = Vec::with_capacity(self.entries.len());
        for (s, (_, _, moments, aggregation)) in prepared.solids.iter().enumerate() {
            let e = &aggregation.extension;
            let ev = vector_extension(e).map_err(|e| internal("generative", e))?;
            let offset = prepared.offsets[s];
            let mut u = vec![0.0; ev.rows()];
            ev.mul_vec(&x[offset..offset + ev.cols()], &mut u)
                .map_err(|e| internal("generative", e))?;
            let own = &self.blocks.own[s];
            let energy = |entry: usize| -> Result<f64, AnalysisError> {
                let k: &Block = own
                    .get(entry)
                    .and_then(|b| b.as_deref())
                    .unwrap_or(&self.blocks.full);
                let nodes = moments
                    .node_keys(grid, entry)
                    .ok_or_else(|| internal("generative", "cell outside grid"))?;
                let mut ue = [0.0; 24];
                for (a, n) in nodes.iter().enumerate() {
                    let at = e
                        .position(*n)
                        .ok_or_else(|| internal("generative", "node without a value"))?;
                    for i in 0..3 {
                        ue[3 * a + i] = u[3 * at + i];
                    }
                }
                Ok((0..24)
                    .map(|r| ue[r] * (0..24).map(|c| k[r][c] * ue[c]).sum::<f64>())
                    .sum())
            };
            for value in crate::parallel::map_range(0..moments.cells.len(), 256, energy) {
                out.push(value?);
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
#[path = "generative_tests.rs"]
mod tests;
