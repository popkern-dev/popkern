//! Linear elastic statics (map §1.2–1.5): condensed EAS cells, Nitsche
//! supports, traction loads, AgFEM, multigrid PCG, a rigid-mode check
//! before solving, and the force balance (X6) after.

pub mod assemble;
mod field;
pub mod from_setup;
pub mod generative;
pub mod infill;
mod lanczos;
pub mod modal;
pub mod modal_restarts;
pub mod rigid;
mod solve;
mod stability;
pub mod submodel;
pub mod support;
mod system;
pub mod thermoelastic;

use crate::error::AnalysisError;
use crate::grid::Grid;
use crate::model::AnalysisModel;
use crate::moments::{SolidMoments, moment_index};
use crate::options::AnalysisOptions;
use crate::report::SolveReport;
use crate::thermal::sample::point_value;

use assemble::FaceMechanics;

/// Largest Poisson's ratio accepted. Measured (PH-034): with the
/// Nitsche flux taken from the enhanced stress, a ν = 0.49 cantilever
/// stays definite at 1–3 cells aligned and rotated, and is within 2.8%
/// (aligned) and 6.5% (rotated 17°) of the 8-cell solution at 2 cells;
/// with the plain Q1 flux it was indefinite. Above 0.49 is unmeasured.
pub const POISSON_LIMIT: f64 = 0.49;

/// Eigenvalues of the rigid-mode matrix below this fraction of its
/// trace are free motions.
pub const RIGID_MODE_TOL: f64 = 1.0e-10;

/// An elastic problem: material (MPa) and a condition per face.
#[derive(Debug, Clone, PartialEq)]
pub struct ElasticProblem {
    /// Young's modulus, MPa.
    pub youngs: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Condition of each face, indexed by `FaceId`.
    pub faces: Vec<FaceMechanics>,
}

/// The result of a static solve.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuralResult {
    /// Displacement at every mesh vertex, per solid, mm.
    pub vertex_displacement: Vec<Vec<[f64; 3]>>,
    /// Largest displacement magnitude over the mesh vertices and the
    /// cell-resolution surface samples, mm.
    pub max_displacement: f64,
    /// Force each face passes to its supports, N, by `FaceId`.
    pub face_reaction: Vec<[f64; 3]>,
    /// The displacement field.
    pub field: DisplacementField,
    /// Recovered stresses at the mesh vertices, per solid.
    pub stresses: Vec<crate::recover::SolidStresses>,
    /// Global relative ZZ error over cells with an estimate.
    pub global_error: Option<f64>,
    /// What the solve did (incl. the force balance).
    pub report: SolveReport,
}

/// The solved displacement field, evaluable anywhere in the part.
#[derive(Debug, Clone, PartialEq)]
pub struct DisplacementField {
    grid: Grid,
    /// Per solid: moments, extension, and nodal values per component.
    solids: Vec<(SolidMoments, crate::aggregate::Extension, [Vec<f64>; 3])>,
}

impl DisplacementField {
    /// Displacement at `p`, mm; `None` outside every solid.
    pub fn displacement_at(&self, p: [f64; 3]) -> Result<Option<[f64; 3]>, AnalysisError> {
        for (moments, extension, values) in &self.solids {
            let mut out = [0.0; 3];
            let mut found = true;
            for (slot, v) in out.iter_mut().zip(values) {
                match point_value(&self.grid, moments, extension, v, p)? {
                    Some(x) => *slot = x,
                    None => found = false,
                }
            }
            if found {
                return Ok(Some(out));
            }
        }
        Ok(None)
    }
}

fn internal(stage: &'static str, e: impl std::fmt::Display) -> AnalysisError {
    AnalysisError::Internal {
        stage,
        detail: e.to_string(),
    }
}

/// Combined global ZZ error over solids (each solid's own estimate is
/// already relative; the worst is reported).
fn global_error(stresses: &[crate::recover::SolidStresses]) -> Option<f64> {
    stresses
        .iter()
        .filter_map(|s| s.global_error)
        .fold(None, |acc: Option<f64>, e| {
            Some(acc.map_or(e, |a| a.max(e)))
        })
}

/// Centroid of a solid from its cell moments.
fn centroid(grid: &Grid, m: &SolidMoments) -> [f64; 3] {
    let h = grid.cell_size();
    let mut first = [0.0; 3];
    for c in &m.cells {
        let centre = grid
            .cell_min(grid.cell_ijk(c.cell).unwrap_or_default())
            .map(|x| x + 0.5 * h);
        let e = [
            moment_index(1, 0, 0),
            moment_index(0, 1, 0),
            moment_index(0, 0, 1),
        ];
        for a in 0..3 {
            first[a] += centre[a] * c.m[0] + 0.5 * h * c.m[e[a]];
        }
    }
    first.map(|x| x / m.volume)
}

/// Solves `problem` on `model` with `options`.
pub fn solve_elastic(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    options: &AnalysisOptions,
) -> Result<StructuralResult, AnalysisError> {
    solve_elastic_with(
        model,
        problem,
        options,
        crate::recover::Fallback::InCell,
        None,
        crate::moments::Splitting::Exact,
    )
}

/// Solves `problem` on `model` with the thermal strain of `thermal` (a
/// staggered thermoelastic solve: the temperature comes first).
pub fn solve_thermoelastic(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    thermal: &thermoelastic::ThermalLoad<'_>,
    options: &AnalysisOptions,
) -> Result<StructuralResult, AnalysisError> {
    solve::solve_core(
        model,
        problem,
        options,
        crate::recover::Fallback::InCell,
        None,
        crate::moments::Splitting::Exact,
        &solve::Extras {
            thermal: Some(thermal),
            ..solve::Extras::default()
        },
    )
}

/// Solves `problem` on `model` with an optional thermal strain and an
/// optional infill (both, either or neither).
pub fn solve_with(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    thermal: Option<&thermoelastic::ThermalLoad<'_>>,
    infill: Option<&infill::Infill>,
    options: &AnalysisOptions,
) -> Result<StructuralResult, AnalysisError> {
    solve::solve_core(
        model,
        problem,
        options,
        crate::recover::Fallback::InCell,
        None,
        crate::moments::Splitting::Exact,
        &solve::Extras {
            thermal,
            infill,
            ..solve::Extras::default()
        },
    )
}

/// Solves `problem` on `model` with an infill: solid walls around a
/// homogenised gyroid core (`infill`).
pub fn solve_infilled(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    infill: &infill::Infill,
    options: &AnalysisOptions,
) -> Result<StructuralResult, AnalysisError> {
    solve::solve_core(
        model,
        problem,
        options,
        crate::recover::Fallback::InCell,
        None,
        crate::moments::Splitting::Exact,
        &solve::Extras {
            infill: Some(infill),
            ..solve::Extras::default()
        },
    )
}

/// [`solve_elastic`] with an explicit recovery fallback (tests use the
/// constant one as a negative control) and fragment splitting mode
/// (welded and lattice moments are CC1's and LF1's controls).
pub(crate) fn solve_elastic_with(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    options: &AnalysisOptions,
    fallback: crate::recover::Fallback,
    exact: Option<crate::recover::ExactStress<'_>>,
    split: crate::moments::Splitting,
) -> Result<StructuralResult, AnalysisError> {
    solve::solve_core(
        model,
        problem,
        options,
        fallback,
        exact,
        split,
        &solve::Extras::default(),
    )
}

#[cfg(test)]
#[path = "bore_tests.rs"]
mod bore_tests;
#[cfg(test)]
#[path = "nitsche_tests.rs"]
mod nitsche_tests;
#[cfg(test)]
#[path = "octant_tests.rs"]
mod octant_tests;
#[cfg(test)]
#[path = "recovery_tests.rs"]
mod recovery_tests;
#[cfg(test)]
#[path = "refine_tests.rs"]
mod refine_tests;
#[cfg(test)]
#[path = "roller_tests.rs"]
mod roller_tests;
#[cfg(test)]
#[path = "slit_tests.rs"]
mod slit_tests;
#[cfg(test)]
#[path = "structural_tests.rs"]
mod structural_tests;
#[cfg(test)]
#[path = "support_tests.rs"]
mod support_tests;
#[cfg(test)]
#[path = "thin_wall_tests.rs"]
mod thin_wall_tests;
#[cfg(test)]
#[path = "two_gap_tests.rs"]
mod two_gap_tests;
#[cfg(test)]
#[path = "vector_mg_tests.rs"]
mod vector_mg_tests;
#[cfg(test)]
#[path = "wall_split_tests.rs"]
mod wall_split_tests;
#[cfg(test)]
#[path = "wrap_tests.rs"]
mod wrap_tests;
