//! Steady heat conduction with fixed temperatures (Nitsche),
//! convection, heat flux and volumetric sources (map Stage 1b).
//!
//! Geometry is in mm; conductivity and coefficients are given in SI and
//! converted. Fixed temperatures are imposed weakly by Nitsche, which
//! is consistent, so linear temperature fields are reproduced exactly.

pub mod assemble;
pub mod from_setup;
mod radiation;
pub mod sample;
mod stability;
pub mod submodel;
mod system;
pub mod transient;

use crate::error::AnalysisError;
use crate::grid::Grid;
use crate::linalg::sparse::Csr;
use crate::model::{AnalysisModel, FaceId, SolidId};
use crate::options::AnalysisOptions;
use crate::report::SolveReport;
use crate::solve::pcg::pcg;

use assemble::MmCondition;

/// Relative tolerance of the X5 volume check on every solve.
pub const VOLUME_PARITY_TOL: f64 = 1.0e-9;

/// A boundary condition on one face, in user units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FaceCondition {
    /// Insulated (no condition).
    Insulated,
    /// T = t0 + gradient·x, °C with x in mm.
    Fixed {
        /// Value at the origin, °C.
        t0_c: f64,
        /// Gradient, °C/mm.
        gradient: [f64; 3],
    },
    /// Convection, W/(m²·K), to an ambient in °C.
    Convection {
        /// Effective coefficient, W/(m²·K).
        h: f64,
        /// Ambient, °C.
        ambient_c: f64,
    },
    /// Uniform heat flux density into the part, W/mm².
    Flux {
        /// Flux density, W/mm².
        w_per_mm2: f64,
    },
    /// Held (Nitsche) to a parent solution's temperature, affine at each
    /// piece from the field's value and gradient at the piece's
    /// centroid: a submodel's cut faces, driven by the coarse field.
    Prescribed,
    /// Convection plus grey radiation to surroundings at the ambient:
    /// heat out h (T − T∞) + εσ(T⁴ − T∞⁴), solved by fixed point on the
    /// radiation coefficient (`radiation`).
    Radiating {
        /// Convective coefficient alone, W/(m²·K).
        h: f64,
        /// Surface emittance ε.
        emittance: f64,
        /// Ambient (and surroundings), °C.
        ambient_c: f64,
    },
}

/// A conduction problem on a model.
#[derive(Debug, Clone, PartialEq)]
pub struct ConductionProblem {
    /// Conductivity, W/(m·K).
    pub conductivity: f64,
    /// Condition of each face, indexed by `FaceId`.
    pub faces: Vec<FaceCondition>,
    /// Total volumetric power, W, spread uniformly over all solids.
    pub source_w: f64,
    /// Power generated uniformly in single solids, W.
    pub body_sources: Vec<(SolidId, f64)>,
    /// Heat deflection temperature, °C, for the service warning.
    pub hdt_c: Option<f64>,
}

/// Heat flows of a solution, W. Face flows are split into gross in and
/// out, so heat passing through a part is visible even when the net is
/// zero.
#[derive(Debug, Clone, PartialEq)]
pub struct HeatBalance {
    /// Imposed flux into the part.
    pub flux_in: f64,
    /// Volumetric generation.
    pub source: f64,
    /// Heat leaving through each face (negative enters), by `FaceId`;
    /// fixed and convective faces.
    pub face_out: Vec<f64>,
}

impl HeatBalance {
    /// Gross heat entering through fixed or convective faces.
    pub fn faces_in(&self) -> f64 {
        self.face_out.iter().filter(|q| **q < 0.0).map(|q| -q).sum()
    }

    /// Gross heat leaving through fixed or convective faces.
    pub fn faces_out(&self) -> f64 {
        self.face_out.iter().filter(|q| **q > 0.0).sum()
    }

    /// In − out; zero to solver tolerance.
    pub fn imbalance(&self) -> f64 {
        self.flux_in + self.source + self.faces_in() - self.faces_out()
    }

    /// The largest of the gross flows: the scale of the imbalance.
    pub fn throughput(&self) -> f64 {
        (self.flux_in.abs() + self.source.abs() + self.faces_in()).max(self.faces_out())
    }
}

/// The solved temperature field, evaluable anywhere in the part.
#[derive(Debug, Clone, PartialEq)]
pub struct TemperatureField {
    grid: Grid,
    solids: Vec<(
        crate::moments::SolidMoments,
        crate::aggregate::Extension,
        Vec<f64>,
    )>,
}

impl TemperatureField {
    /// Temperature at `p` (°C), from the first solid with an active
    /// cell containing it; `None` outside every solid's active cells.
    pub fn temperature_at(&self, p: [f64; 3]) -> Result<Option<f64>, AnalysisError> {
        Ok(self.value_and_gradient_at(p)?.map(|(t, _)| t))
    }

    /// Temperature (°C) and its gradient (°C/mm) at `p`, as
    /// [`TemperatureField::temperature_at`].
    pub fn value_and_gradient_at(
        &self,
        p: [f64; 3],
    ) -> Result<Option<(f64, [f64; 3])>, AnalysisError> {
        for (moments, extension, values) in &self.solids {
            if let Some(t) = sample::point_eval(&self.grid, moments, extension, values, p)? {
                return Ok(Some(t));
            }
        }
        Ok(None)
    }
}

/// The result of a steady thermal solve.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalResult {
    /// Temperature at every mesh vertex, per solid, °C.
    pub vertex_temperature: Vec<Vec<f64>>,
    /// Highest vertex temperature, °C.
    pub max_temperature_c: f64,
    /// Heat flux q = −k∇T at every mesh vertex, per solid, W/m², from
    /// the SPR-recovered gradient (as stresses are recovered).
    pub vertex_heat_flux: Vec<Vec<[f64; 3]>>,
    /// Largest vertex heat flux magnitude, W/m².
    pub max_heat_flux_w_per_m2: f64,
    /// Heat flows.
    pub heat: HeatBalance,
    /// The field, for evaluation away from the vertices.
    pub field: TemperatureField,
    /// What the solve did.
    pub report: SolveReport,
    /// A transient's steps (empty for a steady solve); the rest of the
    /// result is its final state.
    pub transient: Vec<transient::TransientStep>,
}

/// Solves `problem` on `model` with `options`.
pub fn solve_conduction(
    model: &AnalysisModel,
    problem: &ConductionProblem,
    options: &AnalysisOptions,
) -> Result<ThermalResult, AnalysisError> {
    solve_conduction_core(model, problem, options, None)
}

/// [`solve_conduction`], with the field `Prescribed` faces are held to.
pub(crate) fn solve_conduction_core(
    model: &AnalysisModel,
    problem: &ConductionProblem,
    options: &AnalysisOptions,
    prescribed: Option<&TemperatureField>,
) -> Result<ThermalResult, AnalysisError> {
    radiation::solve_radiating(model, problem, options, prescribed)
}

/// One linear solve: `Radiating` faces take their radiation coefficient
/// from `previous` (the ambient where there is none).
fn solve_once(
    model: &AnalysisModel,
    problem: &ConductionProblem,
    options: &AnalysisOptions,
    prescribed: Option<&TemperatureField>,
    previous: Option<&TemperatureField>,
) -> Result<ThermalResult, AnalysisError> {
    let mut prepared = system::prepare(model, problem, options, prescribed, previous)?;
    let x = solve_prepared(&mut prepared, options)?;
    system::finish(model, problem, prepared, &x, (prescribed, previous))
}

/// Solves `prepared`'s assembled system (its operator moves into the
/// hierarchy), recording the iterations and history in its report.
fn solve_prepared(
    prepared: &mut system::Prepared,
    options: &AnalysisOptions,
) -> Result<Vec<f64>, AnalysisError> {
    let empty = Csr::from_triplets(0, 0, &[]).map_err(|e| internal("assemble", e))?;
    let a = std::mem::replace(&mut prepared.a, empty);
    let positions = std::mem::take(&mut prepared.positions);
    let hierarchy = system::hierarchy(a, positions, options, &mut prepared.report)?;
    let report = &mut prepared.report;
    let solution = pcg(
        &hierarchy,
        &prepared.b,
        options.cg_rel_tol,
        options.cg_max_iters,
    )
    .map_err(|e| internal("pcg", e))?;
    report.cg_iterations = solution.iterations;
    report.residual_history = solution.history.clone();
    if !solution.converged {
        return Err(AnalysisError::SolverNotConverged {
            history: solution.history,
        });
    }
    Ok(solution.x)
}

fn internal(stage: &'static str, e: impl std::fmt::Display) -> AnalysisError {
    AnalysisError::Internal {
        stage,
        detail: e.to_string(),
    }
}

fn grid_failure(e: crate::grid::GridError, options: &AnalysisOptions) -> AnalysisError {
    use crate::grid::GridError;
    match e {
        GridError::ZeroExtent { .. } | GridError::NonFiniteBounds => {
            AnalysisError::DegenerateGeometry {
                solid: None,
                detail: crate::error::Degeneracy::ZeroExtent,
            }
        }
        GridError::TooLarge { dims } => AnalysisError::CellBudgetExceeded {
            requested: dims.iter().fold(1usize, |a, &d| a.saturating_mul(d)),
            cap: options.max_active_cells,
        },
        other => internal("grid", other),
    }
}

fn to_mm(c: FaceCondition) -> MmCondition {
    match c {
        FaceCondition::Insulated => MmCondition::Insulated,
        FaceCondition::Fixed { t0_c, gradient } => MmCondition::Fixed { t0: t0_c, gradient },
        FaceCondition::Convection { h, ambient_c } => MmCondition::Convection {
            h: 1e-6 * h,
            ambient: ambient_c,
        },
        FaceCondition::Flux { w_per_mm2 } => MmCondition::Flux { q: w_per_mm2 },
        FaceCondition::Prescribed => MmCondition::Prescribed,
        FaceCondition::Radiating {
            h,
            emittance,
            ambient_c,
        } => MmCondition::Radiating {
            h: 1e-6 * h,
            emittance,
            ambient: ambient_c,
        },
    }
}

fn bounds(model: &AnalysisModel) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for v in model.solids.iter().flat_map(|s| &s.vertices) {
        for a in 0..3 {
            lo[a] = lo[a].min(v[a]);
            hi[a] = hi[a].max(v[a]);
        }
    }
    (lo, hi)
}

/// Sum of the areas of `faces` in `model`, mm².
pub fn faces_area(model: &AnalysisModel, faces: &[FaceId]) -> f64 {
    faces
        .iter()
        .filter_map(|f| model.faces.get(f.0 as usize))
        .map(|f| f.area)
        .sum()
}

#[cfg(test)]
#[path = "source_tests.rs"]
mod source_tests;
#[cfg(test)]
#[path = "thermal_tests.rs"]
mod thermal_tests;
