//! Typed declines and failures of an analysis.
//!
//! Every way an analysis can refuse or fail is a variant here, so the
//! caller can explain it and tests can assert on it. Nothing is
//! reported by printing.

use crate::model::{FaceId, ModelError, SolidId};
use crate::setup::SetupError;

/// Why the analysis declined or failed.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum AnalysisError {
    /// The geometry model failed validation.
    #[error("invalid analysis model: {0}")]
    InvalidModel(#[from] ModelError),

    /// The setup failed validation.
    #[error("invalid setup: {0}")]
    InvalidSetup(#[from] SetupError),

    /// The thermal setup failed validation.
    #[error("invalid thermal setup: {0}")]
    InvalidThermalSetup(#[from] crate::thermal_setup::ThermalSetupError),

    /// A solid's mesh has open edges, so inside and outside are not
    /// defined.
    #[error("solid {solid:?} is not watertight: {open_edges} open edges")]
    NotWatertight {
        /// The solid with open edges.
        solid: SolidId,
        /// Number of open (unmatched) edges after welding.
        open_edges: usize,
    },

    /// No support is defined, so the structure is free to move.
    #[error("no support is defined")]
    NoSupport,

    /// Poisson's ratio is above the limit the Nitsche supports are
    /// measured stable and accurate for (PH-034: 0.49 with the enhanced
    /// flux; higher is unmeasured).
    #[error("Poisson's ratio {poisson} is above the supported limit {limit}")]
    NearlyIncompressible {
        /// The material's ν.
        poisson: f64,
        /// The largest ν accepted.
        limit: f64,
    },

    /// A solid can move freely under its supports.
    #[error("solid {solid:?} is free to move: {motions:?}")]
    FreeMotion {
        /// The solid.
        solid: SolidId,
        /// Its free rigid motions.
        motions: Vec<crate::structural::rigid::RigidMotion>,
    },

    /// A boundary condition names a face that no longer exists.
    #[error("boundary condition refers to face `{name}`, which no longer exists")]
    StaleBoundaryCondition {
        /// The face name in the setup.
        name: String,
    },

    /// A face reference matched no face.
    #[error("face reference `{name}` matched no face")]
    SelectorUnresolved {
        /// The face name in the setup.
        name: String,
    },

    /// A face reference matched a different number of faces than it
    /// expects. Never resolved by picking one.
    #[error("face reference `{name}` matched {matched} faces: {candidates:?}")]
    SelectorAmbiguous {
        /// The face name in the setup.
        name: String,
        /// Number of faces matched.
        matched: usize,
        /// The matched faces, in `FaceId` order.
        candidates: Vec<FaceId>,
    },

    /// A wall is thinner than the grid can resolve.
    #[error("wall of {min_wall} mm is under-resolved by cells of {cell} mm")]
    UnderResolved {
        /// Thinnest wall found, in model units.
        min_wall: f64,
        /// Cell size, in model units.
        cell: f64,
    },

    /// More separate parts of one solid meet in a cell (or around a
    /// node) than the dof encoding allows: features far thinner than
    /// the cell. A smaller cell size resolves them.
    #[error(
        "more than {limit} separate parts of the solid meet at grid cell or node {index} with cells of {cell} mm; use a smaller cell size"
    )]
    TooManyFragments {
        /// The cell or node index.
        index: usize,
        /// Most parts allowed.
        limit: usize,
        /// Cell size, in model units.
        cell: f64,
    },

    /// Multigrid coarsening stalled and a connected block of the
    /// coarsest level is too large to factor densely.
    #[error("coarsest multigrid level has a connected block of {dofs} dofs (cap {cap})")]
    CoarseLevelTooLarge {
        /// Dofs of the largest connected block of the coarsest level.
        dofs: usize,
        /// The cap.
        cap: usize,
    },

    /// The iterative solver did not reach its tolerance in its budget.
    #[error("solver did not converge in {} iterations", history.len())]
    SolverNotConverged {
        /// Relative residual after each iteration.
        history: Vec<f64>,
    },

    /// The requested resolution needs more cells than the budget.
    #[error("{requested} active cells requested, budget is {cap}")]
    CellBudgetExceeded {
        /// Cells the requested resolution needs.
        requested: usize,
        /// The `max_active_cells` budget.
        cap: usize,
    },

    /// The assembled system is singular: nothing fixes the solution
    /// level (no temperature reference, or no support).
    #[error("the system is singular: {detail}")]
    SingularSystem {
        /// What failed.
        detail: String,
    },

    /// The combination asked for is not supported (yet); declined
    /// rather than approximated.
    #[error("not supported: {what}")]
    Unsupported {
        /// What was asked for.
        what: &'static str,
    },

    /// A generative design's options are out of range.
    #[error("generative design: {what}")]
    InvalidGenerative {
        /// Which option, and its valid range.
        what: &'static str,
    },

    /// The material a generative design must keep (every entry with a
    /// supported or loaded piece) is already more than its volume
    /// target.
    #[error(
        "generative design: the kept material ({kept:.3} mm³) exceeds the volume target ({target:.3} mm³)"
    )]
    KeptVolumeExceedsTarget {
        /// Kept volume, mm³.
        kept: f64,
        /// Target volume, mm³.
        target: f64,
    },

    /// An internal invariant failed. Never expected; reported rather
    /// than hidden.
    #[error("internal error in {stage}: {detail}")]
    Internal {
        /// Pipeline stage.
        stage: &'static str,
        /// What failed.
        detail: String,
    },

    /// The geometry is degenerate in a way the analysis cannot use.
    #[error("degenerate geometry in solid {solid:?}: {detail:?}")]
    DegenerateGeometry {
        /// The solid concerned, when known.
        solid: Option<SolidId>,
        /// What is degenerate.
        detail: Degeneracy,
    },
}

impl AnalysisError {
    /// The variant's name, as reports and popkern live spell the kind of
    /// a decline (`"SelectorUnresolved"`, `"CellBudgetExceeded"`, ...).
    pub fn kind(&self) -> &'static str {
        match self {
            AnalysisError::InvalidModel(_) => "InvalidModel",
            AnalysisError::InvalidSetup(_) => "InvalidSetup",
            AnalysisError::InvalidThermalSetup(_) => "InvalidThermalSetup",
            AnalysisError::NotWatertight { .. } => "NotWatertight",
            AnalysisError::NoSupport => "NoSupport",
            AnalysisError::NearlyIncompressible { .. } => "NearlyIncompressible",
            AnalysisError::FreeMotion { .. } => "FreeMotion",
            AnalysisError::StaleBoundaryCondition { .. } => "StaleBoundaryCondition",
            AnalysisError::SelectorUnresolved { .. } => "SelectorUnresolved",
            AnalysisError::SelectorAmbiguous { .. } => "SelectorAmbiguous",
            AnalysisError::UnderResolved { .. } => "UnderResolved",
            AnalysisError::TooManyFragments { .. } => "TooManyFragments",
            AnalysisError::CoarseLevelTooLarge { .. } => "CoarseLevelTooLarge",
            AnalysisError::SolverNotConverged { .. } => "SolverNotConverged",
            AnalysisError::CellBudgetExceeded { .. } => "CellBudgetExceeded",
            AnalysisError::SingularSystem { .. } => "SingularSystem",
            AnalysisError::Unsupported { .. } => "Unsupported",
            AnalysisError::InvalidGenerative { .. } => "InvalidGenerative",
            AnalysisError::KeptVolumeExceedsTarget { .. } => "KeptVolumeExceedsTarget",
            AnalysisError::Internal { .. } => "Internal",
            AnalysisError::DegenerateGeometry { .. } => "DegenerateGeometry",
        }
    }
}

/// A kind of degenerate geometry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Degeneracy {
    /// The shell is oriented inward (negative enclosed volume).
    InsideOut,
    /// The solid encloses zero volume.
    ZeroVolume,
    /// The model's bounding box has zero extent on some axis.
    ZeroExtent,
    /// The rasterised volume disagrees with the mesh's enclosed volume
    /// beyond round-off: the mesh is not closed or the clipping failed.
    VolumeMismatch,
    /// A submodel's box cut failed: the box holds none of the solid, or
    /// a section did not close or triangulate.
    BoxCut,
}
