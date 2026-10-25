//! The solve report: everything a run measured, decided or fell back
//! on, as data.
//!
//! Investigation findings become fields here, never prints, so tests
//! assert on them and the JSON report carries them. Wall-clock timings
//! are not in this struct: they differ run to run, and the report must
//! be bit-identical for identical inputs.

use serde::Serialize;

use crate::model::FaceId;
use crate::options::AnalysisOptions;

/// What a solve did, measured.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SolveReport {
    /// The options the run used (after defaults were filled in).
    pub options: AnalysisOptions,
    /// Cell size the grid used, in model units.
    pub cell_size: Option<f64>,
    /// Face ids each named reference resolved to, in setup order.
    pub resolved_faces: Vec<ResolvedFace>,
    /// Active (non-empty) cells.
    pub cells_active: usize,
    /// Active cells cut by the boundary.
    pub cells_cut: usize,
    /// Small cells merged into an aggregate.
    pub cells_aggregated: usize,
    /// Degrees of freedom after aggregation.
    pub dof: usize,
    /// Multigrid levels.
    pub levels: usize,
    /// Dofs on each multigrid level, finest first.
    pub level_dofs: Vec<usize>,
    /// CG iterations taken.
    pub cg_iterations: usize,
    /// Relative residual after each CG iteration.
    pub residual_history: Vec<f64>,
    /// Rasterised against exact volume, per solid.
    pub volume_parity: Vec<VolumeParity>,
    /// Applied force against reaction, when a structural solve ran.
    pub force_balance: Option<ForceBalance>,
    /// Cut cells where EAS was skipped (volume fraction below η_EAS).
    pub eas_disabled_cells: usize,
    /// Cut cells whose exact moments were inconsistent and were taken
    /// from the sub-cell lattice instead (lattice-accurate only).
    pub moment_fallback_cells: usize,
    /// Cells split into fragments (parts of a solid closer than a cell).
    pub split_cells: usize,
    /// Split cells whose fragments keep lattice moments (no wall section
    /// could be split by fragment), so their moments are lattice-accurate
    /// only.
    pub split_cells_lattice_moments: usize,
    /// Nodes where stress recovery fell back from SPR.
    pub spr_fallbacks: usize,
    /// Per solid (structural): the rigid-mode matrix's eigenvalues over
    /// the diagonal stiffness scale, ascending: how firmly the supports
    /// hold each rigid motion.
    pub rigid_mode_ratios: Vec<[f64; 6]>,
    /// Passes of the radiation fixed point (0 without radiating faces).
    pub radiation_passes: usize,
    /// Lanczos steps of a modal solve, over every run (0 otherwise).
    pub lanczos_steps: usize,
    /// Lanczos runs of a modal solve: the first and its deflated
    /// restarts (0 otherwise).
    pub lanczos_runs: usize,
    /// Volume of the homogenised infill core, model units cubed (0
    /// without infill); the rest of the part is solid walls.
    pub infill_core_volume: f64,
    /// Relative density of the infill core, when the solve had one.
    pub infill_density: Option<f64>,
    /// Warnings, in the order they were raised.
    pub warnings: Vec<Warning>,
}

impl SolveReport {
    /// A singular system as the solver declines it: `UnderResolved` with
    /// the thinnest flagged wall when aggregation found a solid without
    /// a large cell (its Nitsche terms are then not stabilised; measured
    /// on cycloidal_joint, PH-053), else `SingularSystem`.
    pub fn singular(&self, detail: String) -> crate::error::AnalysisError {
        let thinnest = self
            .warnings
            .iter()
            .filter_map(|w| match w {
                Warning::UnderResolved { wall, cell, .. } => Some((*wall, *cell)),
                _ => None,
            })
            .min_by(|a, b| a.0.total_cmp(&b.0));
        match thinnest {
            Some((min_wall, cell)) => crate::error::AnalysisError::UnderResolved { min_wall, cell },
            None => crate::error::AnalysisError::SingularSystem { detail },
        }
    }

    /// An empty report for a run with `options`.
    pub fn new(options: AnalysisOptions) -> Self {
        SolveReport {
            options,
            cell_size: None,
            resolved_faces: Vec::new(),
            cells_active: 0,
            cells_cut: 0,
            cells_aggregated: 0,
            dof: 0,
            levels: 0,
            level_dofs: Vec::new(),
            cg_iterations: 0,
            residual_history: Vec::new(),
            volume_parity: Vec::new(),
            force_balance: None,
            eas_disabled_cells: 0,
            moment_fallback_cells: 0,
            split_cells: 0,
            split_cells_lattice_moments: 0,
            radiation_passes: 0,
            lanczos_steps: 0,
            lanczos_runs: 0,
            infill_core_volume: 0.0,
            infill_density: None,
            spr_fallbacks: 0,
            rigid_mode_ratios: Vec::new(),
            warnings: Vec::new(),
        }
    }
}

/// The faces one named reference resolved to.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct ResolvedFace {
    /// The name in the setup.
    pub name: String,
    /// Matched faces, in `FaceId` order.
    pub faces: Vec<FaceId>,
}

/// Rasterised volume of one solid, and the source kernel's volume (gate
/// X5 holds the rasterised volume to the tessellation's).
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct VolumeParity {
    /// Solid index.
    pub solid: u32,
    /// Sum of the cell volume moments.
    pub rasterised: f64,
    /// The source kernel's volume, when it gave one
    /// (`SolidMesh::exact_volume`: exact for planar solids, a fine
    /// tessellation otherwise).
    pub exact: Option<f64>,
}

/// Applied load against support reaction (gate X6).
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct ForceBalance {
    /// Sum of applied forces, newtons.
    pub applied: [f64; 3],
    /// Sum of support reactions, newtons.
    pub reaction: [f64; 3],
}

/// Something the user should know about a result that is not a decline.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub enum Warning {
    /// A face is referenced by raw id, which is fragile across edits.
    FaceIdReference {
        /// The name in the setup.
        name: String,
    },
    /// A solid has no cell with volume fraction at least γ: walls are
    /// thinner than about one cell, and γ was lowered for it.
    UnderResolved {
        /// Solid index.
        solid: u32,
        /// Wall thickness estimate 2V/A, model units.
        wall: f64,
        /// Cell size, model units.
        cell: f64,
    },
    /// Some curved faces' tessellation does not match the face: its
    /// area exceeds the kernel's upper bound by more than
    /// `AREA_MISMATCH` (PK-787: a sphere corner tessellated as most of
    /// the sphere), or it folds against the surface. Results near those
    /// faces are unreliable, and they take no part in wall readings.
    TessellationMismatch {
        /// Solid index.
        solid: u32,
        /// Mismatched faces, ascending.
        faces: Vec<u32>,
        /// Their tessellated area, model units squared.
        tessellated_area: f64,
        /// Their area folded against the surface, model units squared.
        folded_area: f64,
    },
    /// A solid's tessellated volume differs from the source kernel's
    /// volume by more than `tessellation::VOLUME_DEFICIT`: the solve
    /// sees the tessellation, so mass and stiffness carry the
    /// difference (thin curved features at the display tessellation).
    TessellationVolume {
        /// Solid index.
        solid: u32,
        /// The tessellation's enclosed volume, model units cubed.
        tessellated: f64,
        /// The kernel's volume (`SolidMesh::exact_volume`), model units
        /// cubed.
        kernel: f64,
    },
    /// Small cells sit more than `THIN_AGGREGATE_DISTANCE` cells from
    /// their root: a thin feature.
    ThinAggregates {
        /// Solid index.
        solid: u32,
        /// Largest BFS distance from a root.
        max_distance: usize,
    },
    /// The hottest point is within 10 °C of the heat deflection
    /// temperature: the part may soften, and linear-elastic results are
    /// invalid above about HDT (notes 13).
    ServiceTemperature {
        /// Highest temperature, °C.
        max_c: f64,
        /// The material's HDT, °C.
        hdt_c: f64,
    },
    /// A solid's supported area spans fewer than `SUPPORT_SPAN_CELLS`
    /// cells along some axis: it holds those rotations only through a
    /// patch the grid barely resolves, so the answer depends on the grid
    /// (measured, PH-033). Support a larger area or refine.
    PointLikeSupport {
        /// Solid index.
        solid: u32,
        /// Support span along its weakest axis, model units.
        span: f64,
        /// Cell size, model units.
        cell: f64,
        /// The rotations (about the support centroid) it barely holds.
        rotations: Vec<crate::structural::rigid::RigidMotion>,
    },
    /// Parts of one solid closer than the cell size share cells and are
    /// treated as joined (a slit narrower than a cell).
    FeaturesJoined {
        /// Solid index.
        solid: u32,
        /// Number of cells holding disconnected fragments.
        cells: usize,
    },
    /// A generative design's part is thinner than the filter's
    /// smallest member (2r) over much of its volume: no member there
    /// can be thinner than the part, so the design stays grey or rests
    /// on cut cells (measured, PH-125: a 1 mm plate under r = 0.75 mm
    /// ended half grey at aligned grid offsets).
    ThinDesignSpace {
        /// The share of the part's volume thinner than 2r.
        fraction: f64,
        /// Its median thickness, mm.
        thickness: f64,
        /// The filter's smallest member 2r, mm.
        member: f64,
    },
    /// A generative design's filter radius spans fewer than
    /// `FILTER_CELLS` cells: the design is not converged in the grid,
    /// finer cells give a different design (measured, PH-124: a crisp
    /// cantilever moved 9–11% from r/h = 1.5 and 2 to twice that, 2.2%
    /// from r/h = 3 to 4.5).
    FilterUnderResolved {
        /// Filter radius, mm.
        radius: f64,
        /// Cell size, mm.
        cell: f64,
    },
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::AnalysisError;

    #[test]
    fn a_singular_system_after_under_resolved_walls_declines_as_under_resolved() {
        let mut r = SolveReport::new(AnalysisOptions::default());
        // No under-resolved solid: the singular system is reported as is.
        assert!(matches!(
            r.singular("pivot".into()),
            AnalysisError::SingularSystem { .. }
        ));
        for (solid, wall) in [(0, 0.9), (3, 0.4), (5, 1.2)] {
            r.warnings.push(Warning::UnderResolved {
                solid,
                wall,
                cell: 14.0,
            });
        }
        assert_eq!(
            r.singular("pivot".into()),
            AnalysisError::UnderResolved {
                min_wall: 0.4,
                cell: 14.0
            }
        );
    }
}
