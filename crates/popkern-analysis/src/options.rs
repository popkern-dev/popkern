//! Discretisation and solver options.
//!
//! Every tolerance and budget the analysis uses is a field here, is
//! passed explicitly, and is echoed in the `SolveReport`. There is no
//! global state.

use serde::Serialize;

/// Discretisation and solver options.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct AnalysisOptions {
    /// Grid cell size in model units; derived from the bounding box and
    /// `max_active_cells` when `None`.
    pub target_cell_size: Option<f64>,
    /// Budget of active (non-empty) cells. Default 30,000 for native
    /// runs (up to about 300,000 fit in memory); the browser uses
    /// [`WASM_CELL_CAP`].
    pub max_active_cells: usize,
    /// Sub-cell lattice depth for the fallback classifier: the lattice
    /// is 2^depth per axis. Default 2 (a 4³ lattice).
    pub subcell_depth: u8,
    /// AgFEM threshold γ: a cell is large when its volume fraction is
    /// at least γ. Default 0.25 (implementation map §1.3; tuned by E2).
    pub aggregation_threshold: f64,
    /// Support penalty factor p in β = p·E/h. Default 10³
    /// (implementation map §1.4).
    pub support_penalty_factor: f64,
    /// Nitsche factor γ_D for weakly imposed fixed values (consistent,
    /// unlike a penalty): β = γ_D·k/h for conduction and
    /// β = γ_D·(λ + 2μ)/h for elasticity. Default 10, as in the FEMPAR
    /// AgFEM runs (notes 07).
    pub nitsche_factor: f64,
    /// Relative residual at which CG stops. Default 1e-8 (the gate
    /// tolerance; display runs may use 1e-6).
    pub cg_rel_tol: f64,
    /// CG iteration budget. Default 200 (gates expect fewer than 60).
    pub cg_max_iters: usize,
    /// Shift of the grid origin against the geometry, mm (each ≥ 0).
    /// Default zero. Gates use it to test robustness to where the grid
    /// falls; translating the part instead moves the grid with it.
    pub grid_shift: [f64; 3],
    /// Multigrid coarsening stops once a level has at most this many
    /// dofs; the coarsest level is factored densely. Default 500
    /// (measured, PH-032: on a 30k-cell cube 500 and 1500 tie at 5.8 s
    /// and 4096 takes 8.2 s; on a slender cut beam 150 needs 21 PCG
    /// iterations against 14 at 500).
    pub coarsest_dofs: usize,
}

impl Default for AnalysisOptions {
    fn default() -> Self {
        AnalysisOptions {
            target_cell_size: None,
            max_active_cells: 30_000,
            subcell_depth: 2,
            aggregation_threshold: 0.25,
            support_penalty_factor: 1.0e3,
            nitsche_factor: 10.0,
            cg_rel_tol: 1.0e-8,
            cg_max_iters: 200,
            coarsest_dofs: 500,
            grid_shift: [0.0; 3],
        }
    }
}

/// The cell budget for solves in the browser (gate E3, PH-030, PH-048).
/// Measured 2026-10-11 in headless Chrome on block_two_bores, one solve
/// per page, fastest of repeated runs on a shared machine, after
/// PH-048's work (shared full-cell operator, pattern assembly, solids'
/// matrices stacked without a sort, face forces only where supported):
/// budget 10k → 8,400 cells, 1.95 s, 181 MB; 15k → 12,528 cells, 2.9 s,
/// 192 MB; 20k → 17,739 cells, 3.4 s, 362 MB; 30k → 26,722 cells,
/// 4.7 s, 385 MB (natively 3.6 s).
pub const WASM_CELL_CAP: usize = 30_000;

/// Largest sub-cell depth accepted (a 64³ lattice per cell).
pub const MAX_SUBCELL_DEPTH: u8 = 6;

/// An option outside its valid range.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[error("{field} = {value} is out of range ({range})")]
pub struct OptionError {
    /// The option's field name.
    pub field: &'static str,
    /// The value given, as text.
    pub value: String,
    /// The valid range.
    pub range: &'static str,
}

fn out_of_range(field: &'static str, value: impl ToString, range: &'static str) -> OptionError {
    OptionError {
        field,
        value: value.to_string(),
        range,
    }
}

impl AnalysisOptions {
    /// Checks every option against its valid range.
    pub fn validate(&self) -> Result<(), OptionError> {
        if let Some(h) = self.target_cell_size
            && !(h.is_finite() && h > 0.0)
        {
            return Err(out_of_range("target_cell_size", h, "> 0"));
        }
        if self.max_active_cells == 0 {
            return Err(out_of_range("max_active_cells", 0, ">= 1"));
        }
        if self.subcell_depth > MAX_SUBCELL_DEPTH {
            return Err(out_of_range("subcell_depth", self.subcell_depth, "0..=6"));
        }
        let g = self.aggregation_threshold;
        if !(g.is_finite() && (0.0..1.0).contains(&g)) {
            return Err(out_of_range("aggregation_threshold", g, "[0, 1)"));
        }
        let p = self.support_penalty_factor;
        if !(p.is_finite() && p > 0.0) {
            return Err(out_of_range("support_penalty_factor", p, "> 0"));
        }
        let g = self.nitsche_factor;
        if !(g.is_finite() && g > 0.0) {
            return Err(out_of_range("nitsche_factor", g, "> 0"));
        }
        let t = self.cg_rel_tol;
        if !(t.is_finite() && t > 0.0 && t < 1.0) {
            return Err(out_of_range("cg_rel_tol", t, "(0, 1)"));
        }
        if self.cg_max_iters == 0 {
            return Err(out_of_range("cg_max_iters", 0, ">= 1"));
        }
        if let Some(&s) = self
            .grid_shift
            .iter()
            .find(|s| !(s.is_finite() && **s >= 0.0))
        {
            return Err(out_of_range("grid_shift", s, ">= 0"));
        }
        if self.coarsest_dofs == 0 {
            return Err(out_of_range("coarsest_dofs", 0, ">= 1"));
        }
        Ok(())
    }
}
