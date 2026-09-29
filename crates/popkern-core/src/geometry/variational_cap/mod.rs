//! Variational G2 cap over an n-sided opening (PK-574).
//!
//! Given the n sides of an opening, each a curve lying on a known
//! neighbour surface, this module fits n tensor-product B-spline
//! sectors that meet each other G2 exactly and meet the neighbours
//! G1/G2 to a stated tolerance.
//!
//! # Construction
//!
//! Each side is split at its parameter midpoint `m_s`. Sector `s` is the
//! quad `[o, m_s, c_s, m_{s+1}]`: `S_s(0,0)` is the shared centre `o`,
//! `S_s(u,0)` the spoke to `m_s`, `S_s(0,v)` the spoke to `m_{s+1}`,
//! `S_s(1,v)` the second half of side `s` and `S_s(u,1)` the first half
//! of side `s+1`, reversed. All sectors share one clamped knot vector
//! `K` in both directions, so a spoke is one control row seen from two
//! sectors.
//!
//! - **Spokes (exact).** Loop & Schaefer's correspondence-map
//!   conditions (CGF 27(5), 2008, eqs. (8)–(10)) with the
//!   characteristic map φ_n make adjacent sectors G2 along the whole
//!   spoke, including at the n-valent centre. They are linear and
//!   data-independent, imposed per knot span by collocation (exact for
//!   polynomial identities) and eliminated through a memoised null
//!   space.
//! - **Rim G0 (exact to the approximant).** Each side half becomes a
//!   B-spline row in `K` that interpolates the corners and the side's
//!   2-jet at `m_s`, least-squares elsewhere ([`rim`]).
//! - **Rim G1/G2 (fitted).** Normal and transverse-curvature agreement
//!   with the oriented neighbour, Gauss–Newton over frozen neighbour
//!   jets, weighted against a metric thin-plate energy ([`energy`]).
//! - **Refinement.** Knots are inserted into `K` where the sampled rim
//!   residual exceeds its target, within a budget; the best level is
//!   returned, flagged when it misses the targets.
//!
//! The module knows nothing about fillets or topology: turning sectors
//! into faces, and certifying the result, is the caller's job (the
//! setback rung in popkern-blend).

mod chart;
mod energy;
mod energy_mcv;
mod error;
mod fit;
mod gauss_newton;
mod jets;
mod level;
pub(crate) mod linalg;
mod measure;
mod nullspace;
mod rim;
mod rim_rows;
mod solve;
mod space;
mod spoke;

#[cfg(test)]
mod tests;

pub use error::VariationalCapError;
pub use fit::fit_variational_cap;

use super::Point3;
use super::nurbs_curve::NurbsCurve;
use super::nurbs_surface::NurbsSurface;
use super::surface::SurfaceGeometry;

/// One side of the opening, in loop order: counter-clockwise about the
/// cap's outward normal. `curve` runs from this side's start corner to
/// its end corner; consecutive sides share corners.
#[derive(Debug, Clone)]
pub struct VariationalCapSide {
    /// The side's boundary curve (start corner → end corner).
    pub curve: NurbsCurve,
    /// The surface the cap must meet G1/G2 along this side.
    pub neighbour: SurfaceGeometry,
    /// Whether the neighbour face's outward normal is the reverse of
    /// the surface's natural normal `S_u × S_v`. The G1/G2 terms
    /// compare against the oriented normal and second fundamental form,
    /// so a sign slip cannot read as a curvature match.
    pub neighbour_reversed: bool,
}

/// Residual targets the fit refines towards.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VariationalCapTargets {
    /// Positional gap between the rim and the neighbour (model units).
    pub g0: f64,
    /// Sine of the normal angle across the rim.
    pub g1: f64,
    /// Normal-curvature jump across the rim (1/length); `None` fits G1
    /// only.
    pub g2: Option<f64>,
}

impl VariationalCapTargets {
    /// Targets for a G1 fit.
    pub fn g1(g0: f64, g1: f64) -> Self {
        VariationalCapTargets { g0, g1, g2: None }
    }

    /// Targets for a G2 fit.
    pub fn g2(g0: f64, g1: f64, g2: f64) -> Self {
        VariationalCapTargets {
            g0,
            g1,
            g2: Some(g2),
        }
    }
}

/// Which metric the thin-plate energy uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CapEnergy {
    /// The characteristic map's metric throughout: scale-free, the same
    /// shape for every sector, the regular thin plate at n = 4.
    CharacteristicThinPlate,
    /// `∫ tr(G⁻¹BG⁻¹B)`, B the second derivatives' normal part against a
    /// frozen reference: fourth order, the energy for a G1-only fit.
    Bending,
    /// `∫ |∇_g H|²`, H the mean curvature with the reference normal and
    /// metric frozen (Schneider–Kobbelt's ΔH = 0 fairness): zero on
    /// every constant-mean-curvature surface, so a sphere is an exact
    /// minimiser. Sixth order: it needs the G2 rim rows, so a G1-only fit
    /// (targets without `g2`) uses [`CapEnergy::Bending`] instead.
    MeanCurvatureVariation,
}

/// Construction budget and energy choice.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VariationalCapOptions {
    /// Polynomial degree of every sector in both directions.
    pub degree: usize,
    /// Knot-insertion rounds after the first solve.
    pub max_refinement_levels: usize,
    /// Largest solve (free coefficients over all three coordinates).
    pub max_unknowns: usize,
    /// Gauss–Newton iterations per refinement level.
    pub max_gauss_newton: usize,
    /// Energy metric.
    pub energy: CapEnergy,
    /// Weight of the thin-plate term relative to the rim terms, in
    /// normalised units.
    pub fairness_weight: f64,
    /// Weight of the membrane term relative to the thin-plate term.
    pub membrane_weight: f64,
    /// Weight of the parametric thin plate that regularises the
    /// parameterisation under the bending energy (which sees no
    /// tangential motion).
    pub parameter_regularisation: f64,
    /// Outer passes re-freezing the energy's reference at the current
    /// surface.
    pub reference_passes: usize,
    /// Weight of the tangential anchor to the level's start surface,
    /// relative to the (trace-balanced) regulariser: it holds the
    /// parameterisation without resisting the shape.
    pub anchor_weight: f64,
    /// Weight of the transversality prior on the rim cross derivatives
    /// (relative residual units; the G1 rows outweigh it by ~1e6).
    pub transversal_weight: f64,
}

impl Default for VariationalCapOptions {
    fn default() -> Self {
        VariationalCapOptions {
            degree: 7,
            max_refinement_levels: 4,
            max_unknowns: 2400,
            max_gauss_newton: 16,
            energy: CapEnergy::MeanCurvatureVariation,
            fairness_weight: 1.0e-2,
            membrane_weight: 1.0e-2,
            parameter_regularisation: 1.0e-3,
            reference_passes: 1,
            anchor_weight: 1.0,
            transversal_weight: 1.0e-1,
        }
    }
}

/// Which half of a side a residual belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RimHalf {
    /// From the side's start corner to its midpoint.
    Leading,
    /// From the side's midpoint to its end corner.
    Trailing,
}

/// The worst sampled rim residual on one side half.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RimResidual {
    /// Side index.
    pub side: usize,
    /// Which half.
    pub half: RimHalf,
    /// Side-curve parameter of the worst sample.
    pub parameter: f64,
    /// Worst positional gap to the neighbour on this half.
    pub g0: f64,
    /// Worst normal-angle sine on this half.
    pub g1: f64,
    /// Worst normal-curvature jump (spectral norm of the second
    /// fundamental form difference) on this half.
    pub g2: f64,
}

/// Construction evidence for a fitted cap.
#[derive(Debug, Clone, PartialEq)]
pub struct VariationalCapReport {
    /// The shared knot vector `K` (both directions, every sector).
    pub knots: Vec<f64>,
    /// Knot-insertion rounds used.
    pub refinement_levels: usize,
    /// Free coefficients in the final solve (all coordinates).
    pub unknowns: usize,
    /// Gauss–Newton iterations over all levels.
    pub gauss_newton_iterations: usize,
    /// Worst sampled residual per side half (2n entries).
    pub rim: Vec<RimResidual>,
    /// Worst sampled (g0, g1, g2) across the spokes.
    pub worst_spoke: (f64, f64, f64),
    /// Smallest `|S_u × S_v|` over a dense grid, relative to the mean.
    pub min_regularity: f64,
    /// Final thin-plate energy (normalised units).
    pub energy: f64,
    /// Whether every sampled rim residual met its target. A fit that
    /// ran out of refinement budget returns its best level with this
    /// `false`; the caller decides what the achieved residuals are worth
    /// (a G2 attempt may still certify G1).
    pub targets_met: bool,
}

impl VariationalCapReport {
    /// The worst rim residual against `targets`, by the largest ratio
    /// of residual to target over g0, g1 and (when fitted) g2.
    pub fn worst_rim(&self, targets: &VariationalCapTargets) -> Option<RimResidual> {
        self.rim
            .iter()
            .copied()
            .max_by(|a, b| excess(a, targets).total_cmp(&excess(b, targets)))
    }

    /// Whether every rim half meets `targets`.
    pub fn meets(&self, targets: &VariationalCapTargets) -> bool {
        self.rim.iter().all(|r| excess(r, targets) <= 1.0)
    }
}

/// Largest residual-to-target ratio of one sample.
pub(crate) fn excess(r: &RimResidual, targets: &VariationalCapTargets) -> f64 {
    let mut e = (r.g0 / targets.g0).max(r.g1 / targets.g1);
    if let Some(g2) = targets.g2 {
        e = e.max(r.g2 / g2);
    }
    e
}

/// A fitted cap: one sector surface per side.
#[derive(Debug, Clone)]
pub struct VariationalCap {
    /// Sector `s` (see the module docs for its corners and edges).
    pub sectors: Vec<NurbsSurface>,
    /// The shared centre `o`.
    pub centre: Point3,
    /// `m_s`: the midpoint of side `s`, where spoke `s` meets the rim.
    pub spoke_points: Vec<Point3>,
    /// Construction evidence.
    pub report: VariationalCapReport,
}
