//! A generative design's settings.

use serde::Serialize;

use crate::error::AnalysisError;

use super::print::BuildDirection;
use super::projection::ProjectionKind;

/// A generative design's settings.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct GenerativeOptions {
    /// Material to keep, as a fraction of the part's volume, in (0, 1].
    pub volume_fraction: f64,
    /// Density filter radius r, mm: the smallest member is about 2r.
    pub filter_radius: f64,
    /// SIMP exponent p (top88: 3).
    pub penalty: f64,
    /// Void stiffness ε relative to the material (top88: 1e-9).
    pub void_stiffness: f64,
    /// Largest change of a variable per update (top88: 0.2).
    pub move_limit: f64,
    /// Stop when no physical density changes by more than this (top88:
    /// 0.01, on the variables).
    pub tolerance: f64,
    /// Most updates.
    pub max_iterations: usize,
    /// Material within this distance of a kept face is held solid, mm
    /// (default: the filter radius).
    pub keep_thickness: f64,
    /// The Heaviside projection's final sharpness β (0, the default:
    /// no projection). β starts at 1 and doubles each time the design
    /// settles or after `projection_every` updates, until it reaches
    /// this. Experimental, not offered by the setup: under the
    /// optimality criteria update it stalled on 2 of 5 grid offsets of a
    /// thin cantilever (PH-118, `projection`); a robust update (MMA) is
    /// PH-122.
    pub projection_beta: f64,
    /// Updates at most between doublings of β.
    pub projection_every: usize,
    /// The projection's form.
    pub projection_kind: ProjectionKind,
    /// The update that moves the design.
    pub method: UpdateMethod,
    /// The build direction, when the design must print without
    /// supports (the print filter); `None`: no overhang control.
    pub build_direction: Option<BuildDirection>,
    /// The robust formulation's threshold η in (0, ½) (Wang, Lazarov &
    /// Sigmund 2011 §5), with the threshold projection and MMA: the
    /// eroded design (threshold 1 − η) is solved and its compliance, the
    /// largest of the eroded, intermediate and dilated designs',
    /// minimised; the intermediate design (the projection's own
    /// threshold) is the result and meets the volume. `None`: off.
    /// The setup offers it as `robust` (η = 0.3): on a crisp cantilever
    /// it doubled the thinnest members (0.5 → 1.0 mm) for 16–23% of the
    /// stiffness, and moved 7.8% under refinement from 3 to 4.5 cells
    /// per filter radius where the plain projection moved 2.2% (PH-124). (Wang bound
    /// the dilated design's volume, rescaled every 20 updates; with this
    /// MMA's exact volume the rescale jolted the design and its eroded
    /// design came apart at β = 32.)
    pub robust: Option<f64>,
    /// CG's relative residual while the design moves (the analysis's
    /// `cg_rel_tol` when tighter than it is used); the final design is
    /// solved again at the analysis tolerance.
    pub update_tolerance: f64,
}

/// How each update moves the design.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum UpdateMethod {
    /// The optimality criteria update (top88 eqn 3).
    OptimalityCriteria,
    /// The method of moving asymptotes (Svanberg 1987).
    Mma,
}

impl GenerativeOptions {
    /// top88's settings for a volume fraction and filter radius.
    pub fn new(volume_fraction: f64, filter_radius: f64) -> GenerativeOptions {
        GenerativeOptions {
            volume_fraction,
            filter_radius,
            penalty: 3.0,
            void_stiffness: 1e-9,
            move_limit: 0.2,
            tolerance: 0.01,
            max_iterations: 200,
            keep_thickness: filter_radius,
            projection_beta: 0.0,
            projection_every: 50,
            projection_kind: ProjectionKind::Guest,
            method: UpdateMethod::OptimalityCriteria,
            build_direction: None,
            robust: None,
            update_tolerance: 0.0,
        }
    }

    pub(super) fn validate(&self) -> Result<(), AnalysisError> {
        let bad = |what| Err(AnalysisError::InvalidGenerative { what });
        let finite = |x: f64| x.is_finite();
        if !(finite(self.volume_fraction)
            && self.volume_fraction > 0.0
            && self.volume_fraction <= 1.0)
        {
            return bad("volume_fraction must be in (0, 1]");
        }
        if !(finite(self.filter_radius) && self.filter_radius > 0.0) {
            return bad("filter_radius must be positive");
        }
        if !(finite(self.penalty) && self.penalty >= 1.0) {
            return bad("penalty must be at least 1");
        }
        if !(finite(self.void_stiffness) && self.void_stiffness > 0.0 && self.void_stiffness < 1.0)
        {
            return bad("void_stiffness must be in (0, 1)");
        }
        if !(finite(self.move_limit) && self.move_limit > 0.0 && self.move_limit <= 1.0) {
            return bad("move_limit must be in (0, 1]");
        }
        if !(finite(self.keep_thickness) && self.keep_thickness >= 0.0) {
            return bad("keep_thickness must be at least 0");
        }
        if !(finite(self.projection_beta) && self.projection_beta >= 0.0) {
            return bad("projection_beta must be at least 0");
        }
        if self.projection_every == 0 {
            return bad("projection_every must be at least 1");
        }
        if !(finite(self.tolerance) && self.tolerance > 0.0) {
            return bad("tolerance must be positive");
        }
        if !(finite(self.update_tolerance) && (0.0..1.0).contains(&self.update_tolerance)) {
            return bad("update_tolerance must be in [0, 1)");
        }
        if let Some(eta) = self.robust {
            if !(finite(eta) && eta > 0.0 && eta < 0.5) {
                return bad("the robust threshold must be in (0, 1/2)");
            }
            let threshold = matches!(self.projection_kind, ProjectionKind::Threshold(_));
            if !(threshold && self.projection_beta > 0.0 && self.method == UpdateMethod::Mma) {
                return bad("the robust formulation needs the threshold projection and MMA");
            }
        }
        Ok(())
    }
}
