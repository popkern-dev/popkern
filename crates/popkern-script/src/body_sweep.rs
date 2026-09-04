//! Sweep-family recorder verbs (spec §5): `sweep`, `pipe`, and
//! `helical_sweep` on [`Body`]; `loft` on [`Model`] (sections are
//! sketch bodies, so there is no receiving body).

use popkern_history::feature::{Feature, OperationMode};
use popkern_history::unit::Unit;

use crate::body::Body;
use crate::body_recipes::AxisSpec;
use crate::error::ScriptError;
use crate::model::Model;
use crate::scalar::Scalar;
use crate::selection::EdgeSelector;

/// Helix parameters for [`Body::helical_sweep`]. Radius and pitch are
/// lengths in the model unit; `start_angle` is in degrees; `turns` is
/// a dimensionless count (fractional allowed).
#[derive(Debug, Clone)]
pub struct HelixParams {
    /// Helix radius at the start.
    pub radius: Scalar,
    /// Radius at the end for a conical helix; `None` keeps `radius`.
    pub taper_to: Option<Scalar>,
    /// Axial advance per turn.
    pub pitch: Scalar,
    /// Number of turns.
    pub turns: Scalar,
    /// Phase: angle about the axis at which the helix starts, in
    /// degrees — what phases multiple helical bodies apart (propeller
    /// blades, multi-start threads).
    pub start_angle: Scalar,
    /// Winding sense. `false` = right-hand.
    pub left_handed: bool,
}

impl HelixParams {
    /// A right-hand cylindrical helix with zero phase.
    pub fn new(
        radius: impl Into<Scalar>,
        pitch: impl Into<Scalar>,
        turns: impl Into<Scalar>,
    ) -> HelixParams {
        HelixParams {
            radius: radius.into(),
            taper_to: None,
            pitch: pitch.into(),
            turns: turns.into(),
            start_angle: Scalar::from(0.0),
            left_handed: false,
        }
    }
}

/// Options for [`Model::loft`].
#[derive(Debug, Clone, Default)]
pub struct LoftOptions {
    /// Wrap the loft around so the last section connects to the first.
    pub closed: bool,
    /// New body (default), boss, or cut.
    pub mode: OperationMode,
}

impl Body {
    /// Sweep the body's last sketch along `path` — a single edge of
    /// `path_body` (the profile sketch body has no path edges of its
    /// own). The profile contributes its shape; placement comes from
    /// the path. The path body is referenced, not consumed.
    pub fn sweep(
        self,
        model: &mut Model,
        path_body: Body,
        path: EdgeSelector,
    ) -> Result<Body, ScriptError> {
        self.sweep_with(model, path_body, path, 0.0, OperationMode::NewBody, "")
    }

    /// Sweep with a total twist angle (degrees), mode, and label.
    pub fn sweep_with(
        self,
        model: &mut Model,
        path_body: Body,
        path: EdgeSelector,
        twist: impl Into<Scalar>,
        mode: OperationMode,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let face = self.last_sketch_face(model)?;
        let mut records = Vec::new();
        let path_edge = path_body.one_edge_recorded(model, "path_edge", &path, &mut records)?;
        let feature = Feature::Sweep {
            face,
            path_edge,
            twist: Body::bind(model, "twist", &twist.into(), Unit::Degrees)?,
            mode,
        };
        self.record(model, feature, label, records)
    }

    /// Create a pipe (tube) of `radius` along `path` (a single edge of
    /// this body).
    pub fn pipe(
        self,
        model: &mut Model,
        path: EdgeSelector,
        radius: impl Into<Scalar>,
    ) -> Result<Body, ScriptError> {
        self.pipe_with(model, path, radius, OperationMode::NewBody, "")
    }

    /// Pipe with an explicit mode and label.
    pub fn pipe_with(
        self,
        model: &mut Model,
        path: EdgeSelector,
        radius: impl Into<Scalar>,
        mode: OperationMode,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let unit = model.unit;
        let mut records = Vec::new();
        let path_edge = self.one_edge_recorded(model, "path_edge", &path, &mut records)?;
        let feature = Feature::Pipe {
            path_edge,
            radius: Body::bind(model, "radius", &radius.into(), unit)?,
            mode,
        };
        self.record(model, feature, label, records)
    }

    /// Sweep the body's last sketch along a helix about `axis`
    /// (springs, coils, auger flights, thread forms). The profile
    /// contributes its shape; placement comes from the helix.
    pub fn helical_sweep(
        self,
        model: &mut Model,
        axis: AxisSpec,
        helix: HelixParams,
    ) -> Result<Body, ScriptError> {
        self.helical_sweep_with(model, axis, helix, OperationMode::NewBody, "")
    }

    /// Helical sweep with an explicit mode and label.
    pub fn helical_sweep_with(
        self,
        model: &mut Model,
        axis: AxisSpec,
        helix: HelixParams,
        mode: OperationMode,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let unit = model.unit;
        let face = self.last_sketch_face(model)?;
        let mut records = Vec::new();
        let axis = self.axis(model, axis, &mut records)?;
        let end = helix.taper_to.unwrap_or_else(|| helix.radius.clone());
        let feature = Feature::HelicalSweep {
            face,
            axis,
            start_radius: Body::bind(model, "start_radius", &helix.radius, unit)?,
            end_radius: Body::bind(model, "end_radius", &end, unit)?,
            pitch: Body::bind(model, "pitch", &helix.pitch, unit)?,
            turns: Body::bind(model, "turns", &helix.turns, Unit::Ratio)?,
            start_angle: Body::bind(model, "start_angle", &helix.start_angle, Unit::Degrees)?,
            left_handed: helix.left_handed,
            mode,
        };
        self.record(model, feature, label, records)
    }
}

impl Model {
    /// Loft through `sections` in order — sketch bodies, each
    /// contributing its most recent sketch face. The loft records on
    /// the first section's body; the other section bodies are
    /// consumed building it.
    pub fn loft(&mut self, sections: &[Body]) -> Result<Body, ScriptError> {
        self.loft_with(sections, LoftOptions::default(), "")
    }

    /// Loft with options and a label.
    pub fn loft_with(
        &mut self,
        sections: &[Body],
        options: LoftOptions,
        label: &str,
    ) -> Result<Body, ScriptError> {
        if sections.len() < 2 {
            return Err(ScriptError::Argument {
                argument: "sections".into(),
                reason: format!("loft needs at least two sections (got {})", sections.len()),
            });
        }
        let names = sections
            .iter()
            .map(|b| b.last_sketch_face(self))
            .collect::<Result<Vec<_>, _>>()?;
        let feature = Feature::Loft {
            sections: names,
            closed: options.closed,
            mode: options.mode,
        };
        sections[0].record(self, feature, label, Vec::new())
    }
}
