//! Generative design on a body (PK-956,
//! docs/architecture/generative-features-plan.md): the body is the
//! design space; the recipe's supports, loads and kept faces are face
//! selections on it, recorded so the document can show them back as
//! code and an edit re-resolves them.

use popkern_history::feature::Feature;
use popkern_history::feature_generative::{
    BuildDirectionDesc, GenerativeDesc, GenerativeLoadDesc, GenerativeSupportDesc,
    GenerativeSupportKindDesc,
};
use popkern_history::unit::Unit;

use crate::body::Body;
use crate::error::ScriptError;
use crate::model::Model;
use crate::scalar::Scalar;
use crate::selection::FaceSelector;

/// A total force on the faces a selector picks.
#[derive(Debug, Clone)]
pub struct GenerativeLoad {
    /// The faces.
    pub faces: FaceSelector,
    /// Total force, N.
    pub force: Scalar,
    /// Direction it acts in (non-zero).
    pub direction: [f64; 3],
}

/// A generative design's recipe, as a script states it.
#[derive(Debug, Clone)]
pub struct GenerativeSpec {
    /// Material preset name (`"PLA (FDM)"`, ...).
    pub material: String,
    /// Supports: the faces and what they hold.
    pub supports: Vec<(FaceSelector, GenerativeSupportKindDesc)>,
    /// Loads.
    pub loads: Vec<GenerativeLoad>,
    /// Faces whose nearby material is kept.
    pub keep: Vec<FaceSelector>,
    /// Material kept, as a fraction of the design space.
    pub volume_fraction: Scalar,
    /// Density filter radius, in the model's unit (members are about
    /// twice it).
    pub filter: Scalar,
    /// How deep kept material reaches (default: the filter radius).
    pub keep_depth: Option<Scalar>,
    /// Analysis cell size (default: a third of the filter radius).
    pub cell: Option<Scalar>,
    /// The build direction the design prints along without supports.
    pub build_direction: Option<BuildDirectionDesc>,
    /// Most optimisation updates.
    pub max_iterations: u32,
    /// Thicken the thinnest members (the robust formulation).
    pub robust: bool,
    /// Most active analysis cells (default: the solver's).
    pub max_cells: Option<u32>,
}

/// Millimetres per stored length unit: a mm, cm or inch document stores
/// metres, a default-unit document millimetre magnitudes (PK-462, the
/// rule simulation and the STEP export follow).
fn mm_per_unit(unit: Unit) -> f64 {
    match unit {
        Unit::Millimeters | Unit::Centimeters | Unit::Inches => 1000.0,
        _ => 1.0,
    }
}

impl Body {
    /// Replace this body (the design space) with the stiffest design at
    /// `spec`'s volume fraction, with a label.
    ///
    /// # Errors
    ///
    /// A selection that picks nothing, a bad value, or the design's
    /// decline (typed: `GenerativeDeclined` with the solver's kind).
    pub fn generative(
        self,
        model: &mut Model,
        spec: &GenerativeSpec,
        label: &str,
    ) -> Result<Body, ScriptError> {
        self.require_alive(model)?;
        if spec.supports.is_empty() || spec.loads.is_empty() {
            return Err(ScriptError::Argument {
                argument: if spec.supports.is_empty() {
                    "supports".into()
                } else {
                    "loads".into()
                },
                reason: "a generative design needs at least one support and one load".into(),
            });
        }
        let unit = model.unit;
        let mut records = Vec::new();
        let mut supports = Vec::new();
        for (i, (selector, kind)) in spec.supports.iter().enumerate() {
            let snapshot = self.faces(model, selector)?;
            records.push(Self::selection_record(&format!("supports[{i}]"), &snapshot));
            supports.push(GenerativeSupportDesc {
                faces: snapshot.names,
                kind: *kind,
            });
        }
        let mut loads = Vec::new();
        for (i, load) in spec.loads.iter().enumerate() {
            let snapshot = self.faces(model, &load.faces)?;
            records.push(Self::selection_record(&format!("loads[{i}]"), &snapshot));
            loads.push(GenerativeLoadDesc {
                faces: snapshot.names,
                force: Self::bind(model, "force", &load.force, Unit::Ratio)?,
                direction: load.direction,
            });
        }
        let mut keep = Vec::new();
        for (i, selector) in spec.keep.iter().enumerate() {
            let snapshot = self.faces(model, selector)?;
            records.push(Self::selection_record(&format!("keep[{i}]"), &snapshot));
            keep.push(snapshot.names);
        }
        let desc = GenerativeDesc {
            material: spec.material.clone(),
            supports,
            loads,
            keep,
            volume_fraction: Self::bind(
                model,
                "volume_fraction",
                &spec.volume_fraction,
                Unit::Ratio,
            )?,
            filter: Self::bind(model, "filter", &spec.filter, unit)?,
            keep_depth: spec
                .keep_depth
                .as_ref()
                .map(|s| Self::bind(model, "keep_depth", s, unit))
                .transpose()?,
            cell: spec
                .cell
                .as_ref()
                .map(|s| Self::bind(model, "cell", s, unit))
                .transpose()?,
            build_direction: spec.build_direction,
            max_iterations: spec.max_iterations,
            robust: spec.robust,
            max_cells: spec.max_cells,
            mm_per_unit: mm_per_unit(unit),
        };
        self.record(
            model,
            Feature::Generative {
                spec: Box::new(desc),
            },
            label,
            records,
        )
    }
}
