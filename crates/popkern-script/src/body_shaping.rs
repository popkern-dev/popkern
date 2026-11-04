//! The offset family and draft (spec §5 preview): `thicken`,
//! `offset_solid` and `draft` recorded as `Feature::Thicken`,
//! `Feature::OffsetSolid` and `Feature::Draft`, their scalars bound to
//! feature-scoped parameters like every other operation's.

use popkern_history::feature::{Feature, JoinStyleKind, TrimPolicyKind};
use popkern_history::feature_draft::NeutralPlaneDesc;
use popkern_history::unit::Unit;

use crate::body::Body;
use crate::error::ScriptError;
use crate::model::Model;
use crate::scalar::Scalar;
use crate::selection::FaceSelector;

/// How an offset closes its gaps and treats self-intersecting regions
/// (the options `thicken` and `offset_solid` share with `shell`).
#[derive(Debug, Clone, Copy, Default)]
pub struct OffsetOptions {
    /// Corner/edge join style.
    pub join: JoinStyleKind,
    /// Self-intersection policy.
    pub trim: TrimPolicyKind,
}

/// A draft's neutral plane: the plane each drafted face keeps its
/// intersection with. Its normal is the pull direction.
#[derive(Debug, Clone, PartialEq)]
pub enum NeutralSpec {
    /// A planar face of the body, by query (its outward normal pulls).
    Face(FaceSelector),
    /// A literal plane: a point on it (in the model's unit) and its
    /// normal.
    Plane {
        /// A point on the plane.
        origin: [f64; 3],
        /// The plane normal (the pull direction).
        normal: [f64; 3],
    },
}

impl Body {
    /// Offset every face of the solid by `distance`: positive grows it,
    /// negative shrinks it.
    pub fn offset_solid_with(
        self,
        model: &mut Model,
        distance: impl Into<Scalar>,
        options: OffsetOptions,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let unit = model.unit;
        let distance = Self::bind(model, "distance", &distance.into(), unit)?;
        let feature = Feature::OffsetSolid {
            distance,
            join_style: options.join,
            trim_policy: options.trim,
        };
        self.record(model, feature, label, Vec::new())
    }

    /// Thicken the body to a wall `thickness` thick, half to each side:
    /// a sheet becomes a closed solid, a closed body a hollow one.
    pub fn thicken_with(
        self,
        model: &mut Model,
        thickness: impl Into<Scalar>,
        options: OffsetOptions,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let unit = model.unit;
        let thickness = Self::bind(model, "thickness", &thickness.into(), unit)?;
        let feature = Feature::Thicken {
            thickness,
            join_style: options.join,
            trim_policy: options.trim,
        };
        self.record(model, feature, label, Vec::new())
    }

    /// Draft `faces` by `angle` (degrees) about `neutral`.
    pub fn draft_with(
        self,
        model: &mut Model,
        faces: FaceSelector,
        angle: impl Into<Scalar>,
        neutral: NeutralSpec,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let snapshot = self.faces(model, &faces)?;
        let mut records = vec![Self::selection_record("faces", &snapshot)];
        let neutral_plane = match neutral {
            NeutralSpec::Face(selector) => NeutralPlaneDesc::Face(self.one_face_recorded(
                model,
                "neutral",
                &selector,
                &mut records,
            )?),
            // The origin is a position in the model's unit; the tree
            // stores metres, as every bound length is.
            NeutralSpec::Plane { origin, normal } => NeutralPlaneDesc::Literal {
                origin: origin.map(|c| model.unit.to_si(c)),
                normal,
            },
        };
        let angle = Self::bind(model, "angle", &angle.into(), Unit::Degrees)?;
        let feature = Feature::Draft {
            faces: snapshot.names,
            angle,
            neutral_plane,
        };
        self.record(model, feature, label, records)
    }
}
