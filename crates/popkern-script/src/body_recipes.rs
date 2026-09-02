//! The full fillet/chamfer/shell recipes of spec §5 as scalar-binding
//! arguments: every skin spells a variable-radius law, a two-distance
//! chamfer, or a shell join style through these, and the recorder
//! binds the scalars to feature-scoped parameters exactly like the
//! constant-radius ops do.

use popkern_core::naming::EntityName;
use popkern_history::feature::{Feature, JoinStyleKind, ShellDirectionKind, TrimPolicyKind};
use popkern_history::feature_blend::{
    ChamferKindDesc, FilletRadiusLaw, FilletRecipeDesc, RadiusStation,
};
use popkern_history::unit::Unit;

use crate::body::Body;
use crate::error::{ScriptError, SelectionKind};
use crate::model::Model;
use crate::scalar::{Scalar, bind_scalar};
use crate::selection::{EdgeSelector, FaceSelector, exactly_one, resolve_edges, resolve_faces};

/// How a fillet's radius varies along each edge (spec §5 `law`).
#[derive(Debug, Clone, PartialEq)]
pub enum RadiusLaw {
    /// The feature's `radius` everywhere.
    Constant,
    /// Linear from the feature's `radius` at the edge start to `end`.
    Linear {
        /// Radius at the end vertex.
        end: Scalar,
    },
    /// Radius given at normalised stations `t` in `[0, 1]`.
    ByStations(Vec<(f64, Scalar)>),
}

/// Chamfer geometry beyond the equal-distance bevel (spec §5 `kind`).
#[derive(Debug, Clone, PartialEq)]
pub enum ChamferSpec {
    /// Equal-distance 45° bevel.
    Equal,
    /// The feature's `distance` on `face`, `distance_2` on the other support.
    TwoDistance {
        /// The reference face (must resolve to exactly one face).
        face: FaceSelector,
        /// Distance on the second support.
        distance_2: Scalar,
    },
    /// The feature's `distance` on `face`, at `angle` degrees from it.
    DistanceAngle {
        /// The reference face (must resolve to exactly one face).
        face: FaceSelector,
        /// Angle from the reference face, in degrees.
        angle: Scalar,
    },
}

/// Shell options beyond thickness (spec §5 `shell`): direction
/// defaults to `inward`, join to `arc`, trim to `trim`.
#[derive(Debug, Clone, Copy)]
pub struct ShellOptions {
    /// Offset direction.
    pub direction: ShellDirectionKind,
    /// Corner/edge join style.
    pub join: JoinStyleKind,
    /// Self-intersection policy.
    pub trim: TrimPolicyKind,
}

impl Default for ShellOptions {
    fn default() -> Self {
        ShellOptions {
            direction: ShellDirectionKind::Inward,
            join: JoinStyleKind::default(),
            trim: TrimPolicyKind::default(),
        }
    }
}

impl Model {
    /// Resolve a selection that must name exactly one edge of the
    /// active body (single-entity arguments, spec §6.1).
    pub fn resolve_edge(&mut self, selection: &EdgeSelector) -> Result<EntityName, ScriptError> {
        let body = self.active()?;
        let brep = self.body_brep(body.id)?;
        let labels = |l: &str| self.entity_by_label(l);
        exactly_one(
            resolve_edges(&brep, selection, &labels)?,
            SelectionKind::Edges,
        )
    }

    /// Resolve a selection that must name exactly one face of the
    /// active body.
    pub fn resolve_face(&mut self, selection: &FaceSelector) -> Result<EntityName, ScriptError> {
        let body = self.active()?;
        let brep = self.body_brep(body.id)?;
        let labels = |l: &str| self.entity_by_label(l);
        exactly_one(
            resolve_faces(&brep, selection, &labels)?,
            SelectionKind::Faces,
        )
    }

    fn bind(
        &mut self,
        role: &str,
        s: &Scalar,
        unit: Unit,
    ) -> Result<popkern_history::parameter::ParameterId, ScriptError> {
        bind_scalar(&mut self.tree, &mut self.pending_counter, role, s, unit)
    }

    /// Bind a radius law's scalars into the recipe.
    pub(crate) fn bind_radius_law(
        &mut self,
        law: &RadiusLaw,
        recipe: &mut FilletRecipeDesc,
    ) -> Result<(), ScriptError> {
        let unit = self.unit;
        recipe.radius_law = match law {
            RadiusLaw::Constant => FilletRadiusLaw::Constant,
            RadiusLaw::Linear { end } => FilletRadiusLaw::Linear {
                end: self.bind("radius_end", end, unit)?,
            },
            RadiusLaw::ByStations(stations) => {
                let mut bound = Vec::with_capacity(stations.len());
                for (i, (t, radius)) in stations.iter().enumerate() {
                    if !(0.0..=1.0).contains(t) {
                        return Err(ScriptError::Argument {
                            argument: "law".into(),
                            reason: format!("station {i} parameter {t} is outside [0, 1]"),
                        });
                    }
                    bound.push(RadiusStation {
                        t: *t,
                        radius: self.bind(&format!("radius_{i}"), radius, unit)?,
                    });
                }
                FilletRadiusLaw::ByStations(bound)
            }
        };
        Ok(())
    }
}

impl Body {
    /// Fillet with a full recipe and a radius law.
    pub fn fillet_law(
        self,
        model: &mut Model,
        edges: EdgeSelector,
        radius: impl Into<Scalar>,
        law: RadiusLaw,
        mut recipe: FilletRecipeDesc,
        label: &str,
    ) -> Result<Body, ScriptError> {
        model.bind_radius_law(&law, &mut recipe)?;
        self.fillet_with(model, edges, radius, recipe, label)
    }

    /// Chamfer with an explicit kind and label.
    pub fn chamfer_with(
        self,
        model: &mut Model,
        edges: EdgeSelector,
        distance: impl Into<Scalar>,
        kind: ChamferSpec,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let unit = model.unit;
        let kind = match kind {
            ChamferSpec::Equal => ChamferKindDesc::Equal,
            ChamferSpec::TwoDistance { face, distance_2 } => ChamferKindDesc::TwoDistance {
                face: model.resolve_face(&face)?,
                distance_2: model.bind("distance_2", &distance_2, unit)?,
            },
            ChamferSpec::DistanceAngle { face, angle } => ChamferKindDesc::DistanceAngle {
                face: model.resolve_face(&face)?,
                angle: model.bind("angle", &angle, Unit::Degrees)?,
            },
        };
        self.chamfer_kind(model, edges, distance, kind, label)
    }

    /// Shell with the full option set and a label.
    pub fn shell_with(
        self,
        model: &mut Model,
        faces_to_remove: FaceSelector,
        thickness: impl Into<Scalar>,
        options: ShellOptions,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let feature = |thickness, faces| Feature::Shell {
            faces_to_remove: faces,
            thickness,
            direction: options.direction,
            join_style: options.join,
            trim_policy: options.trim,
        };
        self.shell_feature(model, faces_to_remove, thickness.into(), feature, label)
    }
}
