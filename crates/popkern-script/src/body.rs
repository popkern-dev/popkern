//! Feature operations on a body (spec §5). Each op resolves its
//! selections against the live body, records one `Feature`, replays,
//! and returns the served body.

use popkern_core::naming::{EntityName, EntityRole, RoleCategory};
use popkern_history::feature::{
    AxisDesc, BoolOpKind, ExtrudeEnd, Feature, JoinStyleKind, OperationMode,
    RevolveEnd as TreeRevolveEnd, ShellDirectionKind, TrimPolicyKind, Vec3Desc,
};
use popkern_history::feature_blend::{ChamferKindDesc, ContinuityKind, FilletRecipeDesc};
use popkern_history::feature_transform::RigidMotionDesc;
use popkern_history::node::BodyId;
use popkern_history::unit::Unit;

use crate::error::{ScriptError, SelectionKind};
use crate::model::{FeatureHandle, Model};
use crate::record::record_feature;
use crate::scalar::{Scalar, bind_scalar};
use crate::selection::{EdgeSelector, FaceSelector, exactly_one, resolve_edges, resolve_faces};

/// A body of the document. Cheap to copy; all state lives in the model.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Body {
    /// The tree's body id.
    pub id: BodyId,
}

/// An extrude end condition (spec §5).
#[derive(Debug, Clone, PartialEq)]
pub enum End {
    /// A fixed distance.
    Blind(Scalar),
    /// Through everything.
    ThroughAll,
    /// Up to a face.
    UpToSurface(FaceSelector),
    /// Offset from a face.
    OffsetFromSurface(FaceSelector, Scalar),
    /// Symmetric about the sketch plane.
    MidPlane(Scalar),
}

/// A revolve end condition.
#[derive(Debug, Clone, PartialEq)]
pub enum RevolveEnd {
    /// A fixed angle (degrees at the surface).
    Blind(Scalar),
    /// Full 360°.
    FullCircle,
    /// Symmetric half angle.
    MidPlane(Scalar),
}

impl Body {
    pub(crate) fn new(id: BodyId) -> Body {
        Body { id }
    }

    /// The face of this body's most recent sketch (the default
    /// profile face for extrude/revolve).
    pub fn last_sketch_face(&self, model: &Model) -> Result<EntityName, ScriptError> {
        model
            .tree
            .timeline()
            .iter()
            .rev()
            .filter_map(|&id| model.tree.feature(id))
            .find(|n| n.body == self.id && matches!(n.feature, Feature::Sketch { .. }))
            .map(|n| EntityName {
                operation: n.operation_id,
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: 0,
                },
                ancestors: Vec::new(),
            })
            .ok_or_else(|| ScriptError::Argument {
                argument: "face".into(),
                reason: "body has no sketch to extrude".into(),
            })
    }

    /// The most recent feature of this body.
    pub fn last(&self, model: &Model) -> Option<FeatureHandle> {
        model
            .tree
            .timeline()
            .iter()
            .rev()
            .find(|&&id| model.tree.feature(id).map(|n| n.body) == Some(self.id))
            .and_then(|&id| model.handle(id).ok())
    }

    fn record(self, model: &mut Model, feature: Feature, label: &str) -> Result<Body, ScriptError> {
        let id = record_feature(
            &mut model.tree,
            model.tolerance,
            feature,
            label,
            Some(self.id),
        )?;
        let body = model.tree.feature(id).map(|n| n.body).unwrap_or(self.id);
        model.active = Some(body);
        Ok(Body::new(body))
    }

    fn brep(self, model: &mut Model) -> Result<popkern_core::topology::BRep, ScriptError> {
        model.body_brep(self.id)
    }

    fn edges(
        self,
        model: &mut Model,
        selector: &EdgeSelector,
    ) -> Result<crate::selection::SelectionSnapshot, ScriptError> {
        let brep = self.brep(model)?;
        let labels = |l: &str| model.entity_by_label(l);
        resolve_edges(&brep, selector, &labels)
    }

    fn faces(
        self,
        model: &mut Model,
        selector: &FaceSelector,
    ) -> Result<crate::selection::SelectionSnapshot, ScriptError> {
        let brep = self.brep(model)?;
        let labels = |l: &str| model.entity_by_label(l);
        resolve_faces(&brep, selector, &labels)
    }

    fn one_face(
        self,
        model: &mut Model,
        selector: &FaceSelector,
    ) -> Result<EntityName, ScriptError> {
        exactly_one(self.faces(model, selector)?, SelectionKind::Faces)
    }

    fn bind(
        model: &mut Model,
        role: &str,
        s: &Scalar,
        unit: Unit,
    ) -> Result<popkern_history::parameter::ParameterId, ScriptError> {
        bind_scalar(&mut model.tree, &mut model.pending_counter, role, s, unit)
    }

    fn extrude_end(self, model: &mut Model, end: End) -> Result<ExtrudeEnd, ScriptError> {
        let unit = model.unit;
        Ok(match end {
            End::Blind(d) => ExtrudeEnd::Blind {
                distance: Self::bind(model, "distance", &d, unit)?,
            },
            End::ThroughAll => ExtrudeEnd::ThroughAll,
            End::UpToSurface(sel) => ExtrudeEnd::UpToSurface {
                face: self.one_face(model, &sel)?,
            },
            End::OffsetFromSurface(sel, d) => ExtrudeEnd::OffsetFromSurface {
                face: self.one_face(model, &sel)?,
                offset: Self::bind(model, "offset", &d, unit)?,
            },
            End::MidPlane(d) => ExtrudeEnd::MidPlane {
                half_distance: Self::bind(model, "half_distance", &d, unit)?,
            },
        })
    }

    /// Extrude the body's last sketch.
    pub fn extrude(self, model: &mut Model, end: End) -> Result<Body, ScriptError> {
        self.extrude_with(model, end, OperationMode::NewBody, "")
    }

    /// Extrude with an explicit mode and label.
    pub fn extrude_with(
        self,
        model: &mut Model,
        end: End,
        mode: OperationMode,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let face = self.last_sketch_face(model)?;
        let end_condition = self.extrude_end(model, end)?;
        let feature = Feature::Extrude {
            face,
            direction: None,
            end_condition,
            mode,
            draft_angle: None,
            second_direction: None,
        };
        self.record(model, feature, label)
    }

    /// Revolve the body's last sketch about `axis`.
    pub fn revolve(
        self,
        model: &mut Model,
        axis: AxisDesc,
        end: RevolveEnd,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let face = self.last_sketch_face(model)?;
        let end_condition = match end {
            RevolveEnd::Blind(a) => TreeRevolveEnd::Blind {
                angle: Self::bind(model, "angle", &a, Unit::Degrees)?,
            },
            RevolveEnd::FullCircle => TreeRevolveEnd::FullCircle,
            RevolveEnd::MidPlane(a) => TreeRevolveEnd::MidPlane {
                half_angle: Self::bind(model, "half_angle", &a, Unit::Degrees)?,
            },
        };
        let feature = Feature::Revolve {
            face,
            axis,
            end_condition,
            mode: OperationMode::NewBody,
        };
        self.record(model, feature, label)
    }

    /// Constant-radius fillet with the default recipe (G1, rolling ball).
    pub fn fillet(
        self,
        model: &mut Model,
        edges: EdgeSelector,
        radius: impl Into<Scalar>,
    ) -> Result<Body, ScriptError> {
        self.fillet_with(model, edges, radius, FilletRecipeDesc::default(), "")
    }

    /// Fillet with a full recipe and label.
    pub fn fillet_with(
        self,
        model: &mut Model,
        edges: EdgeSelector,
        radius: impl Into<Scalar>,
        recipe: FilletRecipeDesc,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let snapshot = self.edges(model, &edges)?;
        if snapshot.names.is_empty() {
            return Ok(self);
        }
        let unit = model.unit;
        let radius = Self::bind(model, "radius", &radius.into(), unit)?;
        let feature = Feature::Fillet {
            edges: snapshot.names,
            radius,
            recipe,
        };
        self.record(model, feature, label)
    }

    /// Fillet with a given continuity.
    pub fn fillet_continuity(
        self,
        model: &mut Model,
        edges: EdgeSelector,
        radius: impl Into<Scalar>,
        continuity: ContinuityKind,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let recipe = FilletRecipeDesc {
            continuity,
            ..Default::default()
        };
        self.fillet_with(model, edges, radius, recipe, label)
    }

    /// Equal-distance chamfer.
    pub fn chamfer(
        self,
        model: &mut Model,
        edges: EdgeSelector,
        distance: impl Into<Scalar>,
    ) -> Result<Body, ScriptError> {
        let snapshot = self.edges(model, &edges)?;
        if snapshot.names.is_empty() {
            return Ok(self);
        }
        let unit = model.unit;
        let distance = Self::bind(model, "distance", &distance.into(), unit)?;
        let feature = Feature::Chamfer {
            edges: snapshot.names,
            distance,
            kind: ChamferKindDesc::Equal,
        };
        self.record(model, feature, "")
    }

    /// Shell: remove `faces` and hollow to `thickness`.
    pub fn shell(
        self,
        model: &mut Model,
        faces_to_remove: FaceSelector,
        thickness: impl Into<Scalar>,
        direction: ShellDirectionKind,
    ) -> Result<Body, ScriptError> {
        let faces = self.faces(model, &faces_to_remove)?.names;
        let unit = model.unit;
        let thickness = Self::bind(model, "thickness", &thickness.into(), unit)?;
        let feature = Feature::Shell {
            faces_to_remove: faces,
            thickness,
            direction,
            join_style: JoinStyleKind::default(),
            trim_policy: TrimPolicyKind::default(),
        };
        self.record(model, feature, "")
    }

    fn boolean(self, model: &mut Model, op: BoolOpKind, tool: Body) -> Result<Body, ScriptError> {
        if tool.id == self.id {
            return Err(ScriptError::Argument {
                argument: "tool".into(),
                reason: "a body cannot be its own tool".into(),
            });
        }
        let feature = Feature::Boolean {
            operation: op,
            tool_body: tool.id,
        };
        self.record(model, feature, "")
    }

    /// Boolean union with `tool` (consumes the tool body).
    pub fn union(self, model: &mut Model, tool: Body) -> Result<Body, ScriptError> {
        self.boolean(model, BoolOpKind::Union, tool)
    }

    /// Boolean subtraction of `tool`.
    pub fn subtract(self, model: &mut Model, tool: Body) -> Result<Body, ScriptError> {
        self.boolean(model, BoolOpKind::Subtraction, tool)
    }

    /// Boolean intersection with `tool`.
    pub fn intersect(self, model: &mut Model, tool: Body) -> Result<Body, ScriptError> {
        self.boolean(model, BoolOpKind::Intersection, tool)
    }

    /// Rigid translation.
    pub fn translate(
        self,
        model: &mut Model,
        dx: impl Into<Scalar>,
        dy: impl Into<Scalar>,
        dz: impl Into<Scalar>,
    ) -> Result<Body, ScriptError> {
        let unit = model.unit;
        let motion = RigidMotionDesc::Translate {
            dx: Self::bind(model, "dx", &dx.into(), unit)?,
            dy: Self::bind(model, "dy", &dy.into(), unit)?,
            dz: Self::bind(model, "dz", &dz.into(), unit)?,
        };
        self.record(model, Feature::Transform { motion }, "")
    }

    /// Rigid rotation about `axis` by `angle` degrees.
    pub fn rotate(
        self,
        model: &mut Model,
        axis: AxisDesc,
        angle: impl Into<Scalar>,
    ) -> Result<Body, ScriptError> {
        let motion = RigidMotionDesc::Rotate {
            axis,
            angle: Self::bind(model, "angle", &angle.into(), Unit::Degrees)?,
        };
        self.record(model, Feature::Transform { motion }, "")
    }

    /// The literal direction helper for `Vec3Desc`.
    pub fn dir(v: [f64; 3]) -> Vec3Desc {
        Vec3Desc::Literal(v)
    }
}
