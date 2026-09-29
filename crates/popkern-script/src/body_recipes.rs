//! The full fillet/chamfer/shell recipes of spec §5 as scalar-binding
//! arguments: every skin spells a variable-radius law, a two-distance
//! chamfer, or a shell join style through these, and the recorder
//! binds the scalars to feature-scoped parameters exactly like the
//! constant-radius ops do.

use popkern_core::naming::EntityName;
use popkern_history::feature::{
    AxisDesc, BoolOpKind, Feature, JoinStyleKind, RevolveEnd as TreeRevolveEnd, ShellDirectionKind,
    TrimPolicyKind,
};
use popkern_history::feature_blend::{
    ChamferKindDesc, FilletRadiusLaw, FilletRecipeDesc, RadiusStation,
};
use popkern_history::feature_transform::RigidMotionDesc;
use popkern_history::unit::Unit;

use crate::body::{Body, RevolveEnd};
use crate::error::{ScriptError, SelectionKind};
use crate::model::Model;
use crate::scalar::{Scalar, bind_scalar};
use crate::selection::{
    EdgeSelector, EntityRef, FaceSelector, exactly_one, resolve_edges, resolve_faces,
};

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

/// A revolve/rotation axis (spec §1.4): a literal, or an edge of the
/// active body selected by a query (recorded as the `axis` selection).
#[derive(Debug, Clone, PartialEq)]
pub enum AxisSpec {
    /// `axis_at(origin, direction)`.
    Literal {
        /// A point on the axis.
        origin: [f64; 3],
        /// The axis direction.
        direction: [f64; 3],
    },
    /// `axis(edge_selection)` — must resolve to exactly one edge.
    Edge(EdgeSelector),
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
        let labels = |r: EntityRef| self.resolve_entity_ref(r);
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
        let labels = |r: EntityRef| self.resolve_entity_ref(r);
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

    /// Fillet with several radius groups in one feature (PK-602): each
    /// `(selector, radius)` group's edges take its radius, so the edges
    /// meeting at a corner can carry different radii and the corner is
    /// built as one vertex blend. Groups must not share an edge; a
    /// group that selects nothing is dropped when its selector allows
    /// empty. `recipe.radius_law` must be constant (laws are per fillet).
    pub fn fillet_each(
        self,
        model: &mut Model,
        groups: Vec<(EdgeSelector, Scalar)>,
        mut recipe: FilletRecipeDesc,
        label: &str,
    ) -> Result<Body, ScriptError> {
        if recipe.radius_law != FilletRadiusLaw::Constant {
            return Err(ScriptError::Argument {
                argument: "law".into(),
                reason: "radius laws are per fillet; fillet_each groups take constant radii".into(),
            });
        }
        check_corner_cap(&recipe)?;
        if groups.is_empty() {
            return Err(ScriptError::Argument {
                argument: "groups".into(),
                reason: "fillet_each needs at least one [edge selector, radius] group".into(),
            });
        }
        // Resolve every group before recording anything.
        let mut resolved = Vec::with_capacity(groups.len());
        for (i, (selector, radius)) in groups.into_iter().enumerate() {
            let snapshot = self.edges(model, &selector)?;
            resolved.push((i, snapshot, radius));
        }
        // Every overlapping pair, with its shared-edge count.
        let mut overlaps: Vec<(usize, usize, usize)> = Vec::new();
        for (x, (i, a, _)) in resolved.iter().enumerate() {
            for (j, b, _) in resolved.iter().skip(x + 1) {
                let shared = a.names.iter().filter(|n| b.names.contains(n)).count();
                if shared > 0 {
                    overlaps.push((*i, *j, shared));
                }
            }
        }
        if !overlaps.is_empty() {
            let pairs: Vec<String> = overlaps
                .iter()
                .map(|(i, j, n)| {
                    format!(
                        "groups {i} and {j} ({n} edge{})",
                        if *n == 1 { "" } else { "s" }
                    )
                })
                .collect();
            return Err(ScriptError::Argument {
                argument: "groups".into(),
                reason: format!(
                    "fillet_each: an edge is selected by more than one group — {} — each edge takes one radius",
                    pairs.join(", ")
                ),
            });
        }
        resolved.retain(|(_, snapshot, _)| !snapshot.names.is_empty());
        if resolved.is_empty() {
            return Ok(self);
        }
        let unit = model.unit;
        let mut records = Vec::with_capacity(resolved.len());
        let mut edges: Vec<EntityName> = Vec::new();
        let mut base = None;
        for (k, (_, snapshot, radius)) in resolved.iter().enumerate() {
            let (role, edge_role) = if k == 0 {
                ("radius".to_string(), "edges".to_string())
            } else {
                (format!("radius_{k}"), format!("edges_{k}"))
            };
            let pid = model.bind(&role, radius, unit)?;
            records.push(Self::selection_record(&edge_role, snapshot));
            edges.extend(snapshot.names.iter().cloned());
            if k == 0 {
                base = Some(pid);
            } else {
                recipe
                    .edge_radii
                    .push(popkern_history::feature_blend::EdgeRadiusGroup {
                        edges: snapshot.names.clone(),
                        radius: pid,
                    });
            }
        }
        let radius = base.ok_or_else(|| ScriptError::Argument {
            argument: "groups".into(),
            reason: "no group selected an edge".into(),
        })?;
        let feature = Feature::Fillet {
            edges,
            radius,
            recipe,
        };
        self.record(model, feature, label, records)
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
        let mut records = Vec::new();
        let kind = match kind {
            ChamferSpec::Equal => ChamferKindDesc::Equal,
            ChamferSpec::TwoDistance { face, distance_2 } => ChamferKindDesc::TwoDistance {
                face: self.one_face_recorded(model, "face", &face, &mut records)?,
                distance_2: model.bind("distance_2", &distance_2, unit)?,
            },
            ChamferSpec::DistanceAngle { face, angle } => ChamferKindDesc::DistanceAngle {
                face: self.one_face_recorded(model, "face", &face, &mut records)?,
                angle: model.bind("angle", &angle, Unit::Degrees)?,
            },
        };
        self.chamfer_kind(model, edges, distance, kind, label, records)
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

    pub(crate) fn axis(
        self,
        model: &mut Model,
        axis: AxisSpec,
        records: &mut Vec<popkern_history::SelectionRecord>,
    ) -> Result<AxisDesc, ScriptError> {
        Ok(match axis {
            AxisSpec::Literal { origin, direction } => AxisDesc::Literal { origin, direction },
            AxisSpec::Edge(selection) => {
                AxisDesc::Edge(self.one_edge_recorded(model, "axis", &selection, records)?)
            }
        })
    }

    /// Revolve the body's last sketch about a queried or literal axis,
    /// with a label.
    pub fn revolve_with(
        self,
        model: &mut Model,
        axis: AxisSpec,
        end: RevolveEnd,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let face = self.last_sketch_face(model)?;
        let mut records = Vec::new();
        let axis = self.axis(model, axis, &mut records)?;
        let end_condition = match end {
            RevolveEnd::Blind(a) => TreeRevolveEnd::Blind {
                angle: Body::bind(model, "angle", &a, Unit::Degrees)?,
            },
            RevolveEnd::FullCircle => TreeRevolveEnd::FullCircle,
            RevolveEnd::MidPlane(a) => TreeRevolveEnd::MidPlane {
                half_angle: Body::bind(model, "half_angle", &a, Unit::Degrees)?,
            },
        };
        let feature = Feature::Revolve {
            face,
            axis,
            end_condition,
            mode: popkern_history::feature::OperationMode::NewBody,
        };
        self.record(model, feature, label, records)
    }

    /// Boolean with an explicit label.
    pub fn boolean_with(
        self,
        model: &mut Model,
        operation: BoolOpKind,
        tool: Body,
        label: &str,
    ) -> Result<Body, ScriptError> {
        self.boolean_labeled(model, operation, tool, label)
    }

    /// Rigid translation with a label.
    pub fn translate_with(
        self,
        model: &mut Model,
        dx: impl Into<Scalar>,
        dy: impl Into<Scalar>,
        dz: impl Into<Scalar>,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let unit = model.unit;
        let motion = RigidMotionDesc::Translate {
            dx: Body::bind(model, "dx", &dx.into(), unit)?,
            dy: Body::bind(model, "dy", &dy.into(), unit)?,
            dz: Body::bind(model, "dz", &dz.into(), unit)?,
        };
        self.record(model, Feature::Transform { motion }, label, Vec::new())
    }

    /// Rigid rotation about a queried or literal axis, with a label.
    pub fn rotate_with(
        self,
        model: &mut Model,
        axis: AxisSpec,
        angle: impl Into<Scalar>,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let mut records = Vec::new();
        let axis = self.axis(model, axis, &mut records)?;
        let motion = RigidMotionDesc::Rotate {
            axis,
            angle: Body::bind(model, "angle", &angle.into(), Unit::Degrees)?,
        };
        self.record(model, Feature::Transform { motion }, label, records)
    }
}

/// Refuse a corner-cap preference that cannot apply (PK-603) before
/// anything is recorded.
pub(crate) fn check_corner_cap(recipe: &FilletRecipeDesc) -> Result<(), ScriptError> {
    match recipe.corner_cap_conflict() {
        Some(reason) => Err(ScriptError::Argument {
            argument: "corner_cap".into(),
            reason: reason.into(),
        }),
        None => Ok(()),
    }
}
