//! Hinge recorder verbs (PK-693): `hinge` and `living_hinge` on
//! [`Model`].
//!
//! A hinge records one feature per component (the tree builds one body
//! per feature) sharing one parameter set, owned by the first component
//! recorded. Dimensions left unset follow the length by expression
//! (leaf width L/4, thickness L/30, pin L/15, …), so a length edit
//! rescales the hinge; a dimension given becomes its own slider.

use popkern_history::expression::{BinOp, Expression};
use popkern_history::feature::Feature;
use popkern_history::feature_hinge::{
    HingeComponentDesc, HingeDesc, HingeHolesDesc, HingeKindDesc, LivingHingeDesc,
    LivingHingeStyleDesc,
};
use popkern_history::parameter::ParameterId;
use popkern_history::unit::Unit;

use crate::body::Body;
use crate::body_fastener::ThreadSpecArg;
use crate::error::ScriptError;
use crate::model::Model;
use crate::record::record_feature;
use crate::scalar::Scalar;

/// Screw holes for [`HingeOptions`].
#[derive(Debug, Clone)]
pub struct HingeHolesArg {
    /// Holes per leaf.
    pub per_leaf: Scalar,
    /// The screw.
    pub thread: ThreadSpecArg,
    /// Sink a 90° countersink.
    pub countersink: bool,
}

/// Options for [`Model::hinge`]; `None` follows the length.
#[derive(Debug, Clone, Default)]
pub struct HingeOptions {
    /// Each leaf's width (default L/4).
    pub leaf_width: Option<Scalar>,
    /// Leaf thickness (default L/30).
    pub thickness: Option<Scalar>,
    /// Knuckle count (default 5).
    pub knuckles: Option<Scalar>,
    /// Pin diameter (default L/15).
    pub pin: Option<Scalar>,
    /// Knuckle outer diameter (default: a leaf rolled round the pin).
    pub knuckle_diameter: Option<Scalar>,
    /// Radial clearance (default 0 butt, L/200 print-in-place).
    pub clearance: Option<Scalar>,
    /// Axial knuckle gap (default L/400 butt, L/200 print-in-place).
    pub gap: Option<Scalar>,
    /// Open angle in degrees (default 180, flat).
    pub angle: Option<Scalar>,
    /// Screw holes.
    pub holes: Option<HingeHolesArg>,
    /// Record only this component (an assembly part holds one body).
    pub only: Option<HingeComponentDesc>,
}

/// The bodies [`Model::hinge`] recorded.
#[derive(Debug, Clone, Copy)]
pub struct HingeBodies {
    /// Leaf A (toward −X).
    pub leaf_a: Option<Body>,
    /// Leaf B (toward +X, at the open angle).
    pub leaf_b: Option<Body>,
    /// The separate pin (butt hinges).
    pub pin: Option<Body>,
}

/// Options for [`Model::living_hinge`]; `None` follows the span.
#[derive(Debug, Clone, Default)]
pub struct LivingHingeOptions {
    /// Panel width beyond the flex zone (default 2·span).
    pub panel_width: Option<Scalar>,
    /// Panel thickness (default span/5).
    pub thickness: Option<Scalar>,
    /// Straight: web thickness (default thickness/5).
    pub web: Option<Scalar>,
    /// Straight: root fillet radius (default half the room above the
    /// web, at most span/4).
    pub root_radius: Option<Scalar>,
    /// Slotted: slot pitch (default span/8).
    pub slot_pitch: Option<Scalar>,
    /// Slotted: slot width (default span/24).
    pub slot_width: Option<Scalar>,
    /// Slotted: bridge at each slot's closed end (default length/10).
    pub bridge: Option<Scalar>,
}

impl Model {
    /// `scalar` bound, or — when absent — a parameter driven by
    /// `of × factor` (so it follows `of`).
    /// `scalar` bound, or — absent — a parameter driven by `of × factor`.
    pub(crate) fn bind_or_follow(
        &mut self,
        role: &str,
        scalar: Option<&Scalar>,
        of: ParameterId,
        factor: f64,
        unit: Unit,
    ) -> Result<ParameterId, ScriptError> {
        match scalar {
            Some(s) => self.bind(role, s, unit),
            None => {
                let value = self
                    .tree
                    .parameters()
                    .get(of)
                    .map(|p| p.value)
                    .unwrap_or(0.0);
                let id = self.bind(role, &Scalar::Literal(unit.from_si(value * factor)), unit)?;
                self.tree.set_parameter_expression(
                    id,
                    Expression::BinOp(
                        Box::new(Expression::Param(of)),
                        BinOp::Mul,
                        Box::new(Expression::Literal(factor)),
                    ),
                )?;
                Ok(id)
            }
        }
    }

    fn bind_or(
        &mut self,
        role: &str,
        scalar: Option<&Scalar>,
        default: f64,
        unit: Unit,
    ) -> Result<ParameterId, ScriptError> {
        let s = scalar.cloned().unwrap_or(Scalar::Literal(default));
        self.bind(role, &s, unit)
    }

    /// Record a hinge (PK-693): leaf A, leaf B and — for a butt hinge —
    /// the pin, as bodies in the hinge frame (pin axis +Z from z = 0,
    /// leaf A toward −X, leaf B at the open angle). `options.only`
    /// records one component. Labels: `label` (or `hinge-N`) for the
    /// first component, then `<label>_b`, `<label>_pin`.
    ///
    /// # Errors
    ///
    /// [`ScriptError::Argument`] for a print-in-place hinge asked for its
    /// pin; replay failures (the hinge's typed declines) as for any
    /// feature.
    pub fn hinge(
        &mut self,
        kind: HingeKindDesc,
        length: impl Into<Scalar>,
        options: &HingeOptions,
        label: &str,
    ) -> Result<HingeBodies, ScriptError> {
        if kind == HingeKindDesc::PrintInPlace && options.only == Some(HingeComponentDesc::Pin) {
            return Err(ScriptError::Argument {
                argument: "only".into(),
                reason: "a print-in-place hinge has no separate pin (it is part of leaf A)".into(),
            });
        }
        let unit = self.unit;
        let o = options;
        let print = kind == HingeKindDesc::PrintInPlace;
        let length = self.bind("length", &length.into(), unit)?;
        let spec = HingeDesc {
            length,
            leaf_width: self.bind_or_follow(
                "leaf_width",
                o.leaf_width.as_ref(),
                length,
                0.25,
                unit,
            )?,
            thickness: self.bind_or_follow(
                "thickness",
                o.thickness.as_ref(),
                length,
                1.0 / 30.0,
                unit,
            )?,
            knuckles: self.bind_or("knuckles", o.knuckles.as_ref(), 5.0, Unit::Ratio)?,
            pin_diameter: self.bind_or_follow(
                "pin_diameter",
                o.pin.as_ref(),
                length,
                1.0 / 15.0,
                unit,
            )?,
            knuckle_diameter: o
                .knuckle_diameter
                .as_ref()
                .map(|s| self.bind("knuckle_diameter", s, unit))
                .transpose()?,
            clearance: if print {
                self.bind_or_follow("clearance", o.clearance.as_ref(), length, 1.0 / 200.0, unit)?
            } else {
                self.bind_or("clearance", o.clearance.as_ref(), 0.0, unit)?
            },
            gap: self.bind_or_follow(
                "gap",
                o.gap.as_ref(),
                length,
                if print { 1.0 / 200.0 } else { 1.0 / 400.0 },
                unit,
            )?,
            angle: self.bind_or("angle", o.angle.as_ref(), 180.0, Unit::Degrees)?,
            holes: match &o.holes {
                Some(h) => Some(HingeHolesDesc {
                    per_leaf: self.bind("holes", &h.per_leaf, Unit::Ratio)?,
                    thread: crate::body_fastener::bind_thread_desc(self, &h.thread)?,
                    countersink: h.countersink,
                }),
                None => None,
            },
        };
        let mut components = vec![HingeComponentDesc::LeafA, HingeComponentDesc::LeafB];
        if !print {
            components.push(HingeComponentDesc::Pin);
        }
        if let Some(only) = o.only {
            components.retain(|c| *c == only);
        }
        let mut bodies = HingeBodies {
            leaf_a: None,
            leaf_b: None,
            pin: None,
        };
        let mut owner_label: Option<String> = None;
        for component in components {
            let body = match &owner_label {
                // The first component recorded owns the parameters.
                None => {
                    let feature = Feature::Hinge {
                        kind,
                        component,
                        spec: spec.clone(),
                        owner: true,
                    };
                    let id = record_feature(
                        &mut self.tree,
                        self.tolerance,
                        feature,
                        label,
                        None,
                        Vec::new(),
                        &mut self.geometry_current,
                    )?;
                    self.note_source(id);
                    let node = self
                        .tree
                        .feature(id)
                        .ok_or(popkern_history::HistoryError::FeatureNotFound)?;
                    owner_label = Some(node.label.clone());
                    self.active = Some(node.body);
                    Body::new(node.body)
                }
                Some(owner) => {
                    let owner = owner.clone();
                    self.hinge_component(&owner, component, "")?
                }
            };
            match component {
                HingeComponentDesc::LeafA => bodies.leaf_a = Some(body),
                HingeComponentDesc::LeafB => bodies.leaf_b = Some(body),
                HingeComponentDesc::Pin => bodies.pin = Some(body),
            }
        }
        Ok(bodies)
    }

    /// Record a living hinge (PK-698): one body, two panels and a
    /// flexure across a flex zone `span` wide, `length` along the bend
    /// axis (+Z), mounting face y = 0.
    ///
    /// # Errors
    ///
    /// Replay failures (the hinge's typed declines) as for any feature.
    pub fn living_hinge(
        &mut self,
        style: LivingHingeStyleDesc,
        span: impl Into<Scalar>,
        length: impl Into<Scalar>,
        options: &LivingHingeOptions,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let unit = self.unit;
        let o = options;
        let length = self.bind("length", &length.into(), unit)?;
        let span = self.bind("span", &span.into(), unit)?;
        let panel_width =
            self.bind_or_follow("panel_width", o.panel_width.as_ref(), span, 2.0, unit)?;
        let thickness = self.bind_or_follow("thickness", o.thickness.as_ref(), span, 0.2, unit)?;
        let (mut web, mut root_radius, mut slot_pitch, mut slot_width, mut bridge) =
            (None, None, None, None, None);
        match style {
            LivingHingeStyleDesc::Straight => {
                let w = self.bind_or_follow("web", o.web.as_ref(), thickness, 0.2, unit)?;
                // Half the room above the default web: 0.4 × thickness.
                root_radius = Some(self.bind_or_follow(
                    "root_radius",
                    o.root_radius.as_ref(),
                    thickness,
                    0.4,
                    unit,
                )?);
                web = Some(w);
            }
            LivingHingeStyleDesc::Slotted => {
                slot_pitch = Some(self.bind_or_follow(
                    "slot_pitch",
                    o.slot_pitch.as_ref(),
                    span,
                    0.125,
                    unit,
                )?);
                slot_width = Some(self.bind_or_follow(
                    "slot_width",
                    o.slot_width.as_ref(),
                    span,
                    1.0 / 24.0,
                    unit,
                )?);
                bridge =
                    Some(self.bind_or_follow("bridge", o.bridge.as_ref(), length, 0.1, unit)?);
            }
        }
        let feature = Feature::LivingHinge {
            style,
            spec: LivingHingeDesc {
                length,
                span,
                panel_width,
                thickness,
                web,
                root_radius,
                slot_pitch,
                slot_width,
                bridge,
            },
        };
        let id = record_feature(
            &mut self.tree,
            self.tolerance,
            feature,
            label,
            None,
            Vec::new(),
            &mut self.geometry_current,
        )?;
        self.note_source(id);
        let body = self
            .tree
            .feature(id)
            .map(|n| n.body)
            .ok_or(popkern_history::HistoryError::FeatureNotFound)?;
        self.active = Some(body);
        Ok(Body::new(body))
    }
}

impl Model {
    /// Record another component of the hinge whose owning feature is
    /// labelled `owner` (it shares that hinge's parameters): what a
    /// hinge recorded component by component — and document-as-code —
    /// uses. Labels `<owner>_b` / `<owner>_pin` unless `label` is given.
    ///
    /// # Errors
    ///
    /// [`ScriptError::Argument`] when `owner` labels no hinge's owning
    /// feature, the component already exists, or a print-in-place hinge
    /// is asked for a pin.
    pub fn hinge_component(
        &mut self,
        owner: &str,
        component: HingeComponentDesc,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let bad = |reason: String| ScriptError::Argument {
            argument: "hinge".into(),
            reason,
        };
        let mut found = None;
        let mut taken = Vec::new();
        for &id in self.tree.timeline() {
            let Some(node) = self.tree.feature(id) else {
                continue;
            };
            if let Feature::Hinge {
                kind,
                component: c,
                spec,
                owner: is_owner,
            } = &node.feature
            {
                if *is_owner && node.label == owner {
                    found = Some((*kind, spec.clone()));
                }
                taken.push((spec.length, *c));
            }
        }
        let Some((kind, spec)) = found else {
            return Err(bad(format!(
                "no hinge labelled '{owner}' owns its parameters"
            )));
        };
        if kind == HingeKindDesc::PrintInPlace && component == HingeComponentDesc::Pin {
            return Err(bad("a print-in-place hinge has no separate pin".into()));
        }
        if taken.contains(&(spec.length, component)) {
            return Err(bad(format!(
                "hinge '{owner}' already has its {}",
                component.spelling()
            )));
        }
        let label = if label.is_empty() {
            match component {
                HingeComponentDesc::LeafB => format!("{owner}_b"),
                other => format!("{owner}_{}", other.spelling()),
            }
        } else {
            label.to_string()
        };
        let feature = Feature::Hinge {
            kind,
            component,
            spec,
            owner: false,
        };
        let id = record_feature(
            &mut self.tree,
            self.tolerance,
            feature,
            &label,
            None,
            Vec::new(),
            &mut self.geometry_current,
        )?;
        self.note_source(id);
        let body = self
            .tree
            .feature(id)
            .map(|n| n.body)
            .ok_or(popkern_history::HistoryError::FeatureNotFound)?;
        self.active = Some(body);
        Ok(Body::new(body))
    }
}
