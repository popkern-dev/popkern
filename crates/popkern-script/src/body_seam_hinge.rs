//! Seam hinge recorder verbs (PK-702): knuckles grown on two existing
//! bodies along the seam where they meet (a box and its lid), plus a
//! butt hinge's pin as a new body.
//!
//! The seam spans two bodies and a feature resolves selections against
//! its own body only, so the seam — a point, the seam direction, and the
//! outward and up directions — is resolved when the hinge is recorded
//! and stored. Up points from body A's bounds centre toward body B's;
//! outward is perpendicular to the seam and up, away from both.

use popkern_core::geometry::{Point3, Vec3};
use popkern_history::feature::Feature;
use popkern_history::feature_hinge::{HingeComponentDesc, HingeKindDesc, SeamDesc, SeamHingeDesc};
use popkern_history::unit::Unit;

use crate::body::Body;
use crate::body_hinge::HingeBodies;
use crate::error::ScriptError;
use crate::model::Model;
use crate::record::record_feature;
use crate::scalar::Scalar;

/// Options for [`Model::hinge_bodies`]; `None` follows the length.
#[derive(Debug, Clone, Default)]
pub struct SeamHingeOptions {
    /// Knuckle count (default 5).
    pub knuckles: Option<Scalar>,
    /// Pin diameter (default L/15).
    pub pin: Option<Scalar>,
    /// Knuckle wall round the bore (default L/30).
    pub wall: Option<Scalar>,
    /// Knuckle outer diameter (default pin + 2·clearance + 2·wall).
    pub knuckle_diameter: Option<Scalar>,
    /// Radial clearance (default L/200 print-in-place, L/400 butt — a
    /// seam hinge always needs room for the bodies to turn).
    pub clearance: Option<Scalar>,
    /// Axial gap between knuckles (default as `clearance`).
    pub gap: Option<Scalar>,
}

fn bounds_centre(model: &mut Model, body: Body) -> Result<Point3, ScriptError> {
    let brep = model.brep_of(body)?;
    let mut lo = [f64::MAX; 3];
    let mut hi = [f64::MIN; 3];
    for v in brep.vertices.values() {
        let p = [v.position.x, v.position.y, v.position.z];
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    if lo[0] > hi[0] {
        return Err(ScriptError::Argument {
            argument: "body".into(),
            reason: "the body has no geometry".into(),
        });
    }
    Ok(Point3::new(
        0.5 * (lo[0] + hi[0]),
        0.5 * (lo[1] + hi[1]),
        0.5 * (lo[2] + hi[2]),
    ))
}

impl Model {
    /// Record a hinge on bodies `a` and `b` along the seam line through
    /// `origin` in direction `direction`, `length` long from `origin`:
    /// knuckles grown on both bodies (A's below the seam, B's above) and
    /// — for a butt hinge — a pin body. Labels: `label` (or a minted one)
    /// on body A's feature, `<label>_b`, `<label>_pin`.
    ///
    /// # Errors
    ///
    /// [`ScriptError::Argument`] for a seam direction of zero length or
    /// bodies on the seam line (no up direction); the hinge's typed
    /// declines on replay.
    /// The end points of the straight edge `selector` picks on `body`
    /// (a seam given as an edge: the box's top back edge).
    ///
    /// # Errors
    ///
    /// The selection's errors (not exactly one edge); an edge that is
    /// not straight.
    pub fn seam_edge(
        &mut self,
        body: Body,
        selector: &crate::selection::EdgeSelector,
    ) -> Result<(Point3, Point3), ScriptError> {
        let mut records = Vec::new();
        let name = body.one_edge_recorded(self, "seam", selector, &mut records)?;
        let brep = self.brep_of(body)?;
        let edge = brep
            .edges
            .values()
            .find(|e| e.name.and_then(|n| brep.name_pool.get(n)) == Some(&name))
            .ok_or_else(|| ScriptError::Argument {
                argument: "seam".into(),
                reason: "the selected edge is not on the body".into(),
            })?;
        let straight = matches!(
            edge.curve.as_ref().map(|c| &c.geometry),
            Some(popkern_core::geometry::curve::CurveGeometry::Line { .. })
        );
        if !straight {
            return Err(ScriptError::Argument {
                argument: "seam".into(),
                reason: "a hinge seam is a straight edge".into(),
            });
        }
        let he = &brep.half_edges[edge.half_edge];
        let a = brep.vertices[he.origin].position;
        let b = brep.vertices[brep.half_edges[he.twin].origin].position;
        Ok((a, b))
    }

    /// The seam between bodies `a` and `b` along the line through
    /// `origin` in `direction`: up from A's bounds centre toward B's,
    /// outward away from both.
    ///
    /// # Errors
    ///
    /// [`ScriptError::Argument`] for a zero direction or bodies whose
    /// centres lie along the seam.
    pub fn seam_between(
        &mut self,
        a: Body,
        b: Body,
        origin: Point3,
        direction: Vec3,
    ) -> Result<SeamDesc, ScriptError> {
        let bad = |argument: &str, reason: &str| ScriptError::Argument {
            argument: argument.into(),
            reason: reason.into(),
        };
        let axis = direction
            .normalized()
            .ok_or_else(|| bad("axis", "the seam direction has zero length"))?;
        let (ca, cb) = (bounds_centre(self, a)?, bounds_centre(self, b)?);
        let across = cb - ca;
        let up = (across - axis * across.dot(&axis))
            .normalized()
            .ok_or_else(|| {
                bad(
                    "bodies",
                    "the two bodies' centres lie along the seam: no up side",
                )
            })?;
        let mut outward = axis.cross(&up);
        let mid = Point3::new(
            0.5 * (ca.x + cb.x),
            0.5 * (ca.y + cb.y),
            0.5 * (ca.z + cb.z),
        );
        if (origin - mid).dot(&outward) < 0.0 {
            outward = -outward;
        }
        Ok(SeamDesc {
            origin: [origin.x, origin.y, origin.z],
            axis: [axis.x, axis.y, axis.z],
            outward: [outward.x, outward.y, outward.z],
            up: [up.x, up.y, up.z],
        })
    }

    /// Record a seam hinge's first component — body A's knuckles on
    /// `body` — owning the hinge's parameters (bound from `length` and
    /// `options`). Its siblings follow with
    /// [`Model::seam_hinge_component`].
    ///
    /// # Errors
    ///
    /// Binding errors; the hinge's typed declines on replay.
    pub fn seam_hinge(
        &mut self,
        kind: HingeKindDesc,
        body: Body,
        seam: SeamDesc,
        length: impl Into<Scalar>,
        options: &SeamHingeOptions,
        label: &str,
    ) -> Result<Body, ScriptError> {
        self.seam_hinge_owned(
            kind,
            HingeComponentDesc::LeafA,
            Some(body),
            seam,
            length,
            options,
            label,
        )
    }

    /// Record seam hinge component `component` as the owner of its own
    /// parameters: what an assembly needs, where each part's document
    /// holds one component (PK-712). Leaves take their `body`; the pin
    /// takes none (a new body).
    ///
    /// # Errors
    ///
    /// [`ScriptError::Argument`] for a leaf without a body, a pin with
    /// one, or a print-in-place pin; binding errors.
    #[allow(clippy::too_many_arguments)]
    pub fn seam_hinge_owned(
        &mut self,
        kind: HingeKindDesc,
        component: HingeComponentDesc,
        body: Option<Body>,
        seam: SeamDesc,
        length: impl Into<Scalar>,
        options: &SeamHingeOptions,
        label: &str,
    ) -> Result<Body, ScriptError> {
        if (component == HingeComponentDesc::Pin) != body.is_none()
            || (component == HingeComponentDesc::Pin && kind == HingeKindDesc::PrintInPlace)
        {
            return Err(ScriptError::Argument {
                argument: "component".into(),
                reason: "the leaves grow on a body; a butt hinge's pin is a new body".into(),
            });
        }
        let unit = self.unit;
        let o = options;
        let print = kind == HingeKindDesc::PrintInPlace;
        let length = self.bind("length", &length.into(), unit)?;
        let factor = if print { 1.0 / 200.0 } else { 1.0 / 400.0 };
        let spec = SeamHingeDesc {
            length,
            knuckles: self.bind(
                "knuckles",
                &o.knuckles.clone().unwrap_or(Scalar::Literal(5.0)),
                Unit::Ratio,
            )?,
            pin_diameter: self.bind_or_follow(
                "pin_diameter",
                o.pin.as_ref(),
                length,
                1.0 / 15.0,
                unit,
            )?,
            wall: self.bind_or_follow("wall", o.wall.as_ref(), length, 1.0 / 30.0, unit)?,
            knuckle_diameter: o
                .knuckle_diameter
                .as_ref()
                .map(|s| self.bind("knuckle_diameter", s, unit))
                .transpose()?,
            clearance: self.bind_or_follow(
                "clearance",
                o.clearance.as_ref(),
                length,
                factor,
                unit,
            )?,
            gap: self.bind_or_follow("gap", o.gap.as_ref(), length, factor, unit)?,
        };
        let feature = Feature::SeamHinge {
            kind,
            component,
            seam,
            spec,
            owner: true,
        };
        self.record_seam(feature, body, label)
    }

    /// Record another component of the seam hinge whose owning feature is
    /// labelled `owner`: body B's knuckles on `body`, or (`body` none)
    /// the pin. Labels `<owner>_b` / `<owner>_pin` unless `label` is
    /// given.
    ///
    /// # Errors
    ///
    /// [`ScriptError::Argument`] when `owner` labels no seam hinge, the
    /// component already exists, a leaf has no body or the pin has one,
    /// or a print-in-place hinge is asked for a pin.
    pub fn seam_hinge_component(
        &mut self,
        owner: &str,
        component: HingeComponentDesc,
        body: Option<Body>,
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
            if let Feature::SeamHinge {
                kind,
                component: c,
                seam,
                spec,
                owner: is_owner,
            } = &node.feature
            {
                if *is_owner && node.label == owner {
                    found = Some((*kind, *seam, spec.clone()));
                }
                taken.push((spec.length, *c));
            }
        }
        let Some((kind, seam, spec)) = found else {
            return Err(bad(format!(
                "no seam hinge labelled '{owner}' owns its parameters"
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
        if (component == HingeComponentDesc::Pin) != body.is_none() {
            return Err(bad(
                "the leaves grow on a body; the pin is a new body".into()
            ));
        }
        let label = if label.is_empty() {
            match component {
                HingeComponentDesc::LeafB => format!("{owner}_b"),
                other => format!("{owner}_{}", other.spelling()),
            }
        } else {
            label.to_string()
        };
        let feature = Feature::SeamHinge {
            kind,
            component,
            seam,
            spec,
            owner: false,
        };
        self.record_seam(feature, body, &label)
    }

    fn record_seam(
        &mut self,
        feature: Feature,
        body: Option<Body>,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let id = record_feature(
            &mut self.tree,
            self.tolerance,
            feature,
            label,
            body.map(|b| b.id),
            Vec::new(),
            &mut self.geometry_current,
        )?;
        self.note_source(id);
        let node = self
            .tree
            .feature(id)
            .ok_or(popkern_history::HistoryError::FeatureNotFound)?;
        self.active = Some(node.body);
        Ok(Body::new(node.body))
    }

    /// Record a hinge on bodies `a` and `b` along the seam line through
    /// `origin` in `direction`, `length` long from `origin`: knuckles
    /// grown on both bodies (A's below the seam, B's above) and — for a
    /// butt hinge — a pin body. Labels: `label` (or a minted one) on body
    /// A's feature, `<label>_b`, `<label>_pin`.
    ///
    /// # Errors
    ///
    /// As [`Model::seam_between`] and [`Model::seam_hinge`].
    #[allow(clippy::too_many_arguments)]
    pub fn hinge_bodies(
        &mut self,
        kind: HingeKindDesc,
        a: Body,
        b: Body,
        origin: Point3,
        direction: Vec3,
        length: impl Into<Scalar>,
        options: &SeamHingeOptions,
        label: &str,
    ) -> Result<HingeBodies, ScriptError> {
        let seam = self.seam_between(a, b, origin, direction)?;
        let leaf_a = self.seam_hinge(kind, a, seam, length, options, label)?;
        let owner = self
            .features()
            .last()
            .map(|f| f.label.clone())
            .unwrap_or_default();
        let leaf_b = self.seam_hinge_component(&owner, HingeComponentDesc::LeafB, Some(b), "")?;
        let pin = if kind == HingeKindDesc::Butt {
            Some(self.seam_hinge_component(&owner, HingeComponentDesc::Pin, None, "")?)
        } else {
            None
        };
        Ok(HingeBodies {
            leaf_a: Some(leaf_a),
            leaf_b: Some(leaf_b),
            pin,
        })
    }
}
