//! Primitive recorder verbs (PK-1083): `cuboid` (Rhai `box`),
//! `cylinder`, `cone`, `sphere`, `prism`, `pyramid`, `wedge` and
//! `torus` on [`Model`], each a new
//! body recording `Feature::Primitive`. Lengths bind in the model unit,
//! angles in degrees, so every dimension is a parameter.

use popkern_history::feature::Feature;
use popkern_history::feature_primitive::{PrimitiveDesc, PrimitiveFrame, PrimitiveShape};
use popkern_history::parameter::ParameterId;
use popkern_history::unit::Unit;

use crate::body::Body;
use crate::error::ScriptError;
use crate::model::Model;
use crate::record::record_feature;
use crate::scalar::{Scalar, bind_scalar};

/// Where a primitive stands: its local origin, axis (local z) and
/// reference direction (local x). Defaults to the world frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Placement {
    /// The local origin, in the model unit.
    pub at: [f64; 3],
    /// The local z.
    pub axis: [f64; 3],
    /// The local x (its part along `axis` is dropped).
    pub reference: [f64; 3],
}

impl Default for Placement {
    fn default() -> Self {
        let world = PrimitiveFrame::default();
        Placement {
            at: world.origin,
            axis: world.axis,
            reference: world.reference,
        }
    }
}

/// A regular polygon's side count and how its radius is measured.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RegularPolygon {
    /// Number of sides (at least 3).
    pub sides: u32,
    /// The radius runs to the middle of a side (the inscribed circle)
    /// rather than to a corner.
    pub flats: bool,
}

impl Model {
    fn bind_length(&mut self, role: &str, s: Scalar) -> Result<ParameterId, ScriptError> {
        let unit = self.unit;
        bind_scalar(&mut self.tree, &mut self.pending_counter, role, &s, unit)
    }

    fn bind_degrees(&mut self, role: &str, s: Scalar) -> Result<ParameterId, ScriptError> {
        bind_scalar(
            &mut self.tree,
            &mut self.pending_counter,
            role,
            &s,
            Unit::Degrees,
        )
    }

    fn primitive(
        &mut self,
        shape: PrimitiveShape,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let si = |v: [f64; 3], unit: Unit| v.map(|c| unit.to_si(c));
        let frame = PrimitiveFrame {
            origin: si(placement.at, self.unit),
            axis: placement.axis,
            reference: placement.reference,
        };
        let feature = Feature::Primitive(PrimitiveDesc { shape, frame });
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

    /// A box `x` × `y` × `z`: a corner at the placement's origin, or
    /// centred on it.
    pub fn cuboid(
        &mut self,
        size: [Scalar; 3],
        centred: bool,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let [x, y, z] = size;
        let shape = PrimitiveShape::Box {
            x: self.bind_length("x", x)?,
            y: self.bind_length("y", y)?,
            z: self.bind_length("z", z)?,
            centred,
        };
        self.primitive(shape, placement, label)
    }

    fn bind_angle(&mut self, angle: Option<Scalar>) -> Result<Option<ParameterId>, ScriptError> {
        angle.map(|a| self.bind_degrees("angle", a)).transpose()
    }

    /// A cylinder standing on the placement's origin along its axis,
    /// optionally a sector of `angle` degrees from the reference.
    pub fn cylinder(
        &mut self,
        radius: impl Into<Scalar>,
        height: impl Into<Scalar>,
        angle: Option<Scalar>,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let shape = PrimitiveShape::Cylinder {
            radius: self.bind_length("radius", radius.into())?,
            height: self.bind_length("height", height.into())?,
            angle: self.bind_angle(angle)?,
        };
        self.primitive(shape, placement, label)
    }

    /// A cone or frustum: base radius `r1` on the origin, top radius
    /// `r2` at `height` (either may be 0), optionally a sector.
    pub fn cone(
        &mut self,
        r1: impl Into<Scalar>,
        r2: impl Into<Scalar>,
        height: impl Into<Scalar>,
        angle: Option<Scalar>,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let shape = PrimitiveShape::Cone {
            r1: self.bind_length("r1", r1.into())?,
            r2: self.bind_length("r2", r2.into())?,
            height: self.bind_length("height", height.into())?,
            angle: self.bind_angle(angle)?,
        };
        self.primitive(shape, placement, label)
    }

    /// A sphere centred on the origin, optionally only between the
    /// latitudes `[lo, hi]` (degrees from the equator) and only a
    /// sector of `angle` degrees.
    pub fn sphere(
        &mut self,
        radius: impl Into<Scalar>,
        latitudes: Option<[Scalar; 2]>,
        angle: Option<Scalar>,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let radius = self.bind_length("radius", radius.into())?;
        let latitudes = match latitudes {
            Some([lo, hi]) => Some([
                self.bind_degrees("lat_lo", lo)?,
                self.bind_degrees("lat_hi", hi)?,
            ]),
            None => None,
        };
        let angle = self.bind_angle(angle)?;
        self.primitive(
            PrimitiveShape::Sphere {
                radius,
                latitudes,
                angle,
            },
            placement,
            label,
        )
    }

    /// A regular prism standing on the origin, corner 0 on the
    /// reference direction.
    pub fn prism(
        &mut self,
        polygon: RegularPolygon,
        radius: impl Into<Scalar>,
        height: impl Into<Scalar>,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let shape = PrimitiveShape::Prism {
            sides: polygon.sides,
            radius: self.bind_length("radius", radius.into())?,
            height: self.bind_length("height", height.into())?,
            flats: polygon.flats,
        };
        self.primitive(shape, placement, label)
    }

    /// A regular pyramid: base radius `r1` on the origin, top radius
    /// `r2` at `height` (either may be 0, an apex).
    pub fn pyramid(
        &mut self,
        polygon: RegularPolygon,
        radii: [Scalar; 2],
        height: impl Into<Scalar>,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let [r1, r2] = radii;
        let shape = PrimitiveShape::Pyramid {
            sides: polygon.sides,
            r1: self.bind_length("r1", r1)?,
            r2: self.bind_length("r2", r2)?,
            height: self.bind_length("height", height.into())?,
            flats: polygon.flats,
        };
        self.primitive(shape, placement, label)
    }

    /// A wedge: the box `x` × `y` × `z` on a corner at the origin
    /// whose top runs only `top_x` along x (0: a ramp).
    pub fn wedge(
        &mut self,
        size: [Scalar; 3],
        top_x: impl Into<Scalar>,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let [x, y, z] = size;
        let shape = PrimitiveShape::Wedge {
            x: self.bind_length("x", x)?,
            y: self.bind_length("y", y)?,
            z: self.bind_length("z", z)?,
            top_x: self.bind_length("top_x", top_x.into())?,
        };
        self.primitive(shape, placement, label)
    }

    /// A torus about the axis, centred on the origin, optionally a
    /// sector of `angle` degrees from the reference direction.
    pub fn torus(
        &mut self,
        major: impl Into<Scalar>,
        minor: impl Into<Scalar>,
        angle: Option<Scalar>,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let major = self.bind_length("major", major.into())?;
        let minor = self.bind_length("minor", minor.into())?;
        let angle = self.bind_angle(angle)?;
        self.primitive(
            PrimitiveShape::Torus {
                major,
                minor,
                angle,
            },
            placement,
            label,
        )
    }
}
