//! Wire recorder verbs (PK-1083): `wire_line`, `wire_polyline`,
//! `wire_arc`, `wire_circle`, `wire_helix` and `wire_spline` on
//! [`Model`], each a new wire body (edges on no face) recording
//! `Feature::Primitive`. A sweep or pipe takes its edge as a path.
//! Points are literals in the model unit (stored in metres), placed by
//! the frame like every primitive; radii, pitch and turns are
//! parameters.

use popkern_history::feature_primitive::PrimitiveShape;

use crate::body::Body;
use crate::error::ScriptError;
use crate::model::Model;
use crate::model_primitive::Placement;
use crate::scalar::Scalar;

impl Model {
    fn si_points(&self, points: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let unit = self.unit;
        points.iter().map(|p| p.map(|c| unit.to_si(c))).collect()
    }

    /// The straight wire from `from` to `to`.
    pub fn wire_line(
        &mut self,
        from: [f64; 3],
        to: [f64; 3],
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        self.wire_polyline(&[from, to], false, placement, label)
    }

    /// Straight edges through `points`, joined back to the first when
    /// `closed`.
    pub fn wire_polyline(
        &mut self,
        points: &[[f64; 3]],
        closed: bool,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let shape = PrimitiveShape::WirePolyline {
            points: self.si_points(points),
            closed,
        };
        self.primitive(shape, placement, label)
    }

    /// The circular arc from `start` through `via` to `end`.
    pub fn wire_arc(
        &mut self,
        [start, via, end]: [[f64; 3]; 3],
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let p = self.si_points(&[start, via, end]);
        let shape = PrimitiveShape::WireArc {
            start: p[0],
            via: p[1],
            end: p[2],
        };
        self.primitive(shape, placement, label)
    }

    /// The circle of `radius` about the placement's axis, starting on
    /// its reference direction.
    pub fn wire_circle(
        &mut self,
        radius: impl Into<Scalar>,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let shape = PrimitiveShape::WireCircle {
            radius: self.bind_length("radius", radius.into())?,
        };
        self.primitive(shape, placement, label)
    }

    /// A helix of `radius` about the placement's axis from its
    /// reference direction, rising `pitch` per turn for `turns`.
    pub fn wire_helix(
        &mut self,
        radius: impl Into<Scalar>,
        pitch: impl Into<Scalar>,
        turns: impl Into<Scalar>,
        left_handed: bool,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let radius = self.bind_length("radius", radius.into())?;
        let pitch = self.bind_length("pitch", pitch.into())?;
        let turns = crate::scalar::bind_scalar(
            &mut self.tree,
            &mut self.pending_counter,
            "turns",
            &turns.into(),
            popkern_history::unit::Unit::Ratio,
        )?;
        let shape = PrimitiveShape::WireHelix {
            radius,
            pitch,
            turns,
            left_handed,
        };
        self.primitive(shape, placement, label)
    }

    /// The spline of `degree` interpolating `points`, or with
    /// `control` the B-spline over them as its control polygon.
    pub fn wire_spline(
        &mut self,
        points: &[[f64; 3]],
        degree: u32,
        control: bool,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let shape = PrimitiveShape::WireSpline {
            points: self.si_points(points),
            degree,
            control,
        };
        self.primitive(shape, placement, label)
    }
}
