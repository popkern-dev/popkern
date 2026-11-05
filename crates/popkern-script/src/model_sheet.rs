//! Sheet recorder verbs (PK-1083): `sheet_rect`, `sheet_disc` and
//! `sheet_polygon` on [`Model`], each a one-face sheet body with free
//! edges recording `Feature::Primitive`, and `sheet_fill`, the planar
//! face a closed wire bounds (the wire referenced, not consumed).

use popkern_history::feature_primitive::PrimitiveShape;

use crate::body::Body;
use crate::error::ScriptError;
use crate::model::Model;
use crate::model_primitive::Placement;
use crate::scalar::Scalar;
use crate::selection::EdgeSelector;

impl Model {
    /// The rectangle `x` × `y` in the placement's xy plane, a corner
    /// on its origin (or centred on it), facing its axis.
    pub fn sheet_rect(
        &mut self,
        size: [Scalar; 2],
        centred: bool,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let [x, y] = size;
        let shape = PrimitiveShape::SheetRect {
            x: self.bind_length("x", x)?,
            y: self.bind_length("y", y)?,
            centred,
        };
        self.primitive(shape, placement, label)
    }

    /// The disc of `radius` about the placement's origin.
    pub fn sheet_disc(
        &mut self,
        radius: impl Into<Scalar>,
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let shape = PrimitiveShape::SheetDisc {
            radius: self.bind_length("radius", radius.into())?,
        };
        self.primitive(shape, placement, label)
    }

    /// The polygon through `points` (the placement's xy plane, model
    /// unit).
    pub fn sheet_polygon(
        &mut self,
        points: &[[f64; 2]],
        placement: Placement,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let unit = self.unit;
        let shape = PrimitiveShape::SheetPolygon {
            points: points.iter().map(|p| p.map(|c| unit.to_si(c))).collect(),
        };
        self.primitive(shape, placement, label)
    }

    /// The planar face bounded by the `edges` of `wire` (a closed
    /// chain, any order and direction).
    pub fn sheet_fill(
        &mut self,
        wire: Body,
        edges: &EdgeSelector,
        label: &str,
    ) -> Result<Body, ScriptError> {
        let snapshot = wire.edges(self, edges)?;
        let record = Body::selection_record("edges", &snapshot);
        if snapshot.names.is_empty() {
            return Err(ScriptError::EmptySelection {
                selection: snapshot.describe,
                kind: crate::error::SelectionKind::Edges,
                available: 0,
            });
        }
        let shape = PrimitiveShape::SheetFill {
            edges: snapshot.names,
        };
        self.primitive_selecting(shape, Placement::default(), label, vec![record])
    }
}
