//! Document model for STEP import/export.
//!
//! [`StepDocument`] is the top-level result of importing a STEP file and the
//! starting point for exporting. It captures the product assembly tree, BRep
//! geometry, units, presentation metadata, and any unrecognized entities for
//! round-trip fidelity.

use std::collections::HashMap;
use std::sync::Arc;

use popkern_core::geometry::Transform3;
use popkern_core::topology::BRep;

use crate::step::parser::types::EntityInstance;

/// The top-level result of importing a STEP file.
#[derive(Debug)]
pub struct StepDocument {
    /// All products found in the file, in discovery order.
    pub products: Vec<StepProductNode>,
    /// Indices into `products` for the top-level (root) products.
    pub roots: Vec<usize>,
    /// Unit context extracted from the file.
    pub units: UnitContext,
    /// Unrecognized entities preserved for round-trip export.
    pub opaque_bag: OpaqueBag,
}

/// A single product node in the assembly tree.
#[derive(Debug, Clone)]
pub struct StepProductNode {
    /// Human-readable product name.
    pub name: String,
    /// Product identifier string.
    pub id: String,
    /// Assembly children (sub-components with placement transforms).
    pub children: Vec<AssemblyChild>,
    /// Shape representations owned by this product.
    pub shapes: Vec<Arc<ShapeRep>>,
    /// Presentation metadata (colors, layers, custom properties).
    pub metadata: ProductMetadata,
}

/// An assembly child with a placement transform.
#[derive(Debug, Clone)]
pub struct AssemblyChild {
    /// Index into [`StepDocument::products`].
    pub product_index: usize,
    /// Transform placing the child relative to the parent.
    pub transform: Transform3,
    /// Optional instance name.
    pub instance_name: Option<String>,
}

/// A shape representation containing BRep solids.
#[derive(Debug, Clone)]
pub struct ShapeRep {
    /// Kind of representation.
    pub kind: ShapeRepKind,
    /// BRep solids extracted from this representation.
    pub solids: Vec<BRep>,
    /// The representation context (units, dimension).
    pub context: RepresentationContext,
    /// Non-BRep items preserved for round-trip.
    pub non_brep_items: Vec<EntityInstance>,
}

/// The kind of shape representation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ShapeRepKind {
    /// ADVANCED_BREP_SHAPE_REPRESENTATION
    AdvancedBrep,
    /// FACETED_BREP_SHAPE_REPRESENTATION
    FacetedBrep,
    /// Some other representation type.
    Other(String),
}

/// Unit and context information for a representation.
#[derive(Debug, Clone)]
pub struct RepresentationContext {
    /// Coordinate space dimension (typically 3).
    pub dimension: u32,
    /// Length unit.
    pub length_unit: LengthUnit,
    /// Angle unit.
    pub angle_unit: AngleUnit,
    /// Geometric uncertainty (tolerance for coincidence).
    pub uncertainty: f64,
}

impl Default for RepresentationContext {
    fn default() -> Self {
        Self {
            dimension: 3,
            length_unit: LengthUnit::Millimetre,
            angle_unit: AngleUnit::Radian,
            uncertainty: 1e-6,
        }
    }
}

/// Physical length unit.
#[derive(Debug, Clone, PartialEq)]
pub enum LengthUnit {
    /// SI metre.
    Metre,
    /// SI millimetre.
    Millimetre,
    /// Imperial inch.
    Inch,
    /// Some other unit.
    Other(String),
}

/// Angle unit.
#[derive(Debug, Clone, PartialEq)]
pub enum AngleUnit {
    /// Radians.
    Radian,
    /// Degrees.
    Degree,
    /// Some other unit.
    Other(String),
}

/// File-level unit context extracted from the STEP header/data.
#[derive(Debug, Clone, Default)]
pub struct UnitContext {
    /// The representation context describing units.
    pub context: RepresentationContext,
}

/// Presentation metadata for a product (colors, layers, custom props).
#[derive(Debug, Clone, Default)]
pub struct ProductMetadata {
    /// Solid-level color (RGB, each in \[0,1\]).
    pub solid_color: Option<[f32; 3]>,
    /// Per-face colors keyed by face index within the solid.
    pub face_colors: HashMap<usize, [f32; 3]>,
    /// Layer assignment.
    pub layer: Option<String>,
    /// Arbitrary string properties.
    pub custom: HashMap<String, String>,
    /// Import only: the file's solid representation items (entity ids)
    /// this product took, so a `STYLED_ITEM` colours the product that
    /// owns its target and no other.
    pub solid_items: Vec<u64>,
}

/// Stores unrecognized entities for round-trip preservation.
#[derive(Debug, Clone, Default)]
pub struct OpaqueBag {
    /// Raw entity instances not consumed by the import pipeline.
    pub entities: Vec<EntityInstance>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn construct_step_document() {
        let doc = StepDocument {
            products: vec![StepProductNode {
                name: "Widget".into(),
                id: "W001".into(),
                children: vec![],
                shapes: vec![],
                metadata: ProductMetadata::default(),
            }],
            roots: vec![0],
            units: UnitContext::default(),
            opaque_bag: OpaqueBag::default(),
        };
        assert_eq!(doc.products.len(), 1);
        assert_eq!(doc.products[0].name, "Widget");
        assert_eq!(doc.roots, vec![0]);
    }

    #[test]
    fn default_representation_context() {
        let ctx = RepresentationContext::default();
        assert_eq!(ctx.dimension, 3);
        assert_eq!(ctx.length_unit, LengthUnit::Millimetre);
        assert_eq!(ctx.angle_unit, AngleUnit::Radian);
        assert!((ctx.uncertainty - 1e-6).abs() < f64::EPSILON);
    }

    #[test]
    fn shape_rep_kind_equality() {
        assert_eq!(ShapeRepKind::AdvancedBrep, ShapeRepKind::AdvancedBrep);
        assert_ne!(ShapeRepKind::AdvancedBrep, ShapeRepKind::FacetedBrep);
        assert_eq!(
            ShapeRepKind::Other("foo".into()),
            ShapeRepKind::Other("foo".into())
        );
    }

    #[test]
    fn assembly_child_with_transform() {
        let child = AssemblyChild {
            product_index: 1,
            transform: Transform3::identity(),
            instance_name: Some("bolt_1".into()),
        };
        assert_eq!(child.product_index, 1);
        assert_eq!(child.instance_name.as_deref(), Some("bolt_1"));
    }
}
