//! Public STEP export API and orchestration.
//!
//! The main entry points are [`export_step`] (to bytes) and [`export_step_file`]
//! (to a filesystem path). Both take a [`StepDocument`] and produce valid
//! Part 21 output.

use std::collections::HashMap;
use std::path::Path;

use popkern_core::geometry::Point3;

use crate::error::IoError;
use crate::step::document::StepDocument;
use crate::step::parser::types::{EntityId, EntityInstance, Parameter};
use crate::step::serializer::serialize_part21;

/// Options for STEP export.
#[derive(Debug, Clone)]
pub struct StepExportOptions {
    /// STEP schema to declare in the header.
    pub schema: StepSchema,
    /// Maximum line width for Part 21 output (default 80).
    pub line_width: usize,
    /// Preserve original entity IDs from import (default true).
    pub preserve_entity_ids: bool,
}

impl Default for StepExportOptions {
    fn default() -> Self {
        Self {
            schema: StepSchema::Ap214,
            line_width: 80,
            preserve_entity_ids: true,
        }
    }
}

/// STEP application protocol schema.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StepSchema {
    /// AP 203 — Configuration Controlled 3D Design.
    Ap203,
    /// AP 214 — Core Data for Automotive Mechanical Design Processes.
    Ap214,
}

/// Export a [`StepDocument`] to STEP Part 21 bytes.
///
/// Walks all products, exports their BRep solids as STEP entities, builds
/// the product structure, merges opaque-bag entities, and serializes.
///
/// # Errors
///
/// Returns [`IoError::Export`] if the document contains data that cannot be
/// represented in STEP format.
pub fn export_step(doc: &StepDocument, options: &StepExportOptions) -> Result<Vec<u8>, IoError> {
    let mut ctx = ExportContext::new();

    // Stage 1–2: Export geometry and topology from each BRep solid
    let mut shape_item_ids: Vec<Vec<EntityId>> = Vec::new();
    for product in &doc.products {
        let mut product_shape_item_ids = Vec::new();
        for shape in &product.shapes {
            for brep in &shape.solids {
                product_shape_item_ids.extend(ctx.export_brep_items(brep)?);
            }
        }
        shape_item_ids.push(product_shape_item_ids);
    }

    // Stage 3: Export product structure
    ctx.export_product_structure(doc, &shape_item_ids)?;

    // Stage 4: Merge opaque bag entities
    for entity in &doc.opaque_bag.entities {
        ctx.entities.push(entity.clone());
    }

    // Stage 5: Serialize to Part 21 format
    serialize_part21(&ctx.entities, doc, options)
}

/// Export a [`StepDocument`] to a STEP file on disk.
///
/// # Errors
///
/// Returns [`IoError::Io`] if the file cannot be written, or any error from
/// [`export_step`].
pub fn export_step_file(
    doc: &StepDocument,
    path: &Path,
    options: &StepExportOptions,
) -> Result<(), IoError> {
    let data = export_step(doc, options)?;
    std::fs::write(path, data)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Export context — entity allocation and caching
// ---------------------------------------------------------------------------

/// Mutable context for the export pipeline, tracking emitted entities.
pub(crate) struct ExportContext {
    /// All emitted entity instances, in order.
    pub(crate) entities: Vec<EntityInstance>,
    /// Next entity ID to allocate.
    next_id: EntityId,
    /// Cache of already-exported points, keyed by exact float bits.
    pub(crate) point_cache: HashMap<PointKey, EntityId>,
    /// Cache of already-exported directions, keyed by exact float bits.
    pub(crate) direction_cache: HashMap<PointKey, EntityId>,
    /// Cached ID of the shared 2-D parametric representation context.
    pub(crate) param_repr_context_id: Option<EntityId>,
}

/// Key for deduplicating points/directions by exact float bits.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct PointKey {
    x: u64,
    y: u64,
    z: u64,
}

impl PointKey {
    /// Create a key from a [`Point3`].
    pub(crate) fn from_point(p: &Point3) -> Self {
        Self {
            x: p.x.to_bits(),
            y: p.y.to_bits(),
            z: p.z.to_bits(),
        }
    }

    /// Create a key from a [`Vec3`](popkern_core::geometry::Vec3).
    pub(crate) fn from_vec(v: &popkern_core::geometry::Vec3) -> Self {
        Self {
            x: v.x.to_bits(),
            y: v.y.to_bits(),
            z: v.z.to_bits(),
        }
    }
}

impl ExportContext {
    /// Create a new export context.
    pub(crate) fn new() -> Self {
        Self {
            entities: Vec::new(),
            next_id: 1,
            point_cache: HashMap::new(),
            direction_cache: HashMap::new(),
            param_repr_context_id: None,
        }
    }

    /// Allocate the next entity ID.
    pub(crate) fn alloc_id(&mut self) -> EntityId {
        let id = self.next_id;
        self.next_id += 1;
        id
    }

    /// Emit an entity and return its ID.
    pub(crate) fn emit(&mut self, type_name: &str, params: Vec<Parameter>) -> EntityId {
        let id = self.alloc_id();
        self.entities.push(EntityInstance {
            id,
            type_name: type_name.into(),
            params,
            complex_parts: vec![],
        });
        id
    }
}

// ---------------------------------------------------------------------------
// Parameter helper functions
// ---------------------------------------------------------------------------

/// Build a STEP string parameter.
pub(crate) fn str_param(s: &str) -> Parameter {
    Parameter::String(s.into())
}

/// Build a STEP real parameter.
pub(crate) fn real(f: f64) -> Parameter {
    Parameter::Real(f)
}

/// Build a STEP integer parameter.
pub(crate) fn int(i: i64) -> Parameter {
    Parameter::Int(i)
}

/// Build a STEP entity reference parameter.
pub(crate) fn ref_(id: EntityId) -> Parameter {
    Parameter::EntityRef(id)
}

/// Build a STEP list parameter.
pub(crate) fn list(items: Vec<Parameter>) -> Parameter {
    Parameter::List(items)
}

/// Build a STEP boolean-as-enum parameter (.T. or .F.).
pub(crate) fn bool_(b: bool) -> Parameter {
    Parameter::Enum(if b { "T" } else { "F" }.into())
}

/// Build a STEP enumeration parameter.
pub(crate) fn enum_(s: &str) -> Parameter {
    Parameter::Enum(s.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step::document::*;
    use std::sync::Arc;

    #[test]
    fn default_options() {
        let opts = StepExportOptions::default();
        assert_eq!(opts.schema, StepSchema::Ap214);
        assert_eq!(opts.line_width, 80);
        assert!(opts.preserve_entity_ids);
    }

    #[test]
    fn export_empty_document() {
        let doc = StepDocument {
            products: vec![],
            roots: vec![],
            units: UnitContext::default(),
            opaque_bag: OpaqueBag::default(),
        };
        let opts = StepExportOptions::default();
        let result = export_step(&doc, &opts);
        assert!(result.is_ok());
        let text = String::from_utf8(result.unwrap()).unwrap();
        assert!(text.contains("ISO-10303-21;"));
        assert!(text.contains("END-ISO-10303-21;"));
    }

    #[test]
    fn export_document_with_product_no_shapes() {
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
        let opts = StepExportOptions::default();
        let result = export_step(&doc, &opts);
        assert!(result.is_ok());
        let text = String::from_utf8(result.unwrap()).unwrap();
        assert!(text.contains("PRODUCT("));
    }

    #[test]
    fn export_with_brep_round_trip_structure() {
        use popkern_core::tolerance::Tolerance;
        use popkern_core::topology::BRep;

        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let doc = StepDocument {
            products: vec![StepProductNode {
                name: "Cube".into(),
                id: "C001".into(),
                children: vec![],
                shapes: vec![Arc::new(ShapeRep {
                    kind: ShapeRepKind::AdvancedBrep,
                    solids: vec![brep],
                    context: RepresentationContext::default(),
                    non_brep_items: vec![],
                })],
                metadata: ProductMetadata::default(),
            }],
            roots: vec![0],
            units: UnitContext::default(),
            opaque_bag: OpaqueBag::default(),
        };

        let opts = StepExportOptions::default();
        let result = export_step(&doc, &opts);
        assert!(result.is_ok(), "export failed: {:?}", result.err());

        let text = String::from_utf8(result.unwrap()).unwrap();

        // Verify basic STEP structure
        assert!(text.contains("ISO-10303-21;"));
        assert!(text.contains("HEADER;"));
        assert!(text.contains("DATA;"));
        assert!(text.contains("CARTESIAN_POINT"));
        assert!(text.contains("VERTEX_POINT"));
        assert!(text.contains("EDGE_CURVE"));
        assert!(text.contains("ADVANCED_FACE"));
        assert!(text.contains("CLOSED_SHELL"));
        assert!(text.contains("MANIFOLD_SOLID_BREP"));
        assert!(text.contains("PRODUCT("));
        assert!(text.contains("END-ISO-10303-21;"));
    }

    #[test]
    fn export_step_file_writes_to_disk() {
        let doc = StepDocument {
            products: vec![],
            roots: vec![],
            units: UnitContext::default(),
            opaque_bag: OpaqueBag::default(),
        };
        let dir = std::env::temp_dir().join("popkern_test_export");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("test_export.stp");

        let opts = StepExportOptions::default();
        export_step_file(&doc, &path, &opts).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("ISO-10303-21;"));

        // Cleanup
        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_dir(&dir);
    }

    #[test]
    fn parameter_helpers() {
        assert_eq!(str_param("hello"), Parameter::String("hello".into()));
        assert_eq!(real(1.5), Parameter::Real(1.5));
        assert_eq!(int(42), Parameter::Int(42));
        assert_eq!(ref_(7), Parameter::EntityRef(7));
        assert_eq!(bool_(true), Parameter::Enum("T".into()));
        assert_eq!(bool_(false), Parameter::Enum("F".into()));
        assert_eq!(enum_("BOTH"), Parameter::Enum("BOTH".into()));
        assert_eq!(
            list(vec![Parameter::Int(1)]),
            Parameter::List(vec![Parameter::Int(1)])
        );
    }

    #[test]
    fn point_key_deduplication() {
        let p1 = Point3::new(1.0, 2.0, 3.0);
        let p2 = Point3::new(1.0, 2.0, 3.0);
        let p3 = Point3::new(1.0, 2.0, 3.1);
        assert_eq!(PointKey::from_point(&p1), PointKey::from_point(&p2));
        assert_ne!(PointKey::from_point(&p1), PointKey::from_point(&p3));
    }
}
