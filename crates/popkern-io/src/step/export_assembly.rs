//! Assembly occurrences in STEP: one `NEXT_ASSEMBLY_USAGE_OCCURRENCE`
//! per placed child, with the shape wiring a reader needs to place it
//! (PK-629):
//!
//! ```text
//! NAUO(id, instance name, '', parent PD, child PD, $)
//! PRODUCT_DEFINITION_SHAPE('', '', NAUO)
//! CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(rel, that PDS)
//! rel = ( REPRESENTATION_RELATIONSHIP('', '', child rep, parent rep)
//!         REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(idt)
//!         SHAPE_REPRESENTATION_RELATIONSHIP() )
//! idt = ITEM_DEFINED_TRANSFORMATION('', '', child origin axis, placement axis in parent)
//! ```
//!
//! The placement maps the child's origin frame onto the placement
//! frame: `T = placement ∘ origin⁻¹`, which is `placement` itself
//! because the origin axis is the identity frame. The reader
//! (`mapping::product`) inverts exactly this.

use popkern_core::geometry::{Point3, Transform3, Vec3};

use crate::step::document::StepDocument;
use crate::step::export::{ExportContext, ref_, str_param};
use crate::step::parser::types::{EntityId, EntityInstance};

/// Per-product representation handles the relationship pass wires
/// together, parallel to `StepDocument::products`.
#[derive(Debug, Default)]
pub(crate) struct AssemblyStructure {
    /// The product's shape representation, when it has one.
    pub(crate) reps: Vec<Option<EntityId>>,
    /// The origin axis in the product's representation, when the
    /// product is placed under a parent.
    pub(crate) origins: Vec<Option<EntityId>>,
    /// One placement axis per child, in the parent's representation.
    pub(crate) child_axes: Vec<Vec<EntityId>>,
}

impl ExportContext {
    /// An `AXIS2_PLACEMENT_3D` for the rigid frame of `transform`
    /// (origin = translation column, axis = third column, ref
    /// direction = first column).
    pub(crate) fn export_transform_axis(&mut self, transform: &Transform3) -> EntityId {
        let m = &transform.matrix;
        let origin = Point3::new(m[0][3], m[1][3], m[2][3]);
        let z = Vec3::new(m[0][2], m[1][2], m[2][2]);
        let x = Vec3::new(m[0][0], m[1][0], m[2][0]);
        self.export_placement(&origin, &z, &x)
    }

    /// Emit the occurrence and its placement wiring for every child of
    /// every product.
    pub(crate) fn export_assembly_relationships(
        &mut self,
        doc: &StepDocument,
        product_def_ids: &[EntityId],
        structure: &AssemblyStructure,
    ) {
        let mut occurrence = 0usize;
        for (i, product) in doc.products.iter().enumerate() {
            let Some(&parent_pd) = product_def_ids.get(i) else {
                continue;
            };
            for (k, child) in product.children.iter().enumerate() {
                let Some(&child_pd) = product_def_ids.get(child.product_index) else {
                    continue;
                };
                occurrence += 1;
                let instance_name = child.instance_name.as_deref().unwrap_or("");
                let nauo_id = self.emit(
                    "NEXT_ASSEMBLY_USAGE_OCCURRENCE",
                    vec![
                        str_param(&occurrence.to_string()),
                        str_param(instance_name),
                        str_param(""),
                        ref_(parent_pd),
                        ref_(child_pd),
                        crate::step::parser::types::Parameter::Omitted,
                    ],
                );
                let child_rep = structure.reps.get(child.product_index).copied().flatten();
                let parent_rep = structure.reps.get(i).copied().flatten();
                let child_origin = structure
                    .origins
                    .get(child.product_index)
                    .copied()
                    .flatten();
                let placement = structure.child_axes.get(i).and_then(|a| a.get(k)).copied();
                let (Some(child_rep), Some(parent_rep), Some(child_origin), Some(placement)) =
                    (child_rep, parent_rep, child_origin, placement)
                else {
                    // Unreachable for documents built by the product
                    // pass (every child gets an origin and every parent
                    // a placement axis); the NAUO alone still names the
                    // occurrence.
                    continue;
                };
                let idt = self.emit(
                    "ITEM_DEFINED_TRANSFORMATION",
                    vec![
                        str_param(""),
                        str_param(""),
                        ref_(child_origin),
                        ref_(placement),
                    ],
                );
                let rel = self.alloc_id();
                // Complex instance: partials in alphabetical order.
                self.entities.push(EntityInstance {
                    id: rel,
                    type_name: "REPRESENTATION_RELATIONSHIP".into(),
                    params: vec![
                        str_param(""),
                        str_param(""),
                        ref_(child_rep),
                        ref_(parent_rep),
                    ],
                    complex_parts: vec![
                        (
                            "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION".into(),
                            vec![ref_(idt)],
                        ),
                        ("SHAPE_REPRESENTATION_RELATIONSHIP".into(), vec![]),
                    ],
                });
                let pds = self.emit(
                    "PRODUCT_DEFINITION_SHAPE",
                    vec![str_param(""), str_param(""), ref_(nauo_id)],
                );
                self.emit(
                    "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION",
                    vec![ref_(rel), ref_(pds)],
                );
            }
        }
    }
}
