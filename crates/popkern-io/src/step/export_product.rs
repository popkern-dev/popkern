//! Export product structure to STEP entities.
//!
//! Emits PRODUCT, PRODUCT_DEFINITION_FORMATION, PRODUCT_DEFINITION,
//! PRODUCT_DEFINITION_SHAPE, SHAPE_DEFINITION_REPRESENTATION,
//! ADVANCED_BREP_SHAPE_REPRESENTATION, and assembly occurrence entities.

use popkern_core::geometry::Transform3;

use crate::error::IoError;
use crate::step::document::{
    AngleUnit, LengthUnit, RepresentationContext, ShapeRepKind, StepDocument,
};
use crate::step::export::{ExportContext, enum_, int, list, real, ref_, str_param};
use crate::step::parser::types::EntityId;

impl ExportContext {
    /// Export the full product structure from a [`StepDocument`].
    ///
    /// `solid_ids` is parallel to `doc.products` — each entry lists the
    /// entity IDs of the MANIFOLD_SOLID_BREP entities for that product's
    /// shapes.
    pub(crate) fn export_product_structure(
        &mut self,
        doc: &StepDocument,
        solid_ids: &[Vec<EntityId>],
    ) -> Result<(), IoError> {
        // Emit shared context entities
        let ctx_id = self.export_representation_context(&doc.units.context);

        // Emit each product
        let mut product_def_ids: Vec<EntityId> = Vec::new();
        for (i, product) in doc.products.iter().enumerate() {
            let solids = solid_ids.get(i).map(|v| v.as_slice()).unwrap_or(&[]);

            // PRODUCT
            let prod_id = self.emit(
                "PRODUCT",
                vec![
                    str_param(&product.id),
                    str_param(&product.name),
                    str_param(""),
                    list(vec![]),
                ],
            );

            // PRODUCT_DEFINITION_FORMATION
            let pdf_id = self.emit(
                "PRODUCT_DEFINITION_FORMATION",
                vec![str_param(""), str_param(""), ref_(prod_id)],
            );

            // PRODUCT_DEFINITION
            let pd_id = self.emit(
                "PRODUCT_DEFINITION",
                vec![
                    str_param(""),
                    str_param(""),
                    ref_(pdf_id),
                    // product_definition_context omitted for simplicity
                ],
            );
            product_def_ids.push(pd_id);

            // PRODUCT_DEFINITION_SHAPE
            let pds_id = self.emit(
                "PRODUCT_DEFINITION_SHAPE",
                vec![str_param(""), str_param(""), ref_(pd_id)],
            );

            // Build SHAPE_REPRESENTATION containing the solid items
            if !solids.is_empty() {
                let items: Vec<_> = solids.iter().map(|&id| ref_(id)).collect();

                // Determine the representation type name
                let rep_type = if product
                    .shapes
                    .iter()
                    .any(|s| s.kind == ShapeRepKind::AdvancedBrep)
                {
                    "ADVANCED_BREP_SHAPE_REPRESENTATION"
                } else {
                    "SHAPE_REPRESENTATION"
                };

                let sr_id = self.emit(rep_type, vec![str_param(""), list(items), ref_(ctx_id)]);

                // SHAPE_DEFINITION_REPRESENTATION
                self.emit(
                    "SHAPE_DEFINITION_REPRESENTATION",
                    vec![ref_(pds_id), ref_(sr_id)],
                );
            }
        }

        // Export assembly relationships
        for (i, product) in doc.products.iter().enumerate() {
            let parent_pd = product_def_ids[i];
            for child in &product.children {
                if child.product_index >= product_def_ids.len() {
                    continue;
                }
                let child_pd = product_def_ids[child.product_index];

                let instance_name = child.instance_name.as_deref().unwrap_or("");

                // NEXT_ASSEMBLY_USAGE_OCCURRENCE
                let nauo_id = self.emit(
                    "NEXT_ASSEMBLY_USAGE_OCCURRENCE",
                    vec![
                        str_param(""),
                        str_param(instance_name),
                        str_param(""),
                        ref_(parent_pd),
                        ref_(child_pd),
                    ],
                );

                // If non-identity transform, emit placement
                if !is_identity(&child.transform) {
                    self.export_assembly_transform(nauo_id, &child.transform);
                }
            }
        }

        Ok(())
    }

    /// Export a geometric representation context with units.
    fn export_representation_context(&mut self, ctx: &RepresentationContext) -> EntityId {
        // SI length unit
        let lu_id = match ctx.length_unit {
            LengthUnit::Millimetre => self.emit("SI_UNIT", vec![enum_("MILLI"), enum_("METRE")]),
            LengthUnit::Metre => self.emit("SI_UNIT", vec![enum_("*"), enum_("METRE")]),
            LengthUnit::Inch => {
                // Conversion-based unit
                let si = self.emit("SI_UNIT", vec![enum_("*"), enum_("METRE")]);
                self.emit(
                    "CONVERSION_BASED_UNIT_AND_LENGTH_UNIT",
                    vec![str_param("inch"), ref_(si), real(0.0254)],
                )
            }
            LengthUnit::Other(ref name) => {
                self.emit("SI_UNIT", vec![str_param(name), enum_("METRE")])
            }
        };

        // SI angle unit
        let au_id = match ctx.angle_unit {
            AngleUnit::Radian => self.emit("SI_UNIT", vec![enum_("*"), enum_("RADIAN")]),
            AngleUnit::Degree => {
                let si = self.emit("SI_UNIT", vec![enum_("*"), enum_("RADIAN")]);
                self.emit(
                    "CONVERSION_BASED_UNIT_AND_PLANE_ANGLE_UNIT",
                    vec![
                        str_param("degree"),
                        ref_(si),
                        real(std::f64::consts::PI / 180.0),
                    ],
                )
            }
            AngleUnit::Other(ref name) => {
                self.emit("SI_UNIT", vec![str_param(name), enum_("RADIAN")])
            }
        };

        // Uncertainty measure
        let unc_id = self.emit(
            "UNCERTAINTY_MEASURE_WITH_UNIT",
            vec![
                real(ctx.uncertainty),
                ref_(lu_id),
                str_param("distance_accuracy_value"),
            ],
        );

        // GEOMETRIC_REPRESENTATION_CONTEXT
        self.emit(
            "GEOMETRIC_REPRESENTATION_CONTEXT",
            vec![
                str_param(""),
                str_param("3D"),
                int(ctx.dimension as i64),
                ref_(lu_id),
                ref_(au_id),
                ref_(unc_id),
            ],
        )
    }

    /// Emit entities for an assembly child placement transform.
    fn export_assembly_transform(&mut self, _nauo_id: EntityId, transform: &Transform3) {
        let m = &transform.matrix;

        // Extract origin (translation column)
        let origin = popkern_core::geometry::Point3::new(m[0][3], m[1][3], m[2][3]);
        // Extract Z axis (3rd column of rotation)
        let z = popkern_core::geometry::Vec3::new(m[0][2], m[1][2], m[2][2]);
        // Extract X axis (1st column of rotation)
        let x = popkern_core::geometry::Vec3::new(m[0][0], m[1][0], m[2][0]);

        self.export_placement(&origin, &z, &x);
    }
}

/// Check if a transform is the identity (no rotation, no translation).
fn is_identity(t: &Transform3) -> bool {
    let id = Transform3::identity();
    for r in 0..4 {
        for c in 0..4 {
            if (t.matrix[r][c] - id.matrix[r][c]).abs() > 1e-15 {
                return false;
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step::document::*;

    #[test]
    fn export_single_product_no_shapes() {
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

        let mut ctx = ExportContext::new();
        ctx.export_product_structure(&doc, &[vec![]]).unwrap();

        let count =
            |name: &str| -> usize { ctx.entities.iter().filter(|e| e.type_name == name).count() };

        assert_eq!(count("PRODUCT"), 1);
        assert_eq!(count("PRODUCT_DEFINITION_FORMATION"), 1);
        assert_eq!(count("PRODUCT_DEFINITION"), 1);
        assert_eq!(count("PRODUCT_DEFINITION_SHAPE"), 1);
    }

    #[test]
    fn export_product_with_solid() {
        let doc = StepDocument {
            products: vec![StepProductNode {
                name: "Box".into(),
                id: "B001".into(),
                children: vec![],
                shapes: vec![std::sync::Arc::new(ShapeRep {
                    kind: ShapeRepKind::AdvancedBrep,
                    solids: vec![],
                    context: RepresentationContext::default(),
                    non_brep_items: vec![],
                })],
                metadata: ProductMetadata::default(),
            }],
            roots: vec![0],
            units: UnitContext::default(),
            opaque_bag: OpaqueBag::default(),
        };

        let mut ctx = ExportContext::new();
        // Provide a fake solid ID
        ctx.export_product_structure(&doc, &[vec![999]]).unwrap();

        let count =
            |name: &str| -> usize { ctx.entities.iter().filter(|e| e.type_name == name).count() };

        assert_eq!(count("ADVANCED_BREP_SHAPE_REPRESENTATION"), 1);
        assert_eq!(count("SHAPE_DEFINITION_REPRESENTATION"), 1);
    }

    #[test]
    fn export_assembly_with_children() {
        let doc = StepDocument {
            products: vec![
                StepProductNode {
                    name: "Assembly".into(),
                    id: "A001".into(),
                    children: vec![AssemblyChild {
                        product_index: 1,
                        transform: Transform3::identity(),
                        instance_name: Some("part_1".into()),
                    }],
                    shapes: vec![],
                    metadata: ProductMetadata::default(),
                },
                StepProductNode {
                    name: "Part".into(),
                    id: "P001".into(),
                    children: vec![],
                    shapes: vec![],
                    metadata: ProductMetadata::default(),
                },
            ],
            roots: vec![0],
            units: UnitContext::default(),
            opaque_bag: OpaqueBag::default(),
        };

        let mut ctx = ExportContext::new();
        ctx.export_product_structure(&doc, &[vec![], vec![]])
            .unwrap();

        let count =
            |name: &str| -> usize { ctx.entities.iter().filter(|e| e.type_name == name).count() };

        assert_eq!(count("PRODUCT"), 2);
        assert_eq!(count("NEXT_ASSEMBLY_USAGE_OCCURRENCE"), 1);
    }

    #[test]
    fn is_identity_check() {
        assert!(is_identity(&Transform3::identity()));
        assert!(!is_identity(&Transform3::translation(1.0, 0.0, 0.0)));
    }

    #[test]
    fn export_assembly_with_nonidentity_transform() {
        let doc = StepDocument {
            products: vec![
                StepProductNode {
                    name: "Asm".into(),
                    id: "A".into(),
                    children: vec![AssemblyChild {
                        product_index: 1,
                        transform: Transform3::translation(10.0, 0.0, 0.0),
                        instance_name: None,
                    }],
                    shapes: vec![],
                    metadata: ProductMetadata::default(),
                },
                StepProductNode {
                    name: "Child".into(),
                    id: "C".into(),
                    children: vec![],
                    shapes: vec![],
                    metadata: ProductMetadata::default(),
                },
            ],
            roots: vec![0],
            units: UnitContext::default(),
            opaque_bag: OpaqueBag::default(),
        };

        let mut ctx = ExportContext::new();
        ctx.export_product_structure(&doc, &[vec![], vec![]])
            .unwrap();

        // Should have an AXIS2_PLACEMENT_3D for the transform
        let count =
            |name: &str| -> usize { ctx.entities.iter().filter(|e| e.type_name == name).count() };
        assert!(count("AXIS2_PLACEMENT_3D") >= 1);
    }
}
