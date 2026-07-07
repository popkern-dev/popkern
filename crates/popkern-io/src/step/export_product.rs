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
    /// `shape_item_ids` is parallel to `doc.products` and lists each exported
    /// top-level STEP shape item for that product.
    pub(crate) fn export_product_structure(
        &mut self,
        doc: &StepDocument,
        shape_item_ids: &[Vec<EntityId>],
    ) -> Result<(), IoError> {
        // Emit shared context entities
        let ctx_id = self.export_representation_context(&doc.units.context);

        // APPLICATION_CONTEXT shared by the product contexts (their
        // frame_of_reference attribute is a required entity ref).
        let ac_id = self.emit(
            "APPLICATION_CONTEXT",
            vec![str_param(
                "core data for automotive mechanical design processes",
            )],
        );
        // PRODUCT_DEFINITION_CONTEXT (shared by all product definitions)
        let pdc_id = self.emit(
            "PRODUCT_DEFINITION_CONTEXT",
            vec![
                str_param("part definition"),
                ref_(ac_id),
                str_param("design"),
            ],
        );

        // PRODUCT_CONTEXT for PRODUCT.frame_of_reference
        let pc_id = self.emit(
            "PRODUCT_CONTEXT",
            vec![str_param(""), ref_(ac_id), str_param("mechanical")],
        );

        // Emit each product
        let mut product_def_ids: Vec<EntityId> = Vec::new();
        for (i, product) in doc.products.iter().enumerate() {
            let shape_items = shape_item_ids.get(i).map(|v| v.as_slice()).unwrap_or(&[]);

            // PRODUCT
            let prod_id = self.emit(
                "PRODUCT",
                vec![
                    str_param(&product.id),
                    str_param(&product.name),
                    str_param(""),
                    list(vec![ref_(pc_id)]),
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
                vec![str_param(""), str_param(""), ref_(pdf_id), ref_(pdc_id)],
            );
            product_def_ids.push(pd_id);

            // PRODUCT_DEFINITION_SHAPE
            let pds_id = self.emit(
                "PRODUCT_DEFINITION_SHAPE",
                vec![str_param(""), str_param(""), ref_(pd_id)],
            );

            // Build SHAPE_REPRESENTATION containing the solid items
            if !shape_items.is_empty() {
                let items: Vec<_> = shape_items.iter().map(|&id| ref_(id)).collect();

                // Determine the representation type name
                let rep_type = if product
                    .shapes
                    .iter()
                    .any(|s| s.kind == ShapeRepKind::AdvancedBrep)
                    && !items_include_surface_model(self, shape_items)
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

    /// Emit an SI unit as the standard complex instance
    /// `(<KIND>() NAMED_UNIT(*) SI_UNIT(<prefix>,.<NAME>.))` —
    /// standalone SI_UNIT lacks the inherited NAMED_UNIT.dimensions
    /// and the unit-kind supertype, which strict readers flag
    /// ("Count of Parameters is not 3 for si_unit").
    fn emit_si_unit(&mut self, kind: &str, prefix: Option<&str>, name: &str) -> EntityId {
        use crate::step::parser::types::Parameter;
        let id = self.alloc_id();
        // Components in alphabetical order per Part 21.
        let mut parts: Vec<(String, Vec<Parameter>)> = Vec::new();
        let si_params = vec![prefix.map(enum_).unwrap_or(Parameter::Omitted), enum_(name)];
        let mut all: Vec<(&str, Vec<Parameter>)> = vec![
            (kind, vec![]),
            ("NAMED_UNIT", vec![Parameter::Derived]),
            ("SI_UNIT", si_params),
        ];
        all.sort_by(|a, b| a.0.cmp(b.0));
        let Some((first, rest)) = all.split_first() else {
            // Unreachable: 'all' is a 3-element literal; the branch
            // exists only to keep library code unwrap-free.
            return id;
        };
        let type_name = first.0.to_string();
        let params = first.1.clone();
        for (t, ps) in rest {
            parts.push((t.to_string(), ps.clone()));
        }
        self.entities
            .push(crate::step::parser::types::EntityInstance {
                id,
                type_name,
                params,
                complex_parts: parts,
            });
        id
    }

    /// Export a geometric representation context with units.
    fn export_representation_context(&mut self, ctx: &RepresentationContext) -> EntityId {
        // SI length unit
        let lu_id = match ctx.length_unit {
            LengthUnit::Millimetre => self.emit_si_unit("LENGTH_UNIT", Some("MILLI"), "METRE"),
            LengthUnit::Metre => self.emit_si_unit("LENGTH_UNIT", None, "METRE"),
            LengthUnit::Inch => {
                // Conversion-based unit
                let si = self.emit_si_unit("LENGTH_UNIT", None, "METRE");
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
            AngleUnit::Radian => self.emit_si_unit("PLANE_ANGLE_UNIT", None, "RADIAN"),
            AngleUnit::Degree => {
                let si = self.emit_si_unit("PLANE_ANGLE_UNIT", None, "RADIAN");
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
                crate::step::parser::types::Parameter::Typed(
                    "LENGTH_MEASURE".into(),
                    Box::new(real(ctx.uncertainty)),
                ),
                ref_(lu_id),
                str_param("distance_accuracy_value"),
                str_param(
                    "maximum model space distance between geometric entities considered equal",
                ),
            ],
        );

        // STEP AP203/214/242 mandates that the representation context
        // for shape data be a *complex* entity that multiply-inherits
        // from GEOMETRIC_REPRESENTATION_CONTEXT, GLOBAL_UNIT_ASSIGNED_-
        // CONTEXT, GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT, and
        // REPRESENTATION_CONTEXT. The simple-entity form
        //   GEOMETRIC_REPRESENTATION_CONTEXT('','3D',3,#u,#a,#unc)
        // is non-standard — OCCT-strict importers (Bambu Studio,
        // FreeCAD) ignore the unit assignment and prompt the user to
        // pick units, then often pick the wrong scale. Emit the
        // complex form so units (millimetres) round-trip cleanly.
        let id = self.alloc_id();
        self.entities
            .push(crate::step::parser::types::EntityInstance {
                id,
                type_name: "GEOMETRIC_REPRESENTATION_CONTEXT".into(),
                params: vec![int(ctx.dimension as i64)],
                complex_parts: vec![
                    (
                        "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT".into(),
                        vec![crate::step::export::list(vec![ref_(unc_id)])],
                    ),
                    (
                        "GLOBAL_UNIT_ASSIGNED_CONTEXT".into(),
                        vec![crate::step::export::list(vec![ref_(lu_id), ref_(au_id)])],
                    ),
                    (
                        "REPRESENTATION_CONTEXT".into(),
                        vec![str_param(""), str_param("3D")],
                    ),
                ],
            });
        id
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

fn items_include_surface_model(ctx: &ExportContext, item_ids: &[EntityId]) -> bool {
    item_ids.iter().any(|item_id| {
        ctx.entities
            .iter()
            .any(|entity| entity.id == *item_id && entity.type_name == "SHELL_BASED_SURFACE_MODEL")
    })
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
