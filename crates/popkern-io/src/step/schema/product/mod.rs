//! Typed product-structure entity definitions from the STEP schema.
//!
//! Context and unit entities are in the [`units`] submodule.

mod units;

pub use units::*;

use super::{MappingError, as_entity_ref, as_ref_list, as_string, param_at};
use crate::step::parser::types::{EntityId, Parameter};

/// PRODUCT — the top-level product identification entity.
#[derive(Debug, Clone)]
pub struct StepProduct {
    /// Product identifier.
    pub id: String,
    /// Product name.
    pub name: String,
    /// Product description.
    pub description: String,
    /// References to `PRODUCT_CONTEXT` entities.
    pub frame_of_reference: Vec<EntityId>,
}

impl StepProduct {
    const TYPE: &str = "PRODUCT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            id: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            name: as_string(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            description: as_string(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            frame_of_reference: as_ref_list(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// PRODUCT_DEFINITION_FORMATION — a version/revision of a product.
#[derive(Debug, Clone)]
pub struct ProductDefinitionFormation {
    /// Formation identifier.
    pub id: String,
    /// Description.
    pub description: String,
    /// Reference to the `PRODUCT`.
    pub of_product: EntityId,
}

impl ProductDefinitionFormation {
    const TYPE: &str = "PRODUCT_DEFINITION_FORMATION";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            id: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            description: as_string(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            of_product: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// PRODUCT_DEFINITION — a definition within a specific context.
#[derive(Debug, Clone)]
pub struct ProductDefinition {
    /// Definition identifier.
    pub id: String,
    /// Description.
    pub description: String,
    /// Reference to the `PRODUCT_DEFINITION_FORMATION`.
    pub formation: EntityId,
    /// Reference to the `PRODUCT_DEFINITION_CONTEXT`.
    pub frame_of_reference: EntityId,
}

impl ProductDefinition {
    const TYPE: &str = "PRODUCT_DEFINITION";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            id: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            description: as_string(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            formation: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            frame_of_reference: as_entity_ref(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// PRODUCT_DEFINITION_SHAPE — shape associated with a product definition.
#[derive(Debug, Clone)]
pub struct ProductDefinitionShape {
    /// Shape name.
    pub name: String,
    /// Description.
    pub description: String,
    /// Reference to the `PRODUCT_DEFINITION`.
    pub definition: EntityId,
}

impl ProductDefinitionShape {
    const TYPE: &str = "PRODUCT_DEFINITION_SHAPE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            description: as_string(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            definition: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// SHAPE_DEFINITION_REPRESENTATION — link between shape definition and representation.
#[derive(Debug, Clone)]
pub struct ShapeDefinitionRepresentation {
    /// Reference to the `PRODUCT_DEFINITION_SHAPE`.
    pub definition: EntityId,
    /// Reference to the `SHAPE_REPRESENTATION`.
    pub used_representation: EntityId,
}

impl ShapeDefinitionRepresentation {
    const TYPE: &str = "SHAPE_DEFINITION_REPRESENTATION";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            definition: as_entity_ref(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            used_representation: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// SHAPE_REPRESENTATION — a representation containing geometric items.
#[derive(Debug, Clone)]
pub struct StepShapeRepresentation {
    /// Representation name.
    pub name: String,
    /// References to representation items.
    pub items: Vec<EntityId>,
    /// Reference to the `REPRESENTATION_CONTEXT`.
    pub context_of_items: EntityId,
}

impl StepShapeRepresentation {
    const TYPE: &str = "SHAPE_REPRESENTATION";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            items: as_ref_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            context_of_items: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// ADVANCED_BREP_SHAPE_REPRESENTATION — shape representation for advanced BRep.
#[derive(Debug, Clone)]
pub struct AdvancedBrepShapeRepresentation {
    /// Representation name.
    pub name: String,
    /// References to representation items.
    pub items: Vec<EntityId>,
    /// Reference to the `REPRESENTATION_CONTEXT`.
    pub context_of_items: EntityId,
}

impl AdvancedBrepShapeRepresentation {
    const TYPE: &str = "ADVANCED_BREP_SHAPE_REPRESENTATION";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            items: as_ref_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            context_of_items: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// FACETED_BREP_SHAPE_REPRESENTATION — shape representation for faceted BRep.
#[derive(Debug, Clone)]
pub struct FacetedBrepShapeRepresentation {
    /// Representation name.
    pub name: String,
    /// References to representation items.
    pub items: Vec<EntityId>,
    /// Reference to the `REPRESENTATION_CONTEXT`.
    pub context_of_items: EntityId,
}

impl FacetedBrepShapeRepresentation {
    const TYPE: &str = "FACETED_BREP_SHAPE_REPRESENTATION";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            items: as_ref_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            context_of_items: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// NEXT_ASSEMBLY_USAGE_OCCURRENCE — parent-child assembly relationship.
#[derive(Debug, Clone)]
pub struct NextAssemblyUsageOccurrence {
    /// Relationship identifier.
    pub id: String,
    /// Relationship name.
    pub name: String,
    /// Description.
    pub description: String,
    /// Reference to the parent `PRODUCT_DEFINITION`.
    pub relating: EntityId,
    /// Reference to the child `PRODUCT_DEFINITION`.
    pub related: EntityId,
    /// Optional reference designator string.
    pub reference_designator: Option<String>,
}

impl NextAssemblyUsageOccurrence {
    const TYPE: &str = "NEXT_ASSEMBLY_USAGE_OCCURRENCE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        let ref_des = match param_at(params, 5, Self::TYPE) {
            Ok(Parameter::String(s)) => Some(s.clone()),
            Ok(Parameter::Omitted) | Err(_) => None,
            Ok(other) => {
                return Err(MappingError {
                    entity_type: Self::TYPE.to_string(),
                    entity_id: None,
                    message: format!(
                        "expected String or Omitted for ref designator, got {other:?}"
                    ),
                });
            }
        };
        Ok(Self {
            id: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            name: as_string(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            description: as_string(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            relating: as_entity_ref(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
            related: as_entity_ref(param_at(params, 4, Self::TYPE)?, Self::TYPE)?,
            reference_designator: ref_des,
        })
    }
}

/// REPRESENTATION_RELATIONSHIP — link between two representations.
#[derive(Debug, Clone)]
pub struct RepresentationRelationship {
    /// Relationship name.
    pub name: String,
    /// Description.
    pub description: String,
    /// Reference to the first representation.
    pub rep_1: EntityId,
    /// Reference to the second representation.
    pub rep_2: EntityId,
}

impl RepresentationRelationship {
    const TYPE: &str = "REPRESENTATION_RELATIONSHIP";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            description: as_string(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            rep_1: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            rep_2: as_entity_ref(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION — relationship carrying a transform.
#[derive(Debug, Clone)]
pub struct RepresentationRelationshipWithTransformation {
    /// Relationship name.
    pub name: String,
    /// Description.
    pub description: String,
    /// Reference to the first representation.
    pub rep_1: EntityId,
    /// Reference to the second representation.
    pub rep_2: EntityId,
    /// Reference to the transformation entity.
    pub transformation: EntityId,
}

impl RepresentationRelationshipWithTransformation {
    const TYPE: &str = "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            description: as_string(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            rep_1: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            rep_2: as_entity_ref(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
            transformation: as_entity_ref(param_at(params, 4, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// ITEM_DEFINED_TRANSFORMATION — transformation between two placement items.
#[derive(Debug, Clone)]
pub struct ItemDefinedTransformation {
    /// Transformation name.
    pub name: String,
    /// Description.
    pub description: String,
    /// Reference to the source placement.
    pub transform_item_1: EntityId,
    /// Reference to the target placement.
    pub transform_item_2: EntityId,
}

impl ItemDefinedTransformation {
    const TYPE: &str = "ITEM_DEFINED_TRANSFORMATION";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            description: as_string(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            transform_item_1: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            transform_item_2: as_entity_ref(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p_str(s: &str) -> Parameter {
        Parameter::String(s.into())
    }
    fn p_ref(id: u64) -> Parameter {
        Parameter::EntityRef(id)
    }
    fn p_enum(s: &str) -> Parameter {
        Parameter::Enum(s.into())
    }

    #[test]
    fn step_product_parse() {
        let params = vec![
            p_str("P001"),
            p_str("Widget"),
            p_str("A test widget"),
            Parameter::List(vec![p_ref(10)]),
        ];
        let p = StepProduct::from_params(&params).unwrap();
        assert_eq!(p.id, "P001");
        assert_eq!(p.name, "Widget");
        assert_eq!(p.frame_of_reference, vec![10]);
    }

    #[test]
    fn product_definition_parse() {
        let params = vec![p_str("def1"), p_str("desc"), p_ref(5), p_ref(6)];
        let pd = ProductDefinition::from_params(&params).unwrap();
        assert_eq!(pd.id, "def1");
        assert_eq!(pd.formation, 5);
        assert_eq!(pd.frame_of_reference, 6);
    }

    #[test]
    fn shape_representation_parse() {
        let params = vec![
            p_str("shape"),
            Parameter::List(vec![p_ref(1), p_ref(2)]),
            p_ref(3),
        ];
        let sr = StepShapeRepresentation::from_params(&params).unwrap();
        assert_eq!(sr.name, "shape");
        assert_eq!(sr.items, vec![1, 2]);
        assert_eq!(sr.context_of_items, 3);
    }

    #[test]
    fn next_assembly_usage_occurrence_parse() {
        let params = vec![
            p_str("rel1"),
            p_str("asm"),
            p_str("desc"),
            p_ref(10),
            p_ref(20),
            p_str("R1"),
        ];
        let nauo = NextAssemblyUsageOccurrence::from_params(&params).unwrap();
        assert_eq!(nauo.relating, 10);
        assert_eq!(nauo.related, 20);
        assert_eq!(nauo.reference_designator, Some("R1".to_string()));
    }

    #[test]
    fn product_definition_formation_parse() {
        let params = vec![p_str("v1"), p_str("first version"), p_ref(100)];
        let pdf = ProductDefinitionFormation::from_params(&params).unwrap();
        assert_eq!(pdf.id, "v1");
        assert_eq!(pdf.description, "first version");
        assert_eq!(pdf.of_product, 100);
    }

    #[test]
    fn product_definition_shape_parse() {
        let params = vec![p_str("shape1"), p_str("main shape"), p_ref(50)];
        let pds = ProductDefinitionShape::from_params(&params).unwrap();
        assert_eq!(pds.name, "shape1");
        assert_eq!(pds.description, "main shape");
        assert_eq!(pds.definition, 50);
    }

    #[test]
    fn shape_definition_representation_parse() {
        let params = vec![p_ref(10), p_ref(20)];
        let sdr = ShapeDefinitionRepresentation::from_params(&params).unwrap();
        assert_eq!(sdr.definition, 10);
        assert_eq!(sdr.used_representation, 20);
    }

    #[test]
    fn advanced_brep_shape_representation_parse() {
        let params = vec![
            p_str("abrep"),
            Parameter::List(vec![p_ref(1), p_ref(2)]),
            p_ref(3),
        ];
        let ab = AdvancedBrepShapeRepresentation::from_params(&params).unwrap();
        assert_eq!(ab.name, "abrep");
        assert_eq!(ab.items, vec![1, 2]);
        assert_eq!(ab.context_of_items, 3);
    }

    #[test]
    fn faceted_brep_shape_representation_parse() {
        let params = vec![p_str("fbrep"), Parameter::List(vec![p_ref(4)]), p_ref(5)];
        let fb = FacetedBrepShapeRepresentation::from_params(&params).unwrap();
        assert_eq!(fb.name, "fbrep");
        assert_eq!(fb.items, vec![4]);
        assert_eq!(fb.context_of_items, 5);
    }

    #[test]
    fn nauo_omitted_ref_designator() {
        let params = vec![
            p_str("rel2"),
            p_str("asm"),
            p_str("desc"),
            p_ref(10),
            p_ref(20),
            Parameter::Omitted,
        ];
        let nauo = NextAssemblyUsageOccurrence::from_params(&params).unwrap();
        assert_eq!(nauo.reference_designator, None);
    }

    #[test]
    fn representation_relationship_parse() {
        let params = vec![p_str("rr"), p_str("link"), p_ref(1), p_ref(2)];
        let rr = RepresentationRelationship::from_params(&params).unwrap();
        assert_eq!(rr.name, "rr");
        assert_eq!(rr.description, "link");
        assert_eq!(rr.rep_1, 1);
        assert_eq!(rr.rep_2, 2);
    }

    #[test]
    fn representation_relationship_with_transformation_parse() {
        let params = vec![p_str("rrt"), p_str("xform"), p_ref(1), p_ref(2), p_ref(3)];
        let rrt = RepresentationRelationshipWithTransformation::from_params(&params).unwrap();
        assert_eq!(rrt.name, "rrt");
        assert_eq!(rrt.rep_1, 1);
        assert_eq!(rrt.rep_2, 2);
        assert_eq!(rrt.transformation, 3);
    }

    #[test]
    fn item_defined_transformation_parse() {
        let params = vec![p_str("idt"), p_str("transform"), p_ref(10), p_ref(20)];
        let idt = ItemDefinedTransformation::from_params(&params).unwrap();
        assert_eq!(idt.name, "idt");
        assert_eq!(idt.description, "transform");
        assert_eq!(idt.transform_item_1, 10);
        assert_eq!(idt.transform_item_2, 20);
    }

    #[test]
    fn si_unit_parse() {
        let params = vec![Parameter::Derived, Parameter::Omitted, p_enum("METRE")];
        let su = SiUnit::from_params(&params).unwrap();
        assert!(su.dimensions.is_none());
        assert!(su.prefix.is_none());
        assert_eq!(su.name, "METRE");
    }
}
