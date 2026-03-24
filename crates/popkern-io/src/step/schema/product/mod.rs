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
#[path = "tests.rs"]
mod tests;
