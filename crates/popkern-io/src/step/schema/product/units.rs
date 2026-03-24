//! Context and unit entity definitions from the STEP schema.

use crate::step::parser::types::{EntityId, Parameter};
use crate::step::schema::{
    MappingError, as_entity_ref, as_int, as_optional_ref, as_real, as_ref_list, as_string, param_at,
};

/// GEOMETRIC_REPRESENTATION_CONTEXT — context providing coordinate space dimension.
#[derive(Debug, Clone)]
pub struct GeometricRepresentationContext {
    /// Context identifier.
    pub id: String,
    /// Coordinate space dimension (typically 3).
    pub dimension: i64,
}

impl GeometricRepresentationContext {
    const TYPE: &str = "GEOMETRIC_REPRESENTATION_CONTEXT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            id: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            dimension: as_int(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT — context with assigned uncertainty values.
#[derive(Debug, Clone)]
pub struct GlobalUncertaintyAssignedContext {
    /// Context identifier.
    pub id: String,
    /// Coordinate space dimension.
    pub dimension: i64,
    /// References to `UNCERTAINTY_MEASURE_WITH_UNIT` entities.
    pub uncertainty: Vec<EntityId>,
}

impl GlobalUncertaintyAssignedContext {
    const TYPE: &str = "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            id: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            dimension: as_int(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            uncertainty: as_ref_list(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// UNCERTAINTY_MEASURE_WITH_UNIT — an uncertainty value with associated unit.
#[derive(Debug, Clone)]
pub struct UncertaintyMeasureWithUnit {
    /// The uncertainty value.
    pub value: f64,
    /// Reference to the unit entity.
    pub unit: EntityId,
    /// Descriptive name.
    pub name: String,
}

impl UncertaintyMeasureWithUnit {
    const TYPE: &str = "UNCERTAINTY_MEASURE_WITH_UNIT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            value: as_real(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            unit: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            name: as_string(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// SI_UNIT — an SI unit with optional prefix.
#[derive(Debug, Clone)]
pub struct SiUnit {
    /// Optional reference to dimensional exponents.
    pub dimensions: Option<EntityId>,
    /// Optional SI prefix (e.g. `MILLI`, `CENTI`).
    pub prefix: Option<String>,
    /// Unit name (e.g. `METRE`, `RADIAN`).
    pub name: String,
}

impl SiUnit {
    const TYPE: &str = "SI_UNIT";

    /// Parse from raw STEP parameters.
    ///
    /// SI_UNIT params: dimensions (often `*`), prefix (enum or `$`), name (enum).
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        let dimensions = as_optional_ref(param_at(params, 0, Self::TYPE)?, Self::TYPE)?;
        let prefix = match param_at(params, 1, Self::TYPE)? {
            Parameter::Enum(s) => Some(s.clone()),
            Parameter::Omitted | Parameter::Derived => None,
            other => {
                return Err(MappingError {
                    entity_type: Self::TYPE.to_string(),
                    entity_id: None,
                    message: format!("expected Enum or Omitted for prefix, got {other:?}"),
                });
            }
        };
        let name = match param_at(params, 2, Self::TYPE)? {
            Parameter::Enum(s) => s.clone(),
            other => {
                return Err(MappingError {
                    entity_type: Self::TYPE.to_string(),
                    entity_id: None,
                    message: format!("expected Enum for unit name, got {other:?}"),
                });
            }
        };
        Ok(Self {
            dimensions,
            prefix,
            name,
        })
    }
}

/// CONVERSION_BASED_UNIT — a unit defined by conversion from a base unit.
#[derive(Debug, Clone)]
pub struct ConversionBasedUnit {
    /// Optional reference to dimensional exponents.
    pub dimensions: Option<EntityId>,
    /// Unit name.
    pub name: String,
    /// Reference to the conversion factor entity.
    pub conversion_factor: EntityId,
}

impl ConversionBasedUnit {
    const TYPE: &str = "CONVERSION_BASED_UNIT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            dimensions: as_optional_ref(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            name: as_string(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            conversion_factor: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// LENGTH_MEASURE_WITH_UNIT — a length value with associated unit.
#[derive(Debug, Clone)]
pub struct LengthMeasureWithUnit {
    /// The length value.
    pub value: f64,
    /// Reference to the unit entity.
    pub unit: EntityId,
}

impl LengthMeasureWithUnit {
    const TYPE: &str = "LENGTH_MEASURE_WITH_UNIT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            value: as_real(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            unit: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// PLANE_ANGLE_MEASURE_WITH_UNIT — an angular value with associated unit.
#[derive(Debug, Clone)]
pub struct PlaneAngleMeasureWithUnit {
    /// The angle value.
    pub value: f64,
    /// Reference to the unit entity.
    pub unit: EntityId,
}

impl PlaneAngleMeasureWithUnit {
    const TYPE: &str = "PLANE_ANGLE_MEASURE_WITH_UNIT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            value: as_real(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            unit: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}
