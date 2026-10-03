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
    /// The measure type when typed (`LENGTH_MEASURE`, `PLANE_ANGLE_MEASURE`).
    pub measure: Option<String>,
    /// Reference to the unit entity.
    pub unit: EntityId,
    /// Descriptive name.
    pub name: String,
}

impl UncertaintyMeasureWithUnit {
    const TYPE: &str = "UNCERTAINTY_MEASURE_WITH_UNIT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        // The value is a typed measure (`LENGTH_MEASURE(1.E-06)`), the
        // form every writer uses; a bare real is accepted too.
        let value = param_at(params, 0, Self::TYPE)?;
        let measure = match value {
            Parameter::Typed(name, _) => Some(name.clone()),
            _ => None,
        };
        Ok(Self {
            value: as_measure_value(value).map_or_else(|| as_real(value, Self::TYPE), Ok)?,
            measure,
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

/// A measure's number: plain, or typed (`LENGTH_MEASURE(25.4)`, which
/// the parser holds as a one-element list).
pub fn as_measure_value(param: &Parameter) -> Option<f64> {
    match param {
        Parameter::Real(v) => Some(*v),
        Parameter::Int(v) => Some(*v as f64),
        Parameter::Typed(_, inner) => as_measure_value(inner),
        Parameter::List(items) if items.len() == 1 => as_measure_value(&items[0]),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step::parser::types::Parameter;

    fn p_str(s: &str) -> Parameter {
        Parameter::String(s.into())
    }
    fn p_real(v: f64) -> Parameter {
        Parameter::Real(v)
    }
    fn p_ref(id: u64) -> Parameter {
        Parameter::EntityRef(id)
    }
    fn p_int(v: i64) -> Parameter {
        Parameter::Int(v)
    }
    fn p_enum(s: &str) -> Parameter {
        Parameter::Enum(s.into())
    }

    #[test]
    fn geometric_representation_context_parse() {
        let params = vec![p_str("3d_context"), p_int(3)];
        let grc = GeometricRepresentationContext::from_params(&params).unwrap();
        assert_eq!(grc.id, "3d_context");
        assert_eq!(grc.dimension, 3);
    }

    #[test]
    fn global_uncertainty_assigned_context_parse() {
        let params = vec![p_str("ctx"), p_int(3), Parameter::List(vec![p_ref(99)])];
        let guac = GlobalUncertaintyAssignedContext::from_params(&params).unwrap();
        assert_eq!(guac.id, "ctx");
        assert_eq!(guac.dimension, 3);
        assert_eq!(guac.uncertainty, vec![99]);
    }

    #[test]
    fn uncertainty_measure_with_unit_parse() {
        let params = vec![p_real(0.001), p_ref(50), p_str("tolerance")];
        let um = UncertaintyMeasureWithUnit::from_params(&params).unwrap();
        assert!((um.value - 0.001).abs() < 1e-12);
        assert_eq!(um.unit, 50);
        assert_eq!(um.name, "tolerance");
    }

    #[test]
    fn si_unit_with_prefix() {
        let params = vec![Parameter::Derived, p_enum("MILLI"), p_enum("METRE")];
        let su = SiUnit::from_params(&params).unwrap();
        assert!(su.dimensions.is_none());
        assert_eq!(su.prefix, Some("MILLI".to_string()));
        assert_eq!(su.name, "METRE");
    }

    #[test]
    fn si_unit_with_dimension_ref() {
        let params = vec![p_ref(5), Parameter::Omitted, p_enum("RADIAN")];
        let su = SiUnit::from_params(&params).unwrap();
        assert_eq!(su.dimensions, Some(5));
        assert!(su.prefix.is_none());
        assert_eq!(su.name, "RADIAN");
    }

    #[test]
    fn conversion_based_unit_parse() {
        let params = vec![p_ref(10), p_str("INCH"), p_ref(20)];
        let cbu = ConversionBasedUnit::from_params(&params).unwrap();
        assert_eq!(cbu.dimensions, Some(10));
        assert_eq!(cbu.name, "INCH");
        assert_eq!(cbu.conversion_factor, 20);
    }

    #[test]
    fn conversion_based_unit_no_dimensions() {
        let params = vec![Parameter::Omitted, p_str("DEGREE"), p_ref(30)];
        let cbu = ConversionBasedUnit::from_params(&params).unwrap();
        assert!(cbu.dimensions.is_none());
        assert_eq!(cbu.name, "DEGREE");
        assert_eq!(cbu.conversion_factor, 30);
    }

    #[test]
    fn length_measure_with_unit_parse() {
        let params = vec![p_real(25.4), p_ref(7)];
        let lm = LengthMeasureWithUnit::from_params(&params).unwrap();
        assert!((lm.value - 25.4).abs() < 1e-12);
        assert_eq!(lm.unit, 7);
    }

    #[test]
    fn plane_angle_measure_with_unit_parse() {
        let params = vec![p_real(0.01745), p_ref(8)];
        let pa = PlaneAngleMeasureWithUnit::from_params(&params).unwrap();
        assert!((pa.value - 0.01745).abs() < 1e-12);
        assert_eq!(pa.unit, 8);
    }
}
