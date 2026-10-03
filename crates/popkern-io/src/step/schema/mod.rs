//! Typed STEP schema entity definitions built on top of the raw Part 21 parser.
//!
//! Each submodule defines Rust structs for a family of STEP entity types,
//! with `from_params()` methods that extract typed fields from raw
//! [`Parameter`](super::parser::types::Parameter) slices.

pub mod geometry;
pub mod presentation;
pub mod product;
pub mod registry;
pub mod topology;

use std::fmt;

use super::parser::types::{EntityId, Parameter};

/// Error produced when mapping raw STEP parameters to typed entity structs.
#[derive(Debug, Clone)]
pub struct MappingError {
    /// The STEP entity type being parsed (e.g. `"CARTESIAN_POINT"`).
    pub entity_type: String,
    /// The entity instance id, if known.
    pub entity_id: Option<EntityId>,
    /// Human-readable description of what went wrong.
    pub message: String,
}

impl fmt::Display for MappingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.entity_id {
            Some(id) => write!(f, "{} #{}: {}", self.entity_type, id, self.message),
            None => write!(f, "{}: {}", self.entity_type, self.message),
        }
    }
}

impl std::error::Error for MappingError {}

// ---------------------------------------------------------------------------
// Helper extractors
// ---------------------------------------------------------------------------

/// Retrieve the parameter at `index`, or return a `MappingError`.
pub fn param_at<'a>(
    params: &'a [Parameter],
    index: usize,
    entity_type: &str,
) -> Result<&'a Parameter, MappingError> {
    params.get(index).ok_or_else(|| MappingError {
        entity_type: entity_type.to_string(),
        entity_id: None,
        message: format!("missing parameter at index {index}"),
    })
}

/// Extract an `EntityRef` value.
pub fn as_entity_ref(param: &Parameter, entity_type: &str) -> Result<EntityId, MappingError> {
    match param {
        Parameter::EntityRef(id) => Ok(*id),
        _ => Err(MappingError {
            entity_type: entity_type.to_string(),
            entity_id: None,
            message: format!("expected EntityRef, got {param:?}"),
        }),
    }
}

/// Extract a `String` value.
pub fn as_string(param: &Parameter, entity_type: &str) -> Result<String, MappingError> {
    match param {
        Parameter::String(s) => Ok(s.clone()),
        _ => Err(MappingError {
            entity_type: entity_type.to_string(),
            entity_id: None,
            message: format!("expected String, got {param:?}"),
        }),
    }
}

/// Extract a `Real` value (also accepts `Int` by promotion).
pub fn as_real(param: &Parameter, entity_type: &str) -> Result<f64, MappingError> {
    match param {
        Parameter::Real(v) => Ok(*v),
        Parameter::Int(v) => Ok(*v as f64),
        _ => Err(MappingError {
            entity_type: entity_type.to_string(),
            entity_id: None,
            message: format!("expected Real, got {param:?}"),
        }),
    }
}

/// Extract an `Int` value.
pub fn as_int(param: &Parameter, entity_type: &str) -> Result<i64, MappingError> {
    match param {
        Parameter::Int(v) => Ok(*v),
        _ => Err(MappingError {
            entity_type: entity_type.to_string(),
            entity_id: None,
            message: format!("expected Int, got {param:?}"),
        }),
    }
}

/// Extract a `List` value.
pub fn as_list<'a>(
    param: &'a Parameter,
    entity_type: &str,
) -> Result<&'a [Parameter], MappingError> {
    match param {
        Parameter::List(items) => Ok(items),
        _ => Err(MappingError {
            entity_type: entity_type.to_string(),
            entity_id: None,
            message: format!("expected List, got {param:?}"),
        }),
    }
}

/// Extract an `Enum` value (the string between dots, e.g. `.BOTH.` → `"BOTH"`).
pub fn as_enum(param: &Parameter, entity_type: &str) -> Result<String, MappingError> {
    match param {
        Parameter::Enum(s) => Ok(s.clone()),
        _ => Err(MappingError {
            entity_type: entity_type.to_string(),
            entity_id: None,
            message: format!("expected Enum, got {param:?}"),
        }),
    }
}

/// Extract an optional entity reference — returns `Ok(None)` for `Omitted`/`Derived`.
pub fn as_optional_ref(
    param: &Parameter,
    entity_type: &str,
) -> Result<Option<EntityId>, MappingError> {
    match param {
        Parameter::EntityRef(id) => Ok(Some(*id)),
        Parameter::Omitted | Parameter::Derived => Ok(None),
        _ => Err(MappingError {
            entity_type: entity_type.to_string(),
            entity_id: None,
            message: format!("expected EntityRef or Omitted, got {param:?}"),
        }),
    }
}

/// Extract a STEP LOGICAL (`.T.`, `.F.` or `.U.` unknown) as "asserted
/// true": unknown reads as false. The schema's LOGICAL attributes
/// (a B-spline's `closed_curve`, `self_intersect`, `u_closed`) are
/// writer hints the geometry itself settles; rejecting `.U.` dropped
/// every such curve and surface (measured: 37 entities in
/// jack_in_the_box.stp).
pub fn as_logical(param: &Parameter, entity_type: &str) -> Result<bool, MappingError> {
    match param {
        Parameter::Enum(s) if s == "U" => Ok(false),
        _ => as_bool(param, entity_type),
    }
}

/// Extract a STEP boolean (`.T.` / `.F.` stored as `Enum`).
pub fn as_bool(param: &Parameter, entity_type: &str) -> Result<bool, MappingError> {
    match param {
        Parameter::Enum(s) if s == "T" => Ok(true),
        Parameter::Enum(s) if s == "F" => Ok(false),
        _ => Err(MappingError {
            entity_type: entity_type.to_string(),
            entity_id: None,
            message: format!("expected .T. or .F., got {param:?}"),
        }),
    }
}

/// Extract a list of entity references.
pub fn as_ref_list(param: &Parameter, entity_type: &str) -> Result<Vec<EntityId>, MappingError> {
    let items = as_list(param, entity_type)?;
    items
        .iter()
        .map(|p| as_entity_ref(p, entity_type))
        .collect()
}

/// Extract a list of real values.
pub fn as_real_list(param: &Parameter, entity_type: &str) -> Result<Vec<f64>, MappingError> {
    let items = as_list(param, entity_type)?;
    items.iter().map(|p| as_real(p, entity_type)).collect()
}

/// Extract a list of integer values.
pub fn as_int_list(param: &Parameter, entity_type: &str) -> Result<Vec<i64>, MappingError> {
    let items = as_list(param, entity_type)?;
    items.iter().map(|p| as_int(p, entity_type)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn param_at_success() {
        let params = vec![Parameter::Int(1), Parameter::Real(2.0)];
        assert!(param_at(&params, 0, "TEST").is_ok());
        assert!(param_at(&params, 1, "TEST").is_ok());
    }

    #[test]
    fn param_at_out_of_bounds() {
        let params = vec![Parameter::Int(1)];
        let err = param_at(&params, 5, "TEST").unwrap_err();
        assert!(err.message.contains("index 5"));
    }

    #[test]
    fn as_entity_ref_ok_and_err() {
        assert_eq!(as_entity_ref(&Parameter::EntityRef(42), "T").unwrap(), 42);
        assert!(as_entity_ref(&Parameter::Int(1), "T").is_err());
    }

    #[test]
    fn as_bool_true_false_invalid() {
        assert!(as_bool(&Parameter::Enum("T".into()), "T").unwrap());
        assert!(!as_bool(&Parameter::Enum("F".into()), "T").unwrap());
        assert!(as_bool(&Parameter::Enum("X".into()), "T").is_err());
        assert!(as_bool(&Parameter::Enum("U".into()), "T").is_err());
        assert!(!as_logical(&Parameter::Enum("U".into()), "T").unwrap());
        assert!(as_logical(&Parameter::Enum("T".into()), "T").unwrap());
    }

    #[test]
    fn as_real_promotes_int() {
        assert!((as_real(&Parameter::Int(3), "T").unwrap() - 3.0).abs() < f64::EPSILON);
    }

    #[test]
    fn as_ref_list_works() {
        let list = Parameter::List(vec![Parameter::EntityRef(1), Parameter::EntityRef(2)]);
        let ids = as_ref_list(&list, "T").unwrap();
        assert_eq!(ids, vec![1, 2]);
    }
}
