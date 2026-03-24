//! Core types for the Part 21 parser.

/// A STEP entity instance ID (the `#N` in `#1 = CARTESIAN_POINT(...)`).
pub type EntityId = u64;

/// A single parameter value in a STEP entity instance.
#[derive(Debug, Clone, PartialEq)]
pub enum Parameter {
    /// Integer value.
    Int(i64),
    /// Real (floating-point) value.
    Real(f64),
    /// String value (decoded from Part 21 escaping).
    String(String),
    /// Enumeration value (e.g., `.BOTH.`).
    Enum(String),
    /// Reference to another entity (`#N`).
    EntityRef(EntityId),
    /// Aggregate (list) of parameters.
    List(Vec<Parameter>),
    /// Typed parameter for complex entity instances.
    Typed(String, Box<Parameter>),
    /// Omitted parameter (`$`).
    Omitted,
    /// Derived parameter (`*`).
    Derived,
}

/// A parsed STEP entity instance.
#[derive(Debug, Clone)]
pub struct EntityInstance {
    /// The entity ID (`#N`).
    pub id: EntityId,
    /// The entity type name (e.g., `CARTESIAN_POINT`).
    pub type_name: String,
    /// Positional parameter list.
    pub params: Vec<Parameter>,
    /// Additional type/param pairs for complex (AND-connected) entities.
    pub complex_parts: Vec<(String, Vec<Parameter>)>,
}

/// Parsed content of the HEADER section.
#[derive(Debug, Clone, Default)]
pub struct StepHeader {
    /// File description strings.
    pub description: Vec<String>,
    /// Implementation level (e.g., "2;1").
    pub implementation_level: String,
    /// File name.
    pub name: String,
    /// File schema identifiers (e.g., \["AUTOMOTIVE_DESIGN"\]).
    pub schemas: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parameter_int() {
        let p = Parameter::Int(42);
        assert_eq!(p, Parameter::Int(42));
    }

    #[test]
    fn parameter_entity_ref() {
        let p = Parameter::EntityRef(7);
        assert_eq!(p, Parameter::EntityRef(7));
    }

    #[test]
    fn parameter_list_nested() {
        let p = Parameter::List(vec![Parameter::Real(1.0), Parameter::Real(2.0)]);
        if let Parameter::List(items) = &p {
            assert_eq!(items.len(), 2);
        } else {
            panic!("expected list");
        }
    }

    #[test]
    fn entity_instance_creation() {
        let e = EntityInstance {
            id: 1,
            type_name: "CARTESIAN_POINT".into(),
            params: vec![
                Parameter::String("".into()),
                Parameter::List(vec![
                    Parameter::Real(0.0),
                    Parameter::Real(0.0),
                    Parameter::Real(0.0),
                ]),
            ],
            complex_parts: vec![],
        };
        assert_eq!(e.id, 1);
        assert_eq!(e.type_name, "CARTESIAN_POINT");
    }
}
