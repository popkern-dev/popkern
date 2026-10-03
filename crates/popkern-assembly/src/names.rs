//! Instance paths and entity paths.
//!
//! An [`InstancePath`] is the sibling names from the root joined by
//! `/` (`carrier/pin_2`). An entity of an instance is spelled
//! `<instance path>:<entity name>` (`planet_2:planet/extrude-1/Side/3`):
//! `:` separates because entity spellings already use `/`. Names are
//! therefore non-empty and free of both separators.

use crate::error::AssemblyError;

/// Check that `name` can name an instance or definition.
///
/// # Errors
///
/// [`AssemblyError::InvalidName`] for an empty name, one with `/` or
/// `:`, or one with surrounding whitespace.
pub fn validate_name(name: &str, what: &'static str) -> Result<(), AssemblyError> {
    let reason = if name.is_empty() {
        Some("names cannot be empty")
    } else if name.contains('/') {
        Some("'/' separates instance paths")
    } else if name.contains(':') {
        Some("':' separates an instance path from an entity name")
    } else if name.trim() != name {
        Some("names cannot start or end with whitespace")
    } else {
        None
    };
    match reason {
        Some(reason) => Err(AssemblyError::InvalidName {
            name: name.to_string(),
            what,
            reason,
        }),
        None => Ok(()),
    }
}

/// The names of the instances from the root down to one occurrence.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct InstancePath(pub Vec<String>);

impl InstancePath {
    /// Parse `a/b/c`.
    pub fn parse(spelling: &str) -> InstancePath {
        InstancePath(spelling.split('/').map(str::to_string).collect())
    }

    /// One level deeper.
    pub fn child(&self, name: &str) -> InstancePath {
        let mut names = self.0.clone();
        names.push(name.to_string());
        InstancePath(names)
    }

    /// Depth below the root (a top-level instance is 0).
    pub fn depth(&self) -> usize {
        self.0.len().saturating_sub(1)
    }
}

impl std::fmt::Display for InstancePath {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0.join("/"))
    }
}

/// Spell an entity of the instance at `path`: `<path>:<entity>`.
pub fn entity_path(path: &InstancePath, entity: &str) -> String {
    format!("{path}:{entity}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_refuse_separators_and_blank() {
        for bad in ["", "a/b", "a:b", " a", "a "] {
            assert!(validate_name(bad, "instance").is_err(), "{bad:?}");
        }
        validate_name("planet_1", "instance").expect("ok");
    }

    #[test]
    fn paths_spell_and_parse() {
        let p = InstancePath::parse("carrier").child("pin_2");
        assert_eq!(p.to_string(), "carrier/pin_2");
        assert_eq!(p.depth(), 1);
        assert_eq!(InstancePath::parse(&p.to_string()), p);
        assert_eq!(
            entity_path(&p, "pin/extrude-1/Side/0"),
            "carrier/pin_2:pin/extrude-1/Side/0"
        );
    }
}
