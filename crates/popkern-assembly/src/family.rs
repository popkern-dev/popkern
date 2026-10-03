//! Part families (design §6.1, PK-630): one parametric definition,
//! one built variant per distinct parameter set.
//!
//! A variant records the family it came from and its full, canonical
//! parameter set: sorted keys, typed values. Two requests for the
//! same set are the same variant (the closure runs once). The
//! variant's display name spells only the values that differ from the
//! family defaults: `bolt[size=M4, length=10]`.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// One parameter value of a family.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum ParamValue {
    /// `true` / `false`.
    Bool(bool),
    /// A number (integers are numbers).
    Number(f64),
    /// A string (`"M6"`, `"socket"`).
    Text(String),
}

impl ParamValue {
    /// `"a number"`, `"a string"`, `"a bool"`.
    pub fn type_name(&self) -> &'static str {
        match self {
            ParamValue::Bool(_) => "a bool",
            ParamValue::Number(_) => "a number",
            ParamValue::Text(_) => "a string",
        }
    }

    /// Same value for caching: numbers compare by their bits (so
    /// `-0.0` and `0.0` are distinct sets, as they may build distinct
    /// geometry), strings and bools by value.
    pub fn same(&self, other: &ParamValue) -> bool {
        match (self, other) {
            (ParamValue::Number(a), ParamValue::Number(b)) => a.to_bits() == b.to_bits(),
            _ => self == other,
        }
    }
}

impl std::fmt::Display for ParamValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParamValue::Bool(b) => write!(f, "{b}"),
            ParamValue::Number(n) => write!(f, "{n}"),
            ParamValue::Text(s) => f.write_str(s),
        }
    }
}

/// Which family a definition belongs to and with what parameters.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FamilyVariant {
    /// The family name (the default variant's definition name).
    pub family: String,
    /// The full parameter set this variant was built with.
    pub params: BTreeMap<String, ParamValue>,
}

impl FamilyVariant {
    /// Whether `params` is exactly this variant's set.
    pub fn matches(&self, family: &str, params: &BTreeMap<String, ParamValue>) -> bool {
        self.family == family
            && self.params.len() == params.len()
            && self
                .params
                .iter()
                .zip(params)
                .all(|((ka, va), (kb, vb))| ka == kb && va.same(vb))
    }

    /// The definition name for this variant: the family name, plus the
    /// values that differ from `defaults` in brackets.
    pub fn display_name(&self, defaults: &BTreeMap<String, ParamValue>) -> String {
        let changed: Vec<String> = self
            .params
            .iter()
            .filter(|(k, v)| defaults.get(*k).is_none_or(|d| !d.same(v)))
            .map(|(k, v)| format!("{k}={v}"))
            .collect();
        if changed.is_empty() {
            self.family.clone()
        } else {
            format!("{}[{}]", self.family, changed.join(", "))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn set(pairs: &[(&str, ParamValue)]) -> BTreeMap<String, ParamValue> {
        pairs
            .iter()
            .map(|(k, v)| (k.to_string(), v.clone()))
            .collect()
    }

    #[test]
    fn variants_name_only_their_overrides_and_match_exactly() {
        let defaults = set(&[
            ("size", ParamValue::Text("M6".into())),
            ("length", ParamValue::Number(14.0)),
        ]);
        let m4 = FamilyVariant {
            family: "bolt".into(),
            params: set(&[
                ("size", ParamValue::Text("M4".into())),
                ("length", ParamValue::Number(10.0)),
            ]),
        };
        assert_eq!(m4.display_name(&defaults), "bolt[length=10, size=M4]");
        let default = FamilyVariant {
            family: "bolt".into(),
            params: defaults.clone(),
        };
        assert_eq!(default.display_name(&defaults), "bolt");
        assert!(default.matches("bolt", &defaults));
        assert!(!m4.matches("bolt", &defaults));
        assert!(!default.matches("nut", &defaults));
    }
}
