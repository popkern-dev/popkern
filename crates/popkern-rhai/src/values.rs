//! Serde-backed conversions of reports, evidence, and parameter
//! listings into script values (maps and arrays). Floats cross
//! bit-exactly: `to_dynamic` places an `f64` directly into the value.

use popkern_script::ParamInfo;
use rhai::{Dynamic, Map};
use serde::Serialize;

use crate::args::unit_name;

/// Serialise any report-style record into a script map/array.
pub fn to_value<T: Serialize>(value: &T) -> Dynamic {
    rhai::serde::to_dynamic(value).unwrap_or(Dynamic::UNIT)
}

/// A parameter row as `#{ name, value, unit, expression, bounds, read_only }`.
pub fn param_map(info: &ParamInfo) -> Dynamic {
    let mut map = Map::new();
    map.insert("name".into(), info.name.clone().into());
    map.insert("value".into(), Dynamic::from(info.value));
    map.insert("unit".into(), unit_name(info.unit).into());
    map.insert(
        "expression".into(),
        match &info.expression {
            Some(text) => text.clone().into(),
            None => Dynamic::UNIT,
        },
    );
    map.insert(
        "bounds".into(),
        match info.bounds {
            Some((lo, hi)) => Dynamic::from_array(vec![Dynamic::from(lo), Dynamic::from(hi)]),
            None => Dynamic::UNIT,
        },
    );
    map.insert("read_only".into(), info.read_only.into());
    Dynamic::from_map(map)
}

/// `#{ label, name }` for an entity label row.
pub fn label_map(label: &str, name: &popkern_core::naming::EntityName) -> Dynamic {
    let mut map = Map::new();
    map.insert("label".into(), label.to_string().into());
    map.insert("name".into(), format!("{name:?}").into());
    Dynamic::from_map(map)
}
