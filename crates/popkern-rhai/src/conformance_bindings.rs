//! The skins' binding tables for the conformance audit: how each spec
//! op is spelled in Rhai (registered function, positional parameters,
//! option keys, documented defaults), which preview ops the Rhai skin
//! leaves unbound and why, and the `pub fn` each stable op maps to in
//! `popkern-script`.

use std::sync::LazyLock;

use crate::conformance_audit::Binding;

pub(crate) const fn b(
    op: &'static str,
    function: &'static str,
    positional: &'static [&'static str],
    options: &'static [&'static str],
    defaulted: &'static [&'static str],
) -> Binding {
    Binding {
        op,
        function,
        positional,
        options,
        defaulted,
    }
}

/// The Rhai skin's bindings: the model, profile and body verbs here,
/// the fasteners and later families in
/// [`crate::conformance_bindings_more`].
pub static RHAI_BINDINGS: LazyLock<Vec<Binding>> = LazyLock::new(|| {
    [
        CORE_BINDINGS,
        crate::conformance_bindings_more::MORE_BINDINGS,
    ]
    .concat()
});

const CORE_BINDINGS: &[Binding] = &[
    b(
        "model",
        "model",
        &["options"],
        &["unit", "tolerance"],
        &["name"],
    ),
    b("open", "open", &["path"], &[], &[]),
    b("previous", "previous", &[], &[], &[]),
    b("save", "save", &["path"], &[], &[]),
    b("replay", "replay", &[], &[], &[]),
    b("to_brep", "to_step", &[], &[], &[]),
    b("to_part", "to_step", &[], &[], &[]),
    b("to_step", "to_step", &["path"], &[], &[]),
    b("assert_watertight", "assert_watertight", &[], &[], &[]),
    b(
        "assert_step_roundtrip",
        "assert_step_roundtrip",
        &["path"],
        &[],
        &[],
    ),
    b("bodies", "bodies", &[], &[], &[]),
    b("body", "body", &["id_or_label"], &[], &[]),
    b("active_body", "active", &[], &[], &[]),
    b("set_active_body", "set_active", &["body"], &[], &[]),
    b("profile", "profile", &["plane"], &[], &[]),
    b("rect", "rect", &["w", "h"], &[], &[]),
    b("circle", "circle", &["r"], &[], &[]),
    b("ellipse", "ellipse", &["rx", "ry"], &[], &[]),
    b("polygon", "polygon", &["points"], &[], &[]),
    b("spline", "spline", &["points"], &["degree", "weights"], &[]),
    b("path", "path", &["x0", "y0"], &[], &[]),
    b("line_to", "line_to", &["x", "y"], &[], &[]),
    b("arc_to", "arc_to", &["via", "to"], &[], &[]),
    b(
        "conic_to",
        "conic_to",
        &["control", "to", "weight"],
        &[],
        &[],
    ),
    b(
        "spline_to",
        "spline_to",
        &["points"],
        &["degree", "weights", "knots"],
        &[],
    ),
    b("at", "at", &["cx", "cy"], &[], &[]),
    b("sketch", "sketch", &["label"], &[], &[]),
    b("close", "close", &[], &[], &[]),
    b(
        "param",
        "param",
        &["name", "value", "unit"],
        &[],
        &["bounds"],
    ),
    b(
        "param_expr",
        "param_expr",
        &["name", "text", "unit"],
        &[],
        &[],
    ),
    b(
        "set_param",
        "set_param",
        &["name_or_param", "value"],
        &[],
        &[],
    ),
    b("set_expr", "set_expr", &["name", "text"], &[], &[]),
    b("clear_expr", "clear_expr", &["name"], &[], &[]),
    b(
        "set_bounds",
        "set_bounds",
        &["name", "min", "max"],
        &[],
        &[],
    ),
    b("clear_bounds", "clear_bounds", &["name"], &[], &[]),
    b("rename_param", "rename_param", &["old", "new"], &[], &[]),
    b("params", "params", &[], &[], &[]),
    b(
        "add_config",
        "add_config",
        &["name", "description"],
        &[],
        &[],
    ),
    b(
        "add_derived_config",
        "add_derived_config",
        &["name", "parent"],
        &[],
        &[],
    ),
    b(
        "duplicate_config",
        "duplicate_config",
        &["name", "copy"],
        &[],
        &[],
    ),
    b("delete_config", "delete_config", &["name"], &[], &[]),
    b("rename_config", "rename_config", &["old", "new"], &[], &[]),
    b(
        "set_active_config",
        "set_active_config",
        &["name"],
        &[],
        &[],
    ),
    b(
        "set_override",
        "set_override",
        &["config", "param", "value"],
        &[],
        &[],
    ),
    b(
        "clear_override",
        "clear_override",
        &["config", "param"],
        &[],
        &[],
    ),
    b(
        "set_suppression",
        "set_suppression",
        &["config", "feature", "bool"],
        &[],
        &[],
    ),
    b(
        "clear_suppression",
        "clear_suppression",
        &["config", "feature"],
        &[],
        &[],
    ),
    b("configs", "configs", &[], &[], &[]),
    b("config", "config", &["name"], &[], &[]),
    b(
        "extrude",
        "extrude",
        &["end"],
        &["label", "mode"],
        &["face", "direction", "draft_angle", "second"],
    ),
    b(
        "revolve",
        "revolve",
        &["axis", "end"],
        &["label", "mode"],
        &["face"],
    ),
    b(
        "sweep",
        "sweep",
        &["path_body", "path_edge"],
        &["label", "mode", "twist"],
        &["face"],
    ),
    b(
        "loft",
        "loft",
        &["sections"],
        &["label", "mode", "closed"],
        &[],
    ),
    b(
        "pipe",
        "pipe",
        &["path_edge", "radius"],
        &["label", "mode"],
        &[],
    ),
    b(
        "helical_sweep",
        "helical_sweep",
        &["axis", "radius", "pitch", "turns"],
        &["label", "mode", "taper_to", "start_angle", "left_hand"],
        &["face"],
    ),
    b(
        "box",
        "box",
        &["x", "y", "z"],
        &["label", "at", "axis", "ref", "centred"],
        &[],
    ),
    b(
        "cylinder",
        "cylinder",
        &["r", "h"],
        &["label", "at", "axis", "ref"],
        &[],
    ),
    b(
        "cone",
        "cone",
        &["r1", "r2", "h"],
        &["label", "at", "axis", "ref"],
        &[],
    ),
    b(
        "sphere",
        "sphere",
        &["r"],
        &["label", "at", "axis", "ref", "lat"],
        &[],
    ),
    b(
        "torus",
        "torus",
        &["major", "minor"],
        &["label", "at", "axis", "ref", "angle"],
        &[],
    ),
];

pub use crate::conformance_unbound::RHAI_UNBOUND_PREVIEW;
