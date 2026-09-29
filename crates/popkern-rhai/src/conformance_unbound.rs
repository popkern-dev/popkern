//! The preview ops the Rhai skin leaves unbound, each with the reason
//! (split from `conformance_bindings` for the file budget).

/// Preview ops the Rhai skin deliberately does not bind: `(op,
/// section or "" for any, reason)`.
pub const RHAI_UNBOUND_PREVIEW: &[(&str, &str, &str)] = &[
    (
        "sketch",
        "3.2",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "circle",
        "3.2",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "config_api",
        "4.3",
        "not an op: the history module the section mirrors",
    ),
    (
        "replay_lenient",
        "",
        "requires continue-on-failure replay in popkern-history",
    ),
    (
        "to_stl",
        "",
        "mesh export lands with the CLI (phase 3e) where a viewer needs it",
    ),
    (
        "to_obj",
        "",
        "mesh export lands with the CLI (phase 3e) where a viewer needs it",
    ),
    (
        "point",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "line",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "arc",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "coincident",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "horizontal",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "vertical",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "parallel",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "perpendicular",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "tangent",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "equal",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "fixed",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "distance",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "angle",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "radius",
        "",
        "constraint-sketch builder is not in the Rust skin yet (spec §3.2)",
    ),
    (
        "thicken",
        "",
        "offset-family features are not in the Rust skin yet (spec §5 preview)",
    ),
    (
        "offset_solid",
        "",
        "offset-family features are not in the Rust skin yet (spec §5 preview)",
    ),
    (
        "draft",
        "",
        "draft is not in the Rust skin yet (spec §5 preview)",
    ),
    (
        "transform",
        "",
        "general 4x4 transforms are not in the Rust skin yet (spec §5 preview)",
    ),
    (
        "linear_pattern",
        "",
        "patterns are not in the Rust skin yet (spec §5 preview)",
    ),
    (
        "circular_pattern",
        "",
        "patterns are not in the Rust skin yet (spec §5 preview)",
    ),
    (
        "curve_pattern",
        "",
        "patterns are not in the Rust skin yet (spec §5 preview)",
    ),
    (
        "mirror",
        "",
        "mirror is not in the Rust skin yet (spec §5 preview)",
    ),
    (
        "raw_body",
        "",
        "a BRep has no script value; snapshots enter through the Rust skin",
    ),
    (
        "reference",
        "",
        "measurement references are not in the Rust skin yet (spec §5 preview)",
    ),
    (
        "measure",
        "",
        "measurement queries are not in the Rust skin yet (spec §8 preview)",
    ),
    (
        "max_fillet_radius",
        "",
        "slider bounds are not in the Rust skin yet (spec §8 preview)",
    ),
    (
        "precheck_fillet",
        "",
        "feasibility precheck is not in the Rust skin yet (spec §8 preview)",
    ),
    (
        "dependencies",
        "",
        "dependency query is not in the Rust skin yet (spec §7 preview)",
    ),
];
