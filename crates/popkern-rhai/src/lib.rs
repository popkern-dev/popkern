//! `popkern-rhai` — the embedded Rhai front-end of the scripting
//! surface (design D5/D7, phase 3).
//!
//! A script is a pure function of its source and parameters: the
//! engine carries only pure packages (no clock, RNG, environment,
//! file access, or `eval`), and every operation it exposes is a
//! spelling from `docs/architecture/scripting-surface.md` wrapped
//! over `popkern-script`. Errors cross as maps with the spec §10
//! payload (`kind`, `message`, fields) and the script position of the
//! raising call; `try`/`catch` receives the same map.
//!
//! ```rhai
//! let m = model(#{ unit: "m" });
//! let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate_sketch")
//!              .extrude(blind(0.4));
//! let boss = m.tool(|m| m.profile(xy_at(0.2)).circle(0.35).sketch("boss_sketch")
//!              .extrude(blind(0.8)));
//! plate = plate.union(boss);
//! plate.fillet(
//!     edges().between("cylinder", "plane").at_height(0.4)
//!         .or(edges().straight().at_height(0.4))
//!         .or(edges().straight().direction("z")),
//!     0.08, #{ continuity: "g1", label: "rounds" });
//! m.assert_watertight();
//! m.save("boss_plate.pkn");
//! ```
//!
//! Enum arguments are strings (`"g2"`, `"cylinder"`, `"inward"`);
//! optional arguments ride in a trailing `#{ .. }` map whose unknown
//! keys are rejected; numeric arguments accept literals, `p("name")`
//! parameter references, and `expr("2 * wall")` expressions.

#![forbid(unsafe_code)]

pub mod args;
pub mod body_api;
pub mod conformance;
pub mod conformance_audit;
pub mod conformance_bindings;
pub mod conformance_rust_bindings;
pub mod emit;
pub mod emit_feature;
pub mod emit_spelling;
pub mod engine;
pub mod error;
pub mod feature_api;
pub mod host;
pub mod migration_api;
pub mod model_api;
pub mod outcome;
pub mod profile_api;
pub mod recipe_api;
pub mod selection_api;
pub mod values;

pub use conformance::{Checklist, checklist_from_spec};
pub use conformance_audit::{Finding, audit_error_classes, audit_rhai, audit_rust};
pub use emit::{EmitError, to_rhai};
pub use engine::{FORBIDDEN_FUNCTIONS, build_engine, registered_function_names};
pub use error::{error_map, payload, remediation_lines};
pub use host::{HostOptions, RunOutcome, ScriptHost, ScriptRunError};
pub use model_api::ScriptModel;
pub use outcome::{Outcome, classify, report_value};
