//! Engine construction: a raw Rhai engine with only pure packages
//! (design invariant 9 — no clock, no RNG, no environment, no file
//! access, no `eval`), the popkern dialect registered on top.

use rhai::Engine;
use rhai::packages::{
    ArithmeticPackage, BasicArrayPackage, BasicFnPackage, BasicIteratorPackage, BasicMapPackage,
    BasicMathPackage, BasicStringPackage, BitFieldPackage, LogicPackage, MoreStringPackage,
    Package,
};

/// Names a script may never call: the standard-package functions that
/// reach a clock or the host (`timestamp`, `sleep`, `exit`) and the
/// language's dynamic-evaluation escape (`eval`, disabled at parse).
pub const FORBIDDEN_FUNCTIONS: &[&str] = &["timestamp", "sleep", "exit", "eval"];

/// Build the dialect engine. `print` receives every `print`/`debug`
/// line.
pub fn build_engine(print: impl Fn(&str) + 'static) -> Engine {
    let mut engine = Engine::new_raw();
    // Pin the parse-time expression-depth limits explicitly so every
    // build profile agrees. Rhai's defaults are LOWER in debug (32
    // global / 16 in function bodies) than release (64 / 32), so a
    // computed-geometry script — `a + b*c/d` chains and nested-array
    // point literals, which naturally nest deep — could parse in the
    // release wasm build of popkern live yet be rejected by a debug `cargo test`
    // gate, a portability trap where the test bar is stricter than
    // production. These are Rhai's release defaults: a wasm-safe
    // ceiling (untrusted popkern live already parses at 64/32) that
    // still admits the arithmetic depth such scripts want. Parse-depth
    // is a stack-overflow guard, not a runtime lever — runtime is
    // bounded separately by `set_max_operations`.
    engine.set_max_expr_depths(64, 32);
    // `CorePackage` minus `LanguageCorePackage` (which carries `sleep`,
    // `exit`, `take`, tags, `parse_json`, and metadata introspection).
    engine.register_global_module(ArithmeticPackage::new().as_shared_module());
    engine.register_global_module(BasicStringPackage::new().as_shared_module());
    engine.register_global_module(BasicIteratorPackage::new().as_shared_module());
    engine.register_global_module(BasicFnPackage::new().as_shared_module());
    engine.register_global_module(LogicPackage::new().as_shared_module());
    engine.register_global_module(BasicMathPackage::new().as_shared_module());
    engine.register_global_module(BitFieldPackage::new().as_shared_module());
    engine.register_global_module(BasicArrayPackage::new().as_shared_module());
    engine.register_global_module(BasicMapPackage::new().as_shared_module());
    engine.register_global_module(MoreStringPackage::new().as_shared_module());
    crate::portable_math::register(&mut engine);
    // `eval` is a keyword-level built-in, not a package function.
    engine.disable_symbol("eval");
    let print = std::rc::Rc::new(print);
    let sink = std::rc::Rc::clone(&print);
    engine.on_print(move |line| sink(line));
    engine.on_debug(move |line, _source, _pos| print(line));
    engine.set_fail_on_invalid_map_property(true);
    crate::selection_api::register(&mut engine);
    crate::profile_api::register(&mut engine);
    crate::recipe_api::register(&mut engine);
    crate::model_api::register(&mut engine);
    crate::project_host::register_outside_project(&mut engine);
    crate::body_api::register(&mut engine);
    crate::fastener_api::register(&mut engine);
    crate::feature_api::register(&mut engine);
    crate::assembly_api::register(&mut engine);
    crate::assembly_mate_api::register(&mut engine);
    crate::assembly_param_api::register(&mut engine);
    engine
}

/// Every function name the engine exposes (packages included), for
/// the purity gate.
pub fn registered_function_names(engine: &Engine) -> Vec<String> {
    let mut names: Vec<String> = engine
        .gen_fn_signatures(true)
        .into_iter()
        .filter_map(|signature| {
            signature
                .split('(')
                .next()
                .map(|head| head.trim().rsplit(' ').next().unwrap_or(head).to_string())
        })
        .collect();
    names.sort();
    names.dedup();
    names
}
