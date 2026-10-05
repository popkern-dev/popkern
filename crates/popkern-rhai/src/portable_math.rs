//! Script math with the same bits on every target (PK-674).
//!
//! Rhai's `BasicMathPackage` and its built-in `**` call the platform
//! math library, which differs in the last ulp between a native build
//! and the wasm page — enough for one script to build two different
//! parts (the propeller's blade twist). Every transcendental a script
//! can reach is re-registered here over [`popkern_math::PortableMath`],
//! in the engine's own namespace, which Rhai searches before packages.
//! `**` is a built-in operator; registering it takes effect only with
//! fast operators off, which [`register`] sets.

use popkern_math::PortableMath;
use rhai::{Engine, FLOAT, INT};

/// A one-argument script math function.
type Unary = fn(FLOAT) -> FLOAT;

/// Override the script math functions and the float `**` operators.
pub fn register(engine: &mut Engine) {
    engine.set_fast_operators(false);
    let unary: [(&str, Unary); 15] = [
        ("sin", FLOAT::psin),
        ("cos", FLOAT::pcos),
        ("tan", FLOAT::ptan),
        ("sinh", FLOAT::psinh),
        ("cosh", FLOAT::pcosh),
        ("tanh", FLOAT::ptanh),
        ("asin", FLOAT::pasin),
        ("acos", FLOAT::pacos),
        ("atan", FLOAT::patan),
        ("asinh", FLOAT::pasinh),
        ("acosh", FLOAT::pacosh),
        ("atanh", FLOAT::patanh),
        ("exp", FLOAT::pexp),
        ("ln", FLOAT::pln),
        ("log", FLOAT::plog10),
    ];
    for (name, f) in unary {
        engine.register_fn(name, f);
    }
    engine.register_fn("atan", |y: FLOAT, x: FLOAT| y.patan2(x));
    engine.register_fn("hypot", |x: FLOAT, y: FLOAT| x.phypot(y));
    engine.register_fn("log", |x: FLOAT, base: FLOAT| x.pln() / base.pln());
    engine.register_fn("**", |x: FLOAT, y: FLOAT| x.ppowf(y));
    engine.register_fn("**", |x: FLOAT, y: INT| x.ppowf(y as FLOAT));
    engine.register_fn("**", |x: INT, y: FLOAT| (x as FLOAT).ppowf(y));
    engine.register_fn("**=", |x: &mut FLOAT, y: FLOAT| *x = x.ppowf(y));
    engine.register_fn("**=", |x: &mut FLOAT, y: INT| *x = x.ppowf(y as FLOAT));
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The overrides are what scripts reach: results equal the portable
    /// functions bit for bit, including through `**` and `**=`.
    #[test]
    fn scripts_reach_the_portable_functions() {
        let engine = crate::engine::build_engine(|_| {});
        let eval = |s: &str| engine.eval::<FLOAT>(s).map(f64::to_bits).ok();
        let x = 0.7310585786300049f64;
        assert_eq!(eval("sin(0.7310585786300049)"), Some(x.sin().to_bits()));
        assert_eq!(eval("atan(0.3, -0.9)"), Some(0.3f64.atan2(-0.9).to_bits()));
        assert_eq!(eval("log(100.0)"), Some(100f64.log10().to_bits()));
        assert_eq!(eval("1.37 ** 0.61"), Some(1.37f64.powf(0.61).to_bits()));
        assert_eq!(eval("1.37 ** 3"), Some(1.37f64.powf(3.0).to_bits()));
        assert_eq!(eval("2 ** 0.5"), Some(2f64.powf(0.5).to_bits()));
        assert_eq!(
            eval("let a = 1.37; a **= 0.61; a"),
            Some(1.37f64.powf(0.61).to_bits())
        );
        // Integer powers stay integers.
        assert_eq!(engine.eval::<INT>("2 ** 10").ok(), Some(1024));
    }
}
