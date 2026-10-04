//! PK-437: the parse-time expression-depth limits are pinned so every
//! build profile agrees. Rhai's own defaults are lower in debug (32
//! global / 16 in function bodies) than release (64 / 32); left
//! unpinned, a computed-geometry script with deep arithmetic could
//! parse in the release wasm build of popkern live yet be rejected by a debug
//! `cargo test` gate — a portability trap where the test bar is
//! stricter than production. `build_engine` calls
//! `set_max_expr_depths(64, 32)`; these tests fail (in debug) if that
//! ever regresses.

use popkern_rhai::ScriptHost;

#[test]
fn expr_depth_limits_are_pinned_across_build_profiles() {
    let host = ScriptHost::new();
    let engine = host.engine();
    assert_eq!(
        engine.max_expr_depth(),
        64,
        "global expr depth must be pinned to the release default"
    );
    assert_eq!(
        engine.max_function_expr_depth(),
        32,
        "function-body expr depth must be pinned to the release default"
    );
}

#[test]
fn a_computed_geometry_script_parses_in_every_profile() {
    // The shape that tripped the debug limit (16) during PK-434: an
    // arithmetic chain with a nested function call, plus a nested-array
    // point literal, inside two loops inside a fn. It parses at the
    // pinned function-body depth of 32.
    const SCRIPT: &str = r#"
        fn inv(a) { tan(a) - a }
        fn roll(rb, r) { let x = rb / r; if x >= 1.0 { 0.0 } else { acos(x) } }
        fn toothed(rb, ra, rf, cx, cy) {
            let inv_p = inv(0.5);
            let half_p = 0.1;
            let pts = [];
            for k in 0..3 {
                for s in 0..=6 {
                    let r = rf + (ra - rf) * s / 6;
                    let half = half_p + inv_p - inv(roll(rb, r)) + 0.004;
                    pts.push([cx + r * cos(k + half), cy + r * sin(k + half)]);
                }
            }
            pts
        }
        let m = model();
        toothed(0.2, 0.3, 0.18, 0.5, 0.5);
        m
    "#;
    let host = ScriptHost::new();
    host.compile(SCRIPT, "computed_geometry.rhai")
        .expect("a computed-geometry script must parse under the pinned depth");
}
