//! A chamfer that stops at a tangent junction (PK-429): chamfering only
//! a rounded plate's four top corner arcs, whose ends meet the straight
//! top edges G1. Each arc's chain ends in a tangent cap, and the cap's
//! section used to be authored as a rolling-ball arc: on a chamfer the
//! corner cone meets its cap along a 45-degree LINE, the arc bowed
//! 9.2e-3 off the cone, and check 9 declined all four edges. The plate
//! now serves validate_all clean with the Pappus volume: each quarter
//! cone strip removes (pi/2) * 0.14 * 0.00045 (its section triangle's
//! area and centroid radius).

#[test]
fn chamfering_only_the_corner_arcs_serves_the_exact_volume() {
    let src = r#"let m = model();
let p = m.profile(xy()).rect(1.6, 1.0).sketch("s").extrude(blind(0.08));
p = p.fillet(edges().straight().direction("z"), 0.15);
p = p.chamfer(edges().circular().at_height(0.08), 0.03);
m"#;
    let mut run = popkern_rhai::ScriptHost::new()
        .run(src, "plate.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    let parts = run
        .model
        .take()
        .expect("model")
        .to_parts_all()
        .expect("parts");
    let b = parts[0].brep();
    let errors = b.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    let pi = std::f64::consts::PI;
    let exact = 0.08 * (1.6 * 1.0 - (4.0 - pi) * 0.0225) - 4.0 * (pi / 2.0) * 0.14 * 0.00045;
    let v = popkern_measure::body_volume(b, &b.tolerance).expect("volume");
    assert!((v - exact).abs() < 1e-9, "volume {v} vs Pappus {exact}");
}

/// The whole top loop chamfered (arcs and straight edges, one closed
/// tangent chain): each arc-to-line junction's section is a line on
/// both strips (the cone and the plane), authored as an arc before.
#[test]
fn chamfering_the_whole_tangent_top_loop_serves_the_exact_volume() {
    let src = r#"let m = model();
let p = m.profile(xy()).rect(1.6, 1.0).sketch("s").extrude(blind(0.08));
p = p.fillet(edges().straight().direction("z"), 0.15);
p = p.chamfer(edges().at_height(0.08), 0.03);
m"#;
    let mut run = popkern_rhai::ScriptHost::new()
        .run(src, "plate.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    let parts = run
        .model
        .take()
        .expect("model")
        .to_parts_all()
        .expect("parts");
    let b = parts[0].brep();
    let errors = b.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
    let pi = std::f64::consts::PI;
    let straight = 2.0 * (1.6 - 0.3) + 2.0 * (1.0 - 0.3);
    let exact =
        0.08 * (1.6 - (4.0 - pi) * 0.0225) - straight * 0.00045 - 4.0 * (pi / 2.0) * 0.14 * 0.00045;
    let v = popkern_measure::body_volume(b, &b.tolerance).expect("volume");
    assert!((v - exact).abs() < 1e-9, "volume {v} vs analytic {exact}");
}
