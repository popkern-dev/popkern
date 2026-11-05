//! brep-bench boolean/offset-hole-sphere: a ball of radius 20 bored by
//! a radius-6 hole along y at x = 5. The marching intersection traced
//! one loop distorted and the other as a point, and the subtraction
//! served the ball unchanged (volume 33510 for 29241). The exact
//! ruling-scan pierce serves the bored ball: its area to the bench's
//! OCCT-agreed value (its volume reads 1.2e-5 high in the measure,
//! PK-1092).

use popkern_rhai::ScriptHost;

fn bored(ball: &str) -> (f64, f64) {
    let source = format!(
        "let m = model();\nlet b = {ball};\nlet t = m.tool(|m| m.profile(xz()).circle(6.0).at(5.0, 0.0).sketch(\"bore\").extrude(mid_plane(21.0)));\nb = b.subtract(t);\nm\n"
    );
    let mut host = ScriptHost::with_options(Default::default());
    let mut model = host
        .run(&source, "bored")
        .expect("serves")
        .model
        .expect("model");
    let parts = model.to_parts_all().expect("parts");
    let brep = parts[0].brep();
    (
        popkern_measure::body_volume(brep, &brep.tolerance).expect("volume"),
        popkern_measure::body_surface_area(brep, &brep.tolerance).expect("area"),
    )
}

const AREA: f64 = 6171.861947620873;
const VOLUME: f64 = 29240.844399008805;

#[test]
fn a_ball_bored_off_centre_keeps_its_hole() {
    let (v, a) = bored("m.sphere(20.0)");
    assert!((a - AREA).abs() < 1e-9 * AREA, "area {a}");
    assert!((v - VOLUME).abs() < 2e-5 * VOLUME, "volume {v}");
}
