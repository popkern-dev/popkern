//! Boss-base fillet loops (PK-740): every root edge of a rectangular
//! rib on a plate filleted at once. Adjacent root strips share the
//! plate's REFLEX corner while the rib's vertical corner edge stays
//! sharp; the strips run past the corner and meet on their mutual
//! intersection seam (a mitre), from the plate rails' mitre point up to
//! the corner edge at the fillet radius.
//!
//! Volume is analytic: each strip adds its section `(1 − π/4) r²` along
//! its edge, and at each corner both strips run on into the corner
//! quadrant up to the bisector plane — each adding its section times
//! the section's centroid distance from its own wall,
//! `r (10 − 3π) / (3 (4 − π))`. `body_volume` meshes the strips, so
//! the gate is 5e-3 of the fillet's own volume; the OCC oracle measured
//! both radii exact on STEP (r = 1: 63219.27, r = 2: 63278.61).

use popkern_rhai::ScriptHost;

fn served_volume(source: &str) -> f64 {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(source, "boss.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
    let mut model = outcome.model.expect("model");
    let part = model.to_part().expect("part");
    let brep = part.brep();
    assert!(brep.validate().is_ok(), "{:?}", brep.validate());
    popkern_measure::body_volume(brep, &brep.tolerance).expect("volume")
}

fn rib_loop(r: f64, fillet: bool) -> String {
    let fillet_line = if fillet {
        format!(
            r#"let a = edges().at_vertex([20.0, 2.0, 10.0]).direction("x").or(edges().at_vertex([20.0, 2.0, 10.0]).direction("y"));
let b = edges().at_vertex([-20.0, -2.0, 10.0]).direction("x").or(edges().at_vertex([-20.0, -2.0, 10.0]).direction("y"));
p = p.fillet(a.or(b), {r:?}, #{{ label: "root" }});"#
        )
    } else {
        String::new()
    };
    format!(
        r#"
let m = model();
let p = m.profile(xy()).rect(100.0, 60.0).sketch("plate").extrude(blind(10.0));
let rib = m.tool(|t| t.profile(xy_at(10.0)).rect(40.0, 4.0).sketch("rib").extrude(blind(20.0)));
p = p.union(rib, #{{ label: "rib" }});
{fillet_line}
m.assert_watertight();
m
"#
    )
}

#[test]
fn a_rib_base_loop_fillets_with_mitred_corners() {
    let base = served_volume(&rib_loop(2.0, false));
    assert!((base - (100.0 * 60.0 * 10.0 + 40.0 * 4.0 * 20.0)).abs() < 1e-6);
    for r in [1.0, 2.0] {
        let pi = std::f64::consts::PI;
        let a = (1.0 - pi / 4.0) * r * r;
        let centroid = r * (10.0 - 3.0 * pi) / (3.0 * (4.0 - pi));
        let want = a * (2.0 * 40.0 + 2.0 * 4.0) + 8.0 * a * centroid;
        let change = served_volume(&rib_loop(r, true)) - base;
        assert!(
            (change - want).abs() / want < 5e-3,
            "r = {r}: fillet volume {change} vs analytic {want}"
        );
    }
}
