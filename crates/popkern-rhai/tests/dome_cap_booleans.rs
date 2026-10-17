//! A revolved dome — a sphere cap holding the chart's top pole, bounded
//! by its rim and the two uses of its seam, with no pole edge — against
//! a rod that crosses the dome's untrimmed sphere well off the cap
//! (PK-847: the six-axis arm's base column and connector boss against
//! the j1 housing's dome). The trim containment read every point on the
//! wrong side of the rim, kept the off-cap records, and the subtract
//! declined (4 stranded sides) while the union failed validation; the
//! base/housing contact declined as "interfering" in 62 s.
//!
//! The bodies are disjoint: the union is both, the subtract the dome,
//! the intersection empty, against the analytic volumes.

use popkern_rhai::ScriptHost;

const DOME: &str = r#"
let m = model(#{ unit: "mm" });
let ro = 55.0; let len = 118.0; let dome = 9.9;
let zc = (ro * ro - dome * dome) / (2.0 * dome);
let rho = zc + dome;
let vr = 0.7 * ro;
let vz = zc - sqrt(rho * rho - vr * vr);
let p = m.profile(xz()).path(0.0, -dome).arc_to([vr, vz], [ro, 0.0]);
p = p.line_to(ro, len).line_to(0.0, len).close();
let body = p.sketch("shell").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
let rod = m.tool(|t| t.profile(plane_at([70.0, 0.0, 20.0], [1.0, 0.0, 0.0])).circle(14.0)
    .sketch("rod").extrude(blind(30.0)));
body.OP(rod);
m
"#;

/// The body's volume in mm^3, or `None` when it is empty.
fn volume(op: &str) -> Option<f64> {
    let mut run = ScriptHost::new()
        .run(&DOME.replace("OP", op), "dome.rhai")
        .unwrap_or_else(|e| panic!("{op}: {}", e.render()));
    let parts = run.model.take().expect("model").to_parts_all().ok()?;
    let brep = parts.last().expect("a part").brep();
    let issues = brep.validate_all();
    assert!(issues.is_empty(), "{op}: {issues:?}");
    let v = popkern_measure::body_volume(brep, &brep.tolerance).expect("volume");
    Some(v * 1e9)
}

#[test]
fn a_dome_and_a_rod_off_its_cap_combine_as_disjoint_bodies() {
    use std::f64::consts::PI;
    let (ro, len, d) = (55.0f64, 118.0f64, 9.9f64);
    let rho = (ro * ro - d * d) / (2.0 * d) + d;
    let dome = PI * d * d * (3.0 * rho - d) / 3.0 + PI * ro * ro * len;
    let rod = PI * 14.0 * 14.0 * 30.0;
    let close = |got: f64, want: f64| ((got - want) / want).abs() < 1e-6;
    let union = volume("union").expect("a union");
    assert!(close(union, dome + rod), "union {union} vs {}", dome + rod);
    let subtract = volume("subtract").expect("a subtract");
    assert!(close(subtract, dome), "subtract {subtract} vs {dome}");
    assert_eq!(volume("intersect"), None, "the intersection is empty");
}
