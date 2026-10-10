//! PK-753: the six-axis arm's wrist casting. A web box joins two lumps,
//! a foot disc (axis −y) and a ring with a low web, touching both. The
//! foot's wall carries the ring's hole from the earlier union. Ray
//! classification counted hits inside a curved face's holes, so web
//! fragments on the far side of the hole read Inside and the union
//! declined (the face-split "3.138 away" decline was its downstream
//! symptom). The union serves watertight, and its volume matches
//! OpenCASCADE's own fuse of the four operands.

mod union_oracle;

const PRELUDE: &str = r#"
fn rod_at(t, p, n, r, length, label) {
    t.profile(plane_at(p, n)).circle(r).sketch(label).extrude(blind(length))
}
fn box(t, x0, x1, y0, y1, z0, z1, label) {
    t.profile(xy_at(z0)).rect(x1 - x0, y1 - y0).at(0.5 * (x0 + x1), 0.5 * (y0 + y1))
        .sketch(label).extrude(blind(z1 - z0))
}
"#;

const FOOT: &str = r#"rod_at(m, [-817.0, -27.0, 162.5], [0.0, -1.0, 0.0], 35.25, 12.0, "foot")"#;
const RING: &str = r#"rod_at(m, [-817.0, -133.3, 133.0], [0.0, 0.0, -1.0], 49.0, 10.0, "r")"#;
const LOW: &str = r#"box(m, -829.0, -805.0, -92.0, -40.0, 124.0, 129.0, "wl")"#;
const WEB: &str = r#"box(m, -831.0, -803.0, -86.0, -30.0, 126.0, 150.0, "w")"#;

#[test]
fn a_web_joining_two_lumps_serves_with_the_fused_volume() {
    let tool = |e: &str| e.replacen("(m,", "(t,", 1);
    let body = format!(
        "let b = {FOOT};\n\
         b = b.union(m.tool(|t| {}), #{{ label: \"ring\" }});\n\
         b = b.union(m.tool(|t| {}), #{{ label: \"low\" }});\n\
         b = b.union(m.tool(|t| {}), #{{ label: \"web\" }});\n\
         m.assert_watertight();",
        tool(RING),
        tool(LOW),
        tool(WEB)
    );
    let mut model =
        union_oracle::run(PRELUDE, &body, "web_join.rhai").unwrap_or_else(|e| panic!("{e}"));
    let report = popkern_rhai::report_value(&model);
    assert_eq!(report["failures"], 0, "{report}");
    if let Some((served, fused)) =
        union_oracle::occ_volumes(PRELUDE, &[FOOT, RING, LOW, WEB], &mut model, "web_join")
    {
        let rel = (served - fused).abs() / fused;
        assert!(
            rel < 1e-9,
            "volume {served} vs OCC fuse {fused} (rel {rel:.2e})"
        );
    }
}
