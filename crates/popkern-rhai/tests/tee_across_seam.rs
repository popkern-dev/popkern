//! PK-752: the six-axis arm's root castings. A disc about the y axis
//! and a boss along x whose intersection crosses the disc's seam, in
//! both extrude orientations (the seam toward the boss, or away). The
//! tee serves with OpenCASCADE's own fused volume; the casting's foot
//! bolts then cut cleanly (their crossing of the notched foot face was
//! never searched when the face's rim bulged past its vertices, so a
//! pocket wall below the foot read Inside).

mod union_oracle;

const PRELUDE: &str = r#"
fn rod_at(t, p, n, r, length, label) {
    t.profile(plane_at(p, n)).circle(r).sketch(label).extrude(blind(length))
}
"#;

fn disc(rf: f64, toward_minus_y: bool) -> String {
    if toward_minus_y {
        format!("rod_at(m, [0.0, -1.0, 0.0], [0.0, -1.0, 0.0], {rf:?}, 12.0, \"foot\")")
    } else {
        format!("rod_at(m, [0.0, -13.0, 0.0], [0.0, 1.0, 0.0], {rf:?}, 12.0, \"foot\")")
    }
}

fn boss(rb: f64, zb: f64) -> String {
    format!("rod_at(m, [-90.0, -45.0, {zb:?}], [1.0, 0.0, 0.0], {rb:?}, 90.0, \"boss\")")
}

#[test]
fn the_tee_serves_with_the_fused_volume_whichever_way_the_disc_is_extruded() {
    let mut failures = Vec::new();
    let cases = [
        (50.0, 42.5, 0.0),
        (50.0, 33.0, 3.0),
        (55.0, 42.5, 3.0),
        (55.0, 38.0, 0.0),
    ];
    for (rf, rb, zb) in cases {
        for toward_minus_y in [true, false] {
            let name = format!("tee_{rf}_{rb}_{zb}_{toward_minus_y}");
            let (d, b) = (disc(rf, toward_minus_y), boss(rb, zb));
            let body = format!(
                "let f = {d};\nf.union(m.tool(|t| {}), #{{ label: \"tee\" }});\nm.assert_watertight();",
                b.replacen("(m,", "(t,", 1)
            );
            let mut model = match union_oracle::run(PRELUDE, &body, &format!("{name}.rhai")) {
                Ok(m) => m,
                Err(e) => {
                    failures.push(e);
                    continue;
                }
            };
            if let Some((served, fused)) =
                union_oracle::occ_volumes(PRELUDE, &[d.as_str(), b.as_str()], &mut model, &name)
            {
                let rel = (served - fused).abs() / fused;
                if rel > 1e-9 {
                    failures.push(format!(
                        "{name}: volume {served} vs OCC fuse {fused} (rel {rel:.2e})"
                    ));
                }
            }
        }
    }
    assert!(failures.is_empty(), "{failures:#?}");
}

/// The upper-arm root casting as built in the arm: the natural foot
/// (seam toward the boss), the boss, then a foot counterbore beside the
/// boss notch — watertight, with OCC's own cut volume.
#[test]
fn the_arm_root_foot_bolt_beside_the_boss_notch_cuts_cleanly() {
    let a = 210.0_f64.to_radians();
    let (px, pz) = (42.0 * a.cos(), 162.5 + 42.0 * a.sin());
    let tee = "rod_at(m, [0.0, -1.0, 162.5], [0.0, -1.0, 0.0], 60.0, 12.0, \"foot\").union(m.tool(|t| rod_at(t, [-90.0, -45.0, 162.5], [1.0, 0.0, 0.0], 42.5, 42.0, \"boss\")))".to_string();
    let bolt = format!(
        "rod_at(m, [{px:?}, -18.0, {pz:?}], [0.0, 1.0, 0.0], 5.0, 11.0, \"pocket\").union(m.tool(|t| rod_at(t, [{px:?}, -7.5, {pz:?}], [0.0, 1.0, 0.0], 2.75, 7.0, \"clear\")))"
    );
    let body = format!(
        "let b = {tee};\nb = b.subtract(m.tool(|t| {}), #{{ label: \"bolt\" }});\nm.assert_watertight();",
        bolt.replacen("rod_at(m,", "rod_at(t,", 1)
            .replacen("m.tool(|t|", "t.tool(|t|", 1)
    );
    let mut model =
        union_oracle::run(PRELUDE, &body, "root_bolt.rhai").unwrap_or_else(|e| panic!("{e}"));
    if let Some((served, cut)) = union_oracle::occ_volumes_op(
        PRELUDE,
        &[tee.as_str(), bolt.as_str()],
        &mut model,
        "root_bolt",
        "cut",
    ) {
        let rel = (served - cut).abs() / cut;
        assert!(
            rel < 1e-9,
            "volume {served} vs OCC cut {cut} (rel {rel:.2e})"
        );
    }
}
